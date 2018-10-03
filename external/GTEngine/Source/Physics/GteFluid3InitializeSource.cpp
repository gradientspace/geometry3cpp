// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#include <GTEnginePCH.h>
#include <Physics/GteFluid3InitializeSource.h>
#include <Graphics/GteGraphicsEngine.h>
#include <Graphics/GteProgramFactory.h>
#include <random>
using namespace gte;

Fluid3InitializeSource::Fluid3InitializeSource(std::shared_ptr<ProgramFactory> const& factory,
    int xSize, int ySize, int zSize, int numXThreads, int numYThreads,
    int numZThreads, std::shared_ptr<ConstantBuffer> const& parameters)
    :
    mNumXGroups(xSize/numXThreads),
    mNumYGroups(ySize/numYThreads),
    mNumZGroups(zSize/numZThreads)
{
    // Create the resources for generating velocity from vortices.
    mVortex = std::make_shared<ConstantBuffer>(sizeof(Vortex), true);
    mVelocity0 = std::make_shared<Texture3>(DF_R32G32B32A32_FLOAT, xSize, ySize, zSize);
    mVelocity0->SetUsage(Resource::SHADER_OUTPUT);
    mVelocity1 = std::make_shared<Texture3>(DF_R32G32B32A32_FLOAT, xSize, ySize, zSize);
    mVelocity1->SetUsage(Resource::SHADER_OUTPUT);

    // Create the resources for generating velocity from wind and gravity.
    mExternal = std::make_shared<ConstantBuffer>(sizeof(External), false);
    External& e = *mExternal->Get<External>();
    e.densityProducer = { 0.5f, 0.5f, 0.5f, 0.0f };
    e.densityPData = { 0.01f, 16.0f, 0.0f, 0.0f };
    e.densityConsumer = { 0.75f, 0.75f, 0.75f, 0.0f };
    e.densityCData = { 0.01f, 0.0f, 0.0f, 0.0f };
    e.gravity = { 0.0f, 0.0f, 0.0f, 0.0f };
    e.windData = { 0.001f, 0.0f, 0.0f, 0.0f };
    mSource = std::make_shared<Texture3>(DF_R32G32B32A32_FLOAT, xSize, ySize, zSize);
    mSource->SetUsage(Resource::SHADER_OUTPUT);

    // Create the shader for generating velocity from vortices.
    int i = factory->GetAPI();
    factory->PushDefines();
    factory->defines.Set("NUM_X_THREADS", numXThreads);
    factory->defines.Set("NUM_Y_THREADS", numYThreads);
    factory->defines.Set("NUM_Z_THREADS", numZThreads);
    std::shared_ptr<ComputeShader> cshader;

    mGenerateVortex = factory->CreateFromSource(*msGenerateSource[i]);
    if (mGenerateVortex)
    {
        cshader = mGenerateVortex->GetCShader();
        cshader->Set("Parameters", parameters);
        cshader->Set("Vortex", mVortex);
        cshader->Set("inVelocity", mVelocity0);
        cshader->Set("outVelocity", mVelocity1);
    }

    // Create the shader for generating the sources to the fluid simulation.
    mInitializeSource = factory->CreateFromSource(*msInitializeSource[i]);
    if (mInitializeSource)
    {
        cshader = mInitializeSource->GetCShader();
        cshader->Set("Parameters", parameters);
        cshader->Set("External", mExternal);
        cshader->Set("source", mSource);
    }

    factory->PopDefines();
}

void Fluid3InitializeSource::Execute(std::shared_ptr<GraphicsEngine> const& engine)
{
    // Use a Mersenne twister engine for random numbers.
    std::mt19937 mte;
    std::uniform_real_distribution<float> unirnd(0.0f, 1.0f);
    std::uniform_real_distribution<float> symrnd(-1.0f, 1.0f);
    std::uniform_real_distribution<float> posrnd0(0.001f, 0.01f);
    std::uniform_real_distribution<float> posrnd1(64.0f, 128.0f);

    // Compute the velocity one vortex at a time.  After the loop terminates,
    // the final velocity is stored in mVelocity0.
    std::shared_ptr<ComputeShader> cshader = mGenerateVortex->GetCShader();
    memset(mVelocity0->GetData(), 0, mVelocity0->GetNumBytes());
    Vortex& v = *mVortex->Get<Vortex>();
    for (int i = 0; i < NUM_VORTICES; ++i)
    {
        v.position[0] = unirnd(mte);
        v.position[1] = unirnd(mte);
        v.position[2] = unirnd(mte);
        v.position[3] = 0.0f;
        v.normal[0] = symrnd(mte);
        v.normal[1] = symrnd(mte);
        v.normal[2] = symrnd(mte);
        v.normal[3] = 0.0f;
        Normalize(v.normal);
        v.data[0] = posrnd0(mte);
        v.data[1] = posrnd1(mte);
        v.data[2] = 0.0f;
        v.data[3] = 0.0f;
        engine->Update(mVortex);

        engine->Execute(mGenerateVortex, mNumXGroups, mNumYGroups, mNumZGroups);

        std::swap(mVelocity0, mVelocity1);
        cshader->Set("inVelocity", mVelocity0);
        cshader->Set("outVelocity", mVelocity1);
    }

    // Compute the sources for the fluid simulation.
    cshader = mInitializeSource->GetCShader();
    cshader->Set("vortexVelocity", mVelocity0);
    engine->Execute(mInitializeSource, mNumXGroups, mNumYGroups, mNumZGroups);
}


std::string const Fluid3InitializeSource::msGLSLGenerateSource =
"uniform Parameters\n"
"{\n"
"    vec4 spaceDelta;    // (dx, dy, dz, 0)\n"
"    vec4 halfDivDelta;  // (0.5/dx, 0.5/dy, 0.5/dz, 0)\n"
"    vec4 timeDelta;     // (dt/dx, dt/dy, dt/dz, dt)\n"
"    vec4 viscosityX;    // (velVX, velVX, velVX, denVX)\n"
"    vec4 viscosityY;    // (velVX, velVY, velVY, denVY)\n"
"    vec4 viscosityZ;    // (velVZ, velVZ, velVZ, denVZ)\n"
"    vec4 epsilon;       // (epsilonX, epsilonY, epsilonZ, epsilon0)\n"
"};\n"
"\n"
"uniform Vortex\n"
"{\n"
"    vec4 position;  // (px, py, pz, *)\n"
"    vec4 normal;    // (nx, ny, nz, *)\n"
"    vec4 data;      // (variance, amplitude, *, *)\n"
"};\n"
"\n"
"layout(rgba32f) uniform readonly image3D inVelocity;\n"
"layout(rgba32f) uniform writeonly image3D outVelocity;\n"
"\n"
"layout (local_size_x = NUM_X_THREADS, local_size_y = NUM_Y_THREADS, local_size_z = NUM_Z_THREADS) in;\n"
"void main()\n"
"{\n"
"    ivec3 c = ivec3(gl_GlobalInvocationID.xyz);\n"
"    vec3 location = spaceDelta.xyz*(c + 0.5f);\n"
"    vec3 diff = location - position.xyz;\n"
"    float arg = -dot(diff, diff) / data.x;\n"
"    float magnitude = data.y*exp(arg);\n"
"    vec4 vortexVelocity = vec4(magnitude*cross(normal.xyz, diff), 0.0f);\n"
"    imageStore(outVelocity, c, imageLoad(inVelocity, c) + vortexVelocity);\n"
"}\n";

std::string const Fluid3InitializeSource::msGLSLInitializeSource =
"uniform Parameters\n"
"{\n"
"    vec4 spaceDelta;    // (dx, dy, dz, 0)\n"
"    vec4 halfDivDelta;  // (0.5/dx, 0.5/dy, 0.5/dz, 0)\n"
"    vec4 timeDelta;     // (dt/dx, dt/dy, dt/dz, dt)\n"
"    vec4 viscosityX;    // (velVX, velVX, velVX, denVX)\n"
"    vec4 viscosityY;    // (velVX, velVY, velVY, denVY)\n"
"    vec4 viscosityZ;    // (velVZ, velVZ, velVZ, denVZ)\n"
"    vec4 epsilon;       // (epsilonX, epsilonY, epsilonZ, epsilon0)\n"
"};\n"
"\n"
"uniform External\n"
"{\n"
"    vec4 densityProducer;  // (x, y, z, *)\n"
"    vec4 densityPData;     // (variance, amplitude, *, *)\n"
"    vec4 densityConsumer;  // (x, y, z, *)\n"
"    vec4 densityCData;     // (variance, amplitude, *, *)\n"
"    vec4 gravity;          // (x, y, z, *)\n"
"    vec4 windData;         // (variance, amplitude, *, *)\n"
"};\n"
"\n"
"layout(rgba32f) uniform readonly image3D vortexVelocity;\n"
"layout(rgba32f) uniform writeonly image3D source;\n"
"\n"
"layout (local_size_x = NUM_X_THREADS, local_size_y = NUM_Y_THREADS, local_size_z = NUM_Z_THREADS) in;\n"
"void main()\n"
"{\n"
"    ivec3 c = ivec3(gl_GlobalInvocationID.xyz);\n"
"\n"
"    // Compute the location of the voxel (x,y,z) in normalized [0,1]^3.\n"
"    vec3 location = spaceDelta.xyz*(c + 0.5f);\n"
"\n"
"    // Compute an input to the fluid simulation consisting of a producer of\n"
"    // density and a consumer of density.\n"
"    vec3 diff = location - densityProducer.xyz;\n"
"    float arg = -dot(diff, diff) / densityPData.x;\n"
"    float density = densityPData.y*exp(arg);\n"
"    diff = location - densityConsumer.xyz;\n"
"    arg = -dot(diff, diff) / densityCData.x;\n"
"    density -= densityCData.y*exp(arg);\n"
"\n"
"    // Compute an input to the fluid simulation consisting of gravity,\n"
"    // a single wind source, and vortex impulses.\n"
"    float windArg = -dot(location.xz, location.xz) / windData.x;\n"
"    vec3 windVelocity = vec3(0.0f, windData.y*exp(windArg), 0.0f);\n"
"    vec3 velocity = gravity.xyz + windVelocity + imageLoad(vortexVelocity, c).xyz;\n"
"\n"
"    imageStore(source, c, vec4(velocity.xyz, density));\n"
"}\n";

std::string const Fluid3InitializeSource::msHLSLGenerateSource =
"cbuffer Parameters\n"
"{\n"
"    float4 spaceDelta;    // (dx, dy, dz, 0)\n"
"    float4 halfDivDelta;  // (0.5/dx, 0.5/dy, 0.5/dz, 0)\n"
"    float4 timeDelta;     // (dt/dx, dt/dy, dt/dz, dt)\n"
"    float4 viscosityX;    // (velVX, velVX, velVX, denVX)\n"
"    float4 viscosityY;    // (velVX, velVY, velVY, denVY)\n"
"    float4 viscosityZ;    // (velVZ, velVZ, velVZ, denVZ)\n"
"    float4 epsilon;       // (epsilonX, epsilonY, epsilonZ, epsilon0)\n"
"};\n"
"\n"
"cbuffer Vortex\n"
"{\n"
"    float4 position;  // (px, py, pz, *)\n"
"    float4 normal;    // (nx, ny, nz, *)\n"
"    float4 data;      // (variance, amplitude, *, *)\n"
"};\n"
"\n"
"Texture3D<float4> inVelocity;\n"
"RWTexture3D<float4> outVelocity;\n"
"\n"
"[numthreads(NUM_X_THREADS, NUM_Y_THREADS, NUM_Z_THREADS)]\n"
"void CSMain(uint3 c : SV_DispatchThreadID)\n"
"{\n"
"    float3 location = spaceDelta.xyz*(c + 0.5f);\n"
"    float3 diff = location - position.xyz;\n"
"    float arg = -dot(diff, diff) / data.x;\n"
"    float magnitude = data.y*exp(arg);\n"
"    float4 vortexVelocity = float4(magnitude*cross(normal.xyz, diff), 0.0f);\n"
"    outVelocity[c] = inVelocity[c] + vortexVelocity;\n"
"}\n";

std::string const Fluid3InitializeSource::msHLSLInitializeSource =
"cbuffer Parameters\n"
"{\n"
"    float4 spaceDelta;    // (dx, dy, dz, 0)\n"
"    float4 halfDivDelta;  // (0.5/dx, 0.5/dy, 0.5/dz, 0)\n"
"    float4 timeDelta;     // (dt/dx, dt/dy, dt/dz, dt)\n"
"    float4 viscosityX;    // (velVX, velVX, velVX, denVX)\n"
"    float4 viscosityY;    // (velVX, velVY, velVY, denVY)\n"
"    float4 viscosityZ;    // (velVZ, velVZ, velVZ, denVZ)\n"
"    float4 epsilon;       // (epsilonX, epsilonY, epsilonZ, epsilon0)\n"
"};\n"
"\n"
"cbuffer External\n"
"{\n"
"    float4 densityProducer;  // (x, y, z, *)\n"
"    float4 densityPData;     // (variance, amplitude, *, *)\n"
"    float4 densityConsumer;  // (x, y, z, *)\n"
"    float4 densityCData;     // (variance, amplitude, *, *)\n"
"    float4 gravity;          // (x, y, z, *)\n"
"    float4 windData;         // (variance, amplitude, *, *)\n"
"};\n"
"\n"
"Texture3D<float4> vortexVelocity;\n"
"RWTexture3D<float4> source;\n"
"\n"
"[numthreads(NUM_X_THREADS, NUM_Y_THREADS, NUM_Z_THREADS)]\n"
"void CSMain(uint3 c : SV_DispatchThreadID)\n"
"{\n"
"    // Compute the location of the voxel (x,y,z) in normalized [0,1]^3.\n"
"    float3 location = spaceDelta.xyz*(c + 0.5f);\n"
"\n"
"    // Compute an input to the fluid simulation consisting of a producer of\n"
"    // density and a consumer of density.\n"
"    float3 diff = location - densityProducer.xyz;\n"
"    float arg = -dot(diff, diff) / densityPData.x;\n"
"    float density = densityPData.y*exp(arg);\n"
"    diff = location - densityConsumer.xyz;\n"
"    arg = -dot(diff, diff) / densityCData.x;\n"
"    density -= densityCData.y*exp(arg);\n"
"\n"
"    // Compute an input to the fluid simulation consisting of gravity,\n"
"    // a single wind source, and vortex impulses.\n"
"    float windArg = -dot(location.xz, location.xz) / windData.x;\n"
"    float3 windVelocity = float3(0.0f, windData.y*exp(windArg), 0.0f);\n"
"    float3 velocity = gravity.xyz + windVelocity + vortexVelocity[c].xyz;\n"
"\n"
"    source[c] = float4(velocity.xyz, density);\n"
"}\n";

std::string const* Fluid3InitializeSource::msGenerateSource[] =
{
    &msGLSLGenerateSource,
    &msHLSLGenerateSource
};

std::string const* Fluid3InitializeSource::msInitializeSource[] =
{
    &msGLSLInitializeSource,
    &msHLSLInitializeSource
};
