// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#include <GTEnginePCH.h>
#include <Physics/GteFluid2InitializeSource.h>
#include <Graphics/GteGraphicsEngine.h>
#include <Graphics/GteProgramFactory.h>
#include <random>
using namespace gte;

Fluid2InitializeSource::Fluid2InitializeSource(std::shared_ptr<ProgramFactory> const& factory,
    int xSize, int ySize, int numXThreads, int numYThreads,
    std::shared_ptr<ConstantBuffer> const& parameters)
    :
    mNumXGroups(xSize/numXThreads),
    mNumYGroups(ySize/numYThreads)
{
    // Create the resources for generating velocity from vortices.
    mVortex = std::make_shared<ConstantBuffer>(sizeof(Vortex), true);
    mVelocity0 = std::make_shared<Texture2>(DF_R32G32_FLOAT, xSize, ySize);
    mVelocity0->SetUsage(Resource::SHADER_OUTPUT);
    mVelocity1 = std::make_shared<Texture2>(DF_R32G32_FLOAT, xSize, ySize);
    mVelocity1->SetUsage(Resource::SHADER_OUTPUT);

    // Create the resources for generating velocity from wind and gravity.
    mExternal = std::make_shared<ConstantBuffer>(sizeof(External), false);
    External& e = *mExternal->Get<External>();
    e.densityProducer = { 0.25f, 0.75f, 0.01f, 2.0f };
    e.densityConsumer = { 0.75f, 0.25f, 0.01f, 2.0f };
    e.gravity = { 0.0f, 0.0f, 0.0f, 0.0f };
    e.wind = { 0.0f, 0.5f, 0.001f, 32.0f };
    mSource = std::make_shared<Texture2>(DF_R32G32B32A32_FLOAT, xSize, ySize);
    mSource->SetUsage(Resource::SHADER_OUTPUT);

    // Create the shader for generating velocity from vortices.
    int i = factory->GetAPI();
    factory->PushDefines();
    factory->defines.Set("NUM_X_THREADS", numXThreads);
    factory->defines.Set("NUM_Y_THREADS", numYThreads);
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

void Fluid2InitializeSource::Execute(std::shared_ptr<GraphicsEngine> const& engine)
{
    // Use a Mersenne twister engine for random numbers.
    std::mt19937 mte;
    std::uniform_real_distribution<float> unirnd(0.0f, 1.0f);
    std::uniform_real_distribution<float> symrnd(-1.0f, 1.0f);
    std::uniform_real_distribution<float> posrnd0(0.001f, 0.01f);
    std::uniform_real_distribution<float> posrnd1(128.0f, 256.0f);

    // Compute the velocity one vortex at a time.  After the loop terminates,
    // the final velocity is stored in mVelocity0.
    std::shared_ptr<ComputeShader> cshader = mGenerateVortex->GetCShader();
    memset(mVelocity0->GetData(), 0, mVelocity0->GetNumBytes());
    Vortex& v = *mVortex->Get<Vortex>();
    for (int i = 0; i < NUM_VORTICES; ++i)
    {
        v.data[0] = unirnd(mte);
        v.data[1] = unirnd(mte);
        v.data[2] = posrnd0(mte);
        v.data[3] = posrnd1(mte);
        if (symrnd(mte) < 0.0f)
        {
            v.data[3] = -v.data[3];
        }
        engine->Update(mVortex);

        engine->Execute(mGenerateVortex, mNumXGroups, mNumYGroups, 1);

        std::swap(mVelocity0, mVelocity1);
        cshader->Set("inVelocity", mVelocity0);
        cshader->Set("outVelocity", mVelocity1);
    }

    // Compute the sources for the fluid simulation.
    cshader = mInitializeSource->GetCShader();
    cshader->Set("vortexVelocity", mVelocity0);
    engine->Execute(mInitializeSource, mNumXGroups, mNumYGroups, 1);
}


// TODO:  Write these shaders.
std::string const Fluid2InitializeSource::msGLSLGenerateSource =
"uniform Parameters\n"
"{\n"
"    vec4 spaceDelta;    // (dx, dy, 0, 0)\n"
"    vec4 halfDivDelta;  // (0.5/dx, 0.5/dy, 0, 0)\n"
"    vec4 timeDelta;     // (dt/dx, dt/dy, 0, dt)\n"
"    vec4 viscosityX;    // (velVX, velVX, 0, denVX)\n"
"    vec4 viscosityY;    // (velVX, velVY, 0, denVY)\n"
"    vec4 epsilon;       // (epsilonX, epsilonY, 0, epsilon0)\n"
"};\n"
"\n"
"uniform Vortex\n"
"{\n"
"    vec4 data;      // (x, y, variance, amplitude)\n"
"};\n"
"\n"
"layout(rg32f) uniform readonly image2D inVelocity;\n"
"layout(rg32f) uniform writeonly image2D outVelocity;\n"
"\n"
"layout (local_size_x = NUM_X_THREADS, local_size_y = NUM_Y_THREADS, local_size_z = 1) in;\n"
"void main()\n"
"{\n"
"    ivec2 c = ivec2(gl_GlobalInvocationID.xy);\n"
"    vec2 location = spaceDelta.xy*(c + 0.5f);\n"
"    vec2 diff = location - data.xy;\n"
"    float arg = -dot(diff, diff) / data.z;\n"
"    float magnitude = data.w*exp(arg);\n"
"    vec2 vortexVelocity = magnitude*vec2(diff.y, -diff.x);\n"
"    imageStore(outVelocity, c, vec4(imageLoad(inVelocity, c).xy + vortexVelocity, 0.0f, 0.0f));\n"
"}\n";

std::string const Fluid2InitializeSource::msGLSLInitializeSource =
"uniform Parameters\n"
"{\n"
"    vec4 spaceDelta;    // (dx, dy, 0, 0)\n"
"    vec4 halfDivDelta;  // (0.5/dx, 0.5/dy, 0, 0)\n"
"    vec4 timeDelta;     // (dt/dx, dt/dy, 0, dt)\n"
"    vec4 viscosityX;    // (velVX, velVX, 0, denVX)\n"
"    vec4 viscosityY;    // (velVX, velVY, 0, denVY)\n"
"    vec4 epsilon;       // (epsilonX, epsilonY, 0, epsilon0)\n"
"};\n"
"\n"
"uniform External\n"
"{\n"
"    vec4 densityProducer;  // (x, y, variance, amplitude)\n"
"    vec4 densityConsumer;  // (x, y, variance, amplitude)\n"
"    vec4 gravity;          // (x, y, *, *)\n"
"    vec4 wind;             // (x, y, variance, amplitude)\n"
"};\n"
"\n"
"layout(rg32f) uniform readonly image2D vortexVelocity;\n"
"layout(rgba32f) uniform writeonly image2D source;\n"
"\n"
"layout (local_size_x = NUM_X_THREADS, local_size_y = NUM_Y_THREADS, local_size_z = 1) in;\n"
"void main()\n"
"{\n"
"    ivec2 c = ivec2(gl_GlobalInvocationID.xy);\n"
"\n"
"    // Compute the location of the pixel (x,y) in normalized [0,1]^2.\n"
"    vec2 location = spaceDelta.xy*(c + 0.5f);\n"
"\n"
"    // Compute an input to the fluid simulation consisting of a producer of\n"
"    // density and a consumer of density.\n"
"    vec2 diff = location - densityProducer.xy;\n"
"    float arg = -dot(diff, diff) / densityProducer.z;\n"
"    float density = densityProducer.w*exp(arg);\n"
"    diff = location - densityConsumer.xy;\n"
"    arg = -dot(diff, diff) / densityConsumer.z;\n"
"    density -= densityConsumer.w*exp(arg);\n"
"\n"
"    // Compute an input to the fluid simulation consisting of gravity,\n"
"    // a single wind source, and vortex impulses.\n"
"    float windDiff = location.y - wind.y;\n"
"    float windArg = -windDiff*windDiff / wind.z;\n"
"    vec2 windVelocity = vec2(wind.w*exp(windArg), 0.0f);\n"
"    vec2 velocity = gravity.xy + windVelocity + imageLoad(vortexVelocity, c).xy;\n"
"\n"
"    imageStore(source, c, vec4(velocity, 0.0f, density));\n"
"}\n";

std::string const Fluid2InitializeSource::msHLSLGenerateSource =
"cbuffer Parameters\n"
"{\n"
"    float4 spaceDelta;    // (dx, dy, 0, 0)\n"
"    float4 halfDivDelta;  // (0.5/dx, 0.5/dy, 0, 0)\n"
"    float4 timeDelta;     // (dt/dx, dt/dy, 0, dt)\n"
"    float4 viscosityX;    // (velVX, velVX, 0, denVX)\n"
"    float4 viscosityY;    // (velVX, velVY, 0, denVY)\n"
"    float4 epsilon;       // (epsilonX, epsilonY, 0, epsilon0)\n"
"};\n"
"\n"
"cbuffer Vortex\n"
"{\n"
"    float4 data;  // (x, y, variance, amplitude)\n"
"};\n"
"\n"
"Texture2D<float2> inVelocity;\n"
"RWTexture2D<float2> outVelocity;\n"
"\n"
"[numthreads(NUM_X_THREADS, NUM_Y_THREADS, 1)]\n"
"void CSMain(uint3 c : SV_DispatchThreadID)\n"
"{\n"
"    float2 location = spaceDelta.xy*(c.xy + 0.5f);\n"
"    float2 diff = location - data.xy;\n"
"    float arg = -dot(diff, diff) / data.z;\n"
"    float magnitude = data.w*exp(arg);\n"
"    float2 vortexVelocity = magnitude*float2(diff.y, -diff.x);\n"
"    outVelocity[c.xy] = inVelocity[c.xy] + vortexVelocity;\n"
"}\n";

std::string const Fluid2InitializeSource::msHLSLInitializeSource =
"cbuffer Parameters\n"
"{\n"
"    float4 spaceDelta;    // (dx, dy, 0, 0)\n"
"    float4 halfDivDelta;  // (0.5/dx, 0.5/dy, 0, 0)\n"
"    float4 timeDelta;     // (dt/dx, dt/dy, 0, dt)\n"
"    float4 viscosityX;    // (velVX, velVX, 0, denVX)\n"
"    float4 viscosityY;    // (velVX, velVY, 0, denVY)\n"
"    float4 epsilon;       // (epsilonX, epsilonY, 0, epsilon0)\n"
"};\n"
"\n"
"cbuffer External\n"
"{\n"
"    float4 densityProducer;  // (x, y, variance, amplitude)\n"
"    float4 densityConsumer;  // (x, y, variance, amplitude)\n"
"    float4 gravity;          // (x, y, *, *)\n"
"    float4 wind;             // (x, y, variance, amplitude)\n"
"};\n"
"\n"
"Texture2D<float2> vortexVelocity;\n"
"RWTexture2D<float4> source;\n"
"\n"
"[numthreads(NUM_X_THREADS, NUM_Y_THREADS, 1)]\n"
"void CSMain(uint2 c : SV_DispatchThreadID)\n"
"{\n"
"    // Compute the location of the pixel (x,y) in normalized [0,1]^2.\n"
"    float2 location = spaceDelta.xy*(c + 0.5f);\n"
"\n"
"    // Compute an input to the fluid simulation consisting of a producer of\n"
"    // density and a consumer of density.\n"
"    float2 diff = location - densityProducer.xy;\n"
"    float arg = -dot(diff, diff) / densityProducer.z;\n"
"    float density = densityProducer.w*exp(arg);\n"
"    diff = location - densityConsumer.xy;\n"
"    arg = -dot(diff, diff) / densityConsumer.z;\n"
"    density -= densityConsumer.w*exp(arg);\n"
"\n"
"    // Compute an input to the fluid simulation consisting of gravity,\n"
"    // a single wind source, and vortex impulses.\n"
"    float windDiff = location.y - wind.y;\n"
"    float windArg = -windDiff*windDiff / wind.z;\n"
"    float2 windVelocity = float2(wind.w*exp(windArg), 0.0f);\n"
"    float2 velocity = gravity.xy + windVelocity + vortexVelocity[c];\n"
"\n"
"    source[c] = float4(velocity, 0.0f, density);\n"
"}\n";

std::string const* Fluid2InitializeSource::msGenerateSource[] =
{
    &msGLSLGenerateSource,
    &msHLSLGenerateSource
};

std::string const* Fluid2InitializeSource::msInitializeSource[] =
{
    &msGLSLInitializeSource,
    &msHLSLInitializeSource
};
