// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#include <GTEnginePCH.h>
#include <Physics/GteFluid3SolvePoisson.h>
#include <Graphics/GteGraphicsEngine.h>
#include <Graphics/GteProgramFactory.h>
using namespace gte;

Fluid3SolvePoisson::Fluid3SolvePoisson(std::shared_ptr<ProgramFactory> const& factory,
    int xSize, int ySize, int zSize, int numXThreads, int numYThreads, int numZThreads,
    std::shared_ptr<ConstantBuffer> const& parameters, int numIterations)
    :
    mNumXGroups(xSize/numXThreads),
    mNumYGroups(ySize/numYThreads),
    mNumZGroups(zSize/numZThreads),
    mNumIterations(numIterations)
{
    mPoisson0 = std::make_shared<Texture3>(DF_R32_FLOAT, xSize, ySize, zSize);
    mPoisson0->SetUsage(Resource::SHADER_OUTPUT);
    mPoisson1 = std::make_shared<Texture3>(DF_R32_FLOAT, xSize, ySize, zSize);
    mPoisson1->SetUsage(Resource::SHADER_OUTPUT);

    int i = factory->GetAPI();
    factory->PushDefines();
    factory->defines.Set("NUM_X_THREADS", numXThreads);
    factory->defines.Set("NUM_Y_THREADS", numYThreads);
    factory->defines.Set("NUM_Z_THREADS", numZThreads);

    // For zeroing mPoisson0 on the GPU.
    mZeroPoisson = factory->CreateFromSource(*msZeroSource[i]);
    if (mZeroPoisson)
    {
        mZeroPoisson->GetCShader()->Set("poisson", mPoisson0);
    }

    // Create the shader for generating velocity from vortices.
    mSolvePoisson = factory->CreateFromSource(*msSolveSource[i]);
    if (mSolvePoisson)
    {
        mSolvePoisson->GetCShader()->Set("Parameters", parameters);
    }

    factory->defines.Clear();
    factory->defines.Set("USE_ZERO_X_FACE", 1);
    factory->defines.Set("NUM_Y_THREADS", numYThreads);
    factory->defines.Set("NUM_Z_THREADS", numZThreads);
    mWriteXFace = factory->CreateFromSource(*msEnforceSource[i]);

    factory->defines.Clear();
    factory->defines.Set("USE_ZERO_Y_FACE", 1);
    factory->defines.Set("NUM_X_THREADS", numXThreads);
    factory->defines.Set("NUM_Z_THREADS", numZThreads);
    mWriteYFace = factory->CreateFromSource(*msEnforceSource[i]);

    factory->defines.Clear();
    factory->defines.Set("USE_ZERO_Z_FACE", 1);
    factory->defines.Set("NUM_X_THREADS", numXThreads);
    factory->defines.Set("NUM_Y_THREADS", numYThreads);
    mWriteZFace = factory->CreateFromSource(*msEnforceSource[i]);

    factory->PopDefines();
}

void Fluid3SolvePoisson::Execute(std::shared_ptr<GraphicsEngine> const& engine,
    std::shared_ptr<Texture3> const& divergence)
{
    std::shared_ptr<ComputeShader> solve = mSolvePoisson->GetCShader();
    std::shared_ptr<ComputeShader> xwrite = mWriteXFace->GetCShader();
    std::shared_ptr<ComputeShader> ywrite = mWriteYFace->GetCShader();
    std::shared_ptr<ComputeShader> zwrite = mWriteZFace->GetCShader();

    solve->Set("divergence", divergence);
    engine->Execute(mZeroPoisson, mNumXGroups, mNumYGroups, mNumZGroups);
    for (int i = 0; i < mNumIterations; ++i)
    {
        // Take one step of the Poisson solver.
        solve->Set("poisson", mPoisson0);
        solve->Set("outPoisson", mPoisson1);
        engine->Execute(mSolvePoisson, mNumXGroups, mNumYGroups, mNumZGroups);

        // Set the boundary to zero.
        xwrite->Set("image", mPoisson1);
        engine->Execute(mWriteXFace, 1, mNumYGroups, mNumZGroups);
        ywrite->Set("image", mPoisson1);
        engine->Execute(mWriteYFace, mNumXGroups, 1, mNumZGroups);
        zwrite->Set("image", mPoisson1);
        engine->Execute(mWriteZFace, mNumXGroups, mNumYGroups, 1);

        std::swap(mPoisson0, mPoisson1);
    }
}


std::string const Fluid3SolvePoisson::msGLSLZeroSource =
"layout(r32f) uniform writeonly image3D poisson;\n"
"\n"
"layout (local_size_x = NUM_X_THREADS, local_size_y = NUM_Y_THREADS, local_size_z = NUM_Z_THREADS) in;\n"
"void main()\n"
"{\n"
"    ivec3 c = ivec3(gl_GlobalInvocationID.xyz);\n"
"    imageStore(poisson, c, vec4(0.0f, 0.0f, 0.0f, 0.0f));\n"
"}\n";

std::string const Fluid3SolvePoisson::msGLSLSolveSource =
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
"layout(r32f) uniform readonly image3D divergence;\n"
"layout(r32f) uniform readonly image3D poisson;\n"
"layout(r32f) uniform writeonly image3D outPoisson;\n"
"\n"
"layout (local_size_x = NUM_X_THREADS, local_size_y = NUM_Y_THREADS, local_size_z = NUM_Z_THREADS) in;\n"
"void main()\n"
"{\n"
"    ivec3 c = ivec3(gl_GlobalInvocationID.xyz);\n"
"    ivec3 dim = imageSize(divergence);\n"
"\n"
"    int x = int(c.x);\n"
"    int y = int(c.y);\n"
"    int z = int(c.z);\n"
"    int xm = max(x - 1, 0);\n"
"    int xp = min(x + 1, dim.x - 1);\n"
"    int ym = max(y - 1, 0);\n"
"    int yp = min(y + 1, dim.y - 1);\n"
"    int zm = max(z - 1, 0);\n"
"    int zp = min(z + 1, dim.z - 1);\n"
"\n"
"    // Sample the divergence at (x,y,z).\n"
"    float div = imageLoad(divergence, c).x;\n"
"\n"
"    // Sample Poisson values at immediate neighbors of (x,y,z).\n"
"    float poisPZZ = imageLoad(poisson, ivec3(xp, y, z)).x;\n"
"    float poisMZZ = imageLoad(poisson, ivec3(xm, y, z)).x;\n"
"    float poisZPZ = imageLoad(poisson, ivec3(x, yp, z)).x;\n"
"    float poisZMZ = imageLoad(poisson, ivec3(x, ym, z)).x;\n"
"    float poisZZP = imageLoad(poisson, ivec3(x, y, zp)).x;\n"
"    float poisZZM = imageLoad(poisson, ivec3(x, y, zm)).x;\n"
"\n"
"    vec4 temp = vec4(poisPZZ + poisMZZ, poisZPZ + poisZMZ, poisZZP + poisZZM, div);\n"
"    float outPoissonValue = dot(epsilon, temp);\n"
"    imageStore(outPoisson, c, vec4(outPoissonValue, 0.0f, 0.0f, 0.0f));\n"
"}\n";

std::string const Fluid3SolvePoisson::msGLSLEnforceSource =
"#if USE_ZERO_X_FACE\n"
"layout(r32f) uniform writeonly image3D image;\n"
"\n"
"layout (local_size_x = 1, local_size_y = NUM_Y_THREADS, local_size_z = NUM_Z_THREADS) in;\n"
"void main()\n"
"{\n"
"    ivec3 c = ivec3(gl_GlobalInvocationID.xyz);\n"
"    ivec3 dim = imageSize(image);\n"
"    imageStore(image, ivec3(0, c.y, c.z), vec4(0.0f));\n"
"    imageStore(image, ivec3(dim.x - 1, c.y, c.z), vec4(0.0f));\n"
"}\n"
"#endif\n"
"\n"
"#if USE_ZERO_Y_FACE\n"
"layout(r32f) uniform writeonly image3D image;\n"
"\n"
"layout (local_size_x = NUM_X_THREADS, local_size_y = 1, local_size_z = NUM_Z_THREADS) in;\n"
"void main()\n"
"{\n"
"    ivec3 c = ivec3(gl_GlobalInvocationID.xyz);\n"
"    ivec3 dim = imageSize(image);\n"
"    imageStore(image, ivec3(c.x, 0, c.z), vec4(0.0f));\n"
"    imageStore(image, ivec3(c.x, dim.y - 1, c.z), vec4(0.0f));\n"
"}\n"
"#endif\n"
"\n"
"#if USE_ZERO_Z_FACE\n"
"layout(r32f) uniform writeonly image3D image;\n"
"\n"
"layout (local_size_x = NUM_X_THREADS, local_size_y = NUM_Y_THREADS, local_size_z = 1) in;\n"
"void main()\n"
"{\n"
"    ivec3 c = ivec3(gl_GlobalInvocationID.xyz);\n"
"    ivec3 dim = imageSize(image);\n"
"    imageStore(image, ivec3(c.x, c.y, 0), vec4(0.0f));\n"
"    imageStore(image, ivec3(c.x, c.y, dim.z - 1), vec4(0.0f));\n"
"}\n"
"#endif\n";

std::string const Fluid3SolvePoisson::msHLSLZeroSource =
"RWTexture3D<float> poisson;\n"
"\n"
"[numthreads(NUM_X_THREADS, NUM_Y_THREADS, NUM_Z_THREADS)]\n"
"void CSMain(uint3 c : SV_DispatchThreadID)\n"
"{\n"
"    poisson[c.xyz] = 0.0f;\n"
"}\n";

std::string const Fluid3SolvePoisson::msHLSLSolveSource =
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
"Texture3D<float> divergence;\n"
"Texture3D<float> poisson;\n"
"RWTexture3D<float> outPoisson;\n"
"\n"
"[numthreads(NUM_X_THREADS, NUM_Y_THREADS, NUM_Z_THREADS)]\n"
"void CSMain(uint3 c : SV_DispatchThreadID)\n"
"{\n"
"    uint3 dim;\n"
"    divergence.GetDimensions(dim.x, dim.y, dim.z);\n"
"\n"
"    int x = int(c.x);\n"
"    int y = int(c.y);\n"
"    int z = int(c.z);\n"
"    int xm = max(x - 1, 0);\n"
"    int xp = min(x + 1, dim.x - 1);\n"
"    int ym = max(y - 1, 0);\n"
"    int yp = min(y + 1, dim.y - 1);\n"
"    int zm = max(z - 1, 0);\n"
"    int zp = min(z + 1, dim.z - 1);\n"
"\n"
"    // Sample the divergence at (x,y,z).\n"
"    float div = divergence[c];\n"
"\n"
"    // Sample Poisson values at immediate neighbors of (x,y,z).\n"
"    float poisPZZ = poisson[int3(xp, y, z)];\n"
"    float poisMZZ = poisson[int3(xm, y, z)];\n"
"    float poisZPZ = poisson[int3(x, yp, z)];\n"
"    float poisZMZ = poisson[int3(x, ym, z)];\n"
"    float poisZZP = poisson[int3(x, y, zp)];\n"
"    float poisZZM = poisson[int3(x, y, zm)];\n"
"\n"
"    float4 temp = float4(poisPZZ + poisMZZ, poisZPZ + poisZMZ, poisZZP + poisZZM, div);\n"
"    outPoisson[c] = dot(epsilon, temp);\n"
"}\n";

std::string const Fluid3SolvePoisson::msHLSLEnforceSource =
"#if USE_ZERO_X_FACE\n"
"RWTexture3D<float> image;\n"
"\n"
"[numthreads(1, NUM_Y_THREADS, NUM_Z_THREADS)]\n"
"void CSMain(uint3 c : SV_DispatchThreadID)\n"
"{\n"
"    uint3 dim;\n"
"    image.GetDimensions(dim.x, dim.y, dim.z);\n"
"    image[uint3(0, c.y, c.z)] = 0.0f;\n"
"    image[uint3(dim.x - 1, c.y, c.z)] = 0.0f;\n"
"}\n"
"#endif\n"
"\n"
"#if USE_ZERO_Y_FACE\n"
"RWTexture3D<float> image;\n"
"\n"
"[numthreads(NUM_X_THREADS, 1, NUM_Z_THREADS)]\n"
"void CSMain(uint3 c : SV_DispatchThreadID)\n"
"{\n"
"    uint3 dim;\n"
"    image.GetDimensions(dim.x, dim.y, dim.z);\n"
"    image[uint3(c.x, 0, c.z)] = 0.0f;\n"
"    image[uint3(c.x, dim.y - 1, c.z)] = 0.0f;\n"
"}\n"
"#endif\n"
"\n"
"#if USE_ZERO_Z_FACE\n"
"RWTexture3D<float> image;\n"
"\n"
"[numthreads(NUM_X_THREADS, NUM_Y_THREADS, 1)]\n"
"void CSMain(uint3 c : SV_DispatchThreadID)\n"
"{\n"
"    uint3 dim;\n"
"    image.GetDimensions(dim.x, dim.y, dim.z);\n"
"    image[uint3(c.x, c.y, 0)] = 0.0f;\n"
"    image[uint3(c.x, c.y, dim.z - 1)] = 0.0f;\n"
"}\n"
"#endif\n";

std::string const* Fluid3SolvePoisson::msZeroSource[] =
{
    &msGLSLZeroSource,
    &msHLSLZeroSource
};

std::string const* Fluid3SolvePoisson::msSolveSource[] =
{
    &msGLSLSolveSource,
    &msHLSLSolveSource
};

std::string const* Fluid3SolvePoisson::msEnforceSource[] =
{
    &msGLSLEnforceSource,
    &msHLSLEnforceSource
};
