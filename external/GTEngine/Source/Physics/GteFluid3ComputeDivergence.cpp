// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#include <GTEnginePCH.h>
#include <Physics/GteFluid3ComputeDivergence.h>
#include <Graphics/GteGraphicsEngine.h>
#include <Graphics/GteProgramFactory.h>
using namespace gte;

Fluid3ComputeDivergence::Fluid3ComputeDivergence(std::shared_ptr<ProgramFactory> const& factory,
    int xSize, int ySize, int zSize, int numXThreads, int numYThreads,
    int numZThreads, std::shared_ptr<ConstantBuffer> const& parameters)
    :
    mNumXGroups(xSize/numXThreads),
    mNumYGroups(ySize/numYThreads),
    mNumZGroups(zSize/numZThreads)
{
    mDivergence = std::make_shared<Texture3>(DF_R32_FLOAT, xSize, ySize, zSize);
    mDivergence->SetUsage(Resource::SHADER_OUTPUT);

    int i = factory->GetAPI();
    factory->PushDefines();
    factory->defines.Set("NUM_X_THREADS", numXThreads);
    factory->defines.Set("NUM_Y_THREADS", numYThreads);
    factory->defines.Set("NUM_Z_THREADS", numZThreads);
    mComputeDivergence = factory->CreateFromSource(*msSource[i]);
    if (mComputeDivergence)
    {
        mComputeDivergence->GetCShader()->Set("Parameters", parameters);
        mComputeDivergence->GetCShader()->Set("divergence", mDivergence);
    }

    factory->PopDefines();
}

void Fluid3ComputeDivergence::Execute(std::shared_ptr<GraphicsEngine> const& engine,
    std::shared_ptr<Texture3> const& state)
{
    std::shared_ptr<ComputeShader> cshader = mComputeDivergence->GetCShader();
    cshader->Set("state", state);
    engine->Execute(mComputeDivergence, mNumXGroups, mNumYGroups, mNumZGroups);
}


std::string const Fluid3ComputeDivergence::msGLSLSource =
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
"layout(rgba32f) uniform readonly image3D state;\n"
"layout(r32f) uniform writeonly image3D divergence;\n"
"\n"
"layout (local_size_x = NUM_X_THREADS, local_size_y = NUM_Y_THREADS, local_size_z = NUM_Z_THREADS) in;\n"
"void main()\n"
"{\n"
"    ivec3 c = ivec3(gl_GlobalInvocationID.xyz);\n"
"    ivec3 dim = imageSize(state);\n"
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
"    vec3 velocityGradient = vec3\n"
"    (\n"
"        imageLoad(state, ivec3(xp, y, z)).x - imageLoad(state, ivec3(xm, y, z)).x,\n"
"        imageLoad(state, ivec3(x, yp, z)).y - imageLoad(state, ivec3(x, ym, z)).y,\n"
"        imageLoad(state, ivec3(x, y, zp)).z - imageLoad(state, ivec3(x, y, zm)).z\n"
"    );\n"
"\n"
"    float divergenceValue = dot(halfDivDelta.xyz, velocityGradient);\n"
"    imageStore(divergence, c, vec4(divergenceValue, 0.0f, 0.0f, 0.0f));\n"
"}\n";

std::string const Fluid3ComputeDivergence::msHLSLSource =
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
"Texture3D<float4> state;\n"
"RWTexture3D<float> divergence;\n"
"\n"
"[numthreads(NUM_X_THREADS, NUM_Y_THREADS, NUM_Z_THREADS)]\n"
"void CSMain(uint3 c : SV_DispatchThreadID)\n"
"{\n"
"    uint3 dim;\n"
"    state.GetDimensions(dim.x, dim.y, dim.z);\n"
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
"    float3 velocityGradient = float3\n"
"    (\n"
"        state[int3(xp, y, z)].x - state[int3(xm, y, z)].x,\n"
"        state[int3(x, yp, z)].y - state[int3(x, ym, z)].y,\n"
"        state[int3(x, y, zp)].z - state[int3(x, y, zm)].z\n"
"    );\n"
"\n"
"    divergence[c] = dot(halfDivDelta.xyz, velocityGradient);\n"
"}\n";

std::string const* Fluid3ComputeDivergence::msSource[] =
{
    &msGLSLSource,
    &msHLSLSource
};
