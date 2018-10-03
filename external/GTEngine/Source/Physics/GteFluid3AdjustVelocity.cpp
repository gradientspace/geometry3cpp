// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#include <GTEnginePCH.h>
#include <Physics/GteFluid3AdjustVelocity.h>
#include <Graphics/GteGraphicsEngine.h>
#include <Graphics/GteProgramFactory.h>
using namespace gte;

Fluid3AdjustVelocity::Fluid3AdjustVelocity(std::shared_ptr<ProgramFactory> const& factory,
    int xSize, int ySize, int zSize, int numXThreads, int numYThreads, int numZThreads,
    std::shared_ptr<ConstantBuffer> const& parameters)
    :
    mNumXGroups(xSize/numXThreads),
    mNumYGroups(ySize/numYThreads),
    mNumZGroups(zSize/numZThreads)
{
    int i = factory->GetAPI();
    factory->PushDefines();
    factory->defines.Set("NUM_X_THREADS", numXThreads);
    factory->defines.Set("NUM_Y_THREADS", numYThreads);
    factory->defines.Set("NUM_Z_THREADS", numZThreads);

    mAdjustVelocity = factory->CreateFromSource(*msSource[i]);
    if (mAdjustVelocity)
    {
        mAdjustVelocity->GetCShader()->Set("Parameters", parameters);
    }

    factory->PopDefines();
}

void Fluid3AdjustVelocity::Execute(std::shared_ptr<GraphicsEngine> const& engine,
    std::shared_ptr<Texture3> const& inState,
    std::shared_ptr<Texture3> const& poisson,
    std::shared_ptr<Texture3> const& outState)
{
    std::shared_ptr<ComputeShader> cshader = mAdjustVelocity->GetCShader();
    cshader->Set("inState", inState);
    cshader->Set("poisson", poisson);
    cshader->Set("outState", outState);
    engine->Execute(mAdjustVelocity, mNumXGroups, mNumYGroups, mNumZGroups);
}


std::string const Fluid3AdjustVelocity::msGLSLSource =
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
"layout(rgba32f) uniform readonly image3D inState;\n"
"layout(r32f) uniform readonly image3D poisson;\n"
"layout(rgba32f) uniform writeonly image3D outState;\n"
"\n"
"layout (local_size_x = NUM_X_THREADS, local_size_y = NUM_Y_THREADS, local_size_z = NUM_Z_THREADS) in;\n"
"void main()\n"
"{\n"
"    ivec3 c = ivec3(gl_GlobalInvocationID.xyz);\n"
"    ivec3 dim = imageSize(inState);\n"
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
"    // Sample the state at (x,y,z).\n"
"    vec4 state = imageLoad(inState, c);\n"
"\n"
"    // Sample Poisson values at immediate neighbors of (x,y,z).\n"
"    float poisPZZ = imageLoad(poisson, ivec3(xp, y, z)).x;\n"
"    float poisMZZ = imageLoad(poisson, ivec3(xm, y, z)).x;\n"
"    float poisZPZ = imageLoad(poisson, ivec3(x, yp, z)).x;\n"
"    float poisZMZ = imageLoad(poisson, ivec3(x, ym, z)).x;\n"
"    float poisZZP = imageLoad(poisson, ivec3(x, y, zp)).x;\n"
"    float poisZZM = imageLoad(poisson, ivec3(x, y, zm)).x;\n"
"\n"
"    vec4 diff = vec4(poisPZZ - poisMZZ, poisZPZ - poisZMZ, poisZZP - poisZZM, 0.0f);\n"
"    imageStore(outState, c, state + halfDivDelta*diff);\n"
"}\n";

std::string const Fluid3AdjustVelocity::msHLSLSource =
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
"Texture3D<float4> inState;\n"
"Texture3D<float> poisson;\n"
"RWTexture3D<float4> outState;\n"
"\n"
"[numthreads(NUM_X_THREADS, NUM_Y_THREADS, NUM_Z_THREADS)]\n"
"void CSMain(uint3 c : SV_DispatchThreadID)\n"
"{\n"
"    uint3 dim;\n"
"    inState.GetDimensions(dim.x, dim.y, dim.z);\n"
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
"    // Sample the state at (x,y,z).\n"
"    float4 state = inState[c];\n"
"\n"
"    // Sample Poisson values at immediate neighbors of (x,y,z).\n"
"    float poisPZZ = poisson[int3(xp, y, z)];\n"
"    float poisMZZ = poisson[int3(xm, y, z)];\n"
"    float poisZPZ = poisson[int3(x, yp, z)];\n"
"    float poisZMZ = poisson[int3(x, ym, z)];\n"
"    float poisZZP = poisson[int3(x, y, zp)];\n"
"    float poisZZM = poisson[int3(x, y, zm)];\n"
"\n"
"    float4 diff = float4(poisPZZ - poisMZZ, poisZPZ - poisZMZ, poisZZP - poisZZM, 0.0f);\n"
"    outState[c] = state + halfDivDelta*diff;\n"
"}\n";

std::string const* Fluid3AdjustVelocity::msSource[] =
{
    &msGLSLSource,
    &msHLSLSource
};
