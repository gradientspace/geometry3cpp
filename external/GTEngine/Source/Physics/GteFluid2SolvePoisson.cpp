// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#include <GTEnginePCH.h>
#include <Physics/GteFluid2SolvePoisson.h>
#include <Graphics/GteGraphicsEngine.h>
#include <Graphics/GteProgramFactory.h>
using namespace gte;

Fluid2SolvePoisson::Fluid2SolvePoisson(std::shared_ptr<ProgramFactory> const& factory,
    int xSize, int ySize, int numXThreads, int numYThreads,
    std::shared_ptr<ConstantBuffer> const& parameters, int numIterations)
    :
    mNumXGroups(xSize/numXThreads),
    mNumYGroups(ySize/numYThreads),
    mNumIterations(numIterations)
{
    mPoisson0 = std::make_shared<Texture2>(DF_R32_FLOAT, xSize, ySize);
    mPoisson0->SetUsage(Resource::SHADER_OUTPUT);
    mPoisson1 = std::make_shared<Texture2>(DF_R32_FLOAT, xSize, ySize);
    mPoisson1->SetUsage(Resource::SHADER_OUTPUT);

    int i = factory->GetAPI();
    factory->PushDefines();
    factory->defines.Set("NUM_X_THREADS", numXThreads);
    factory->defines.Set("NUM_Y_THREADS", numYThreads);

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
    factory->defines.Set("USE_ZERO_X_EDGE", 1);
    factory->defines.Set("NUM_Y_THREADS", numYThreads);
    mWriteXEdge = factory->CreateFromSource(*msEnforceSource[i]);

    factory->defines.Clear();
    factory->defines.Set("USE_ZERO_Y_EDGE", 1);
    factory->defines.Set("NUM_X_THREADS", numXThreads);
    mWriteYEdge = factory->CreateFromSource(*msEnforceSource[i]);

    factory->PopDefines();
}

void Fluid2SolvePoisson::Execute(std::shared_ptr<GraphicsEngine> const& engine,
    std::shared_ptr<Texture2> const& divergence)
{
    std::shared_ptr<ComputeShader> solve = mSolvePoisson->GetCShader();
    std::shared_ptr<ComputeShader> xwrite = mWriteXEdge->GetCShader();
    std::shared_ptr<ComputeShader> ywrite = mWriteYEdge->GetCShader();

    solve->Set("divergence", divergence);
    engine->Execute(mZeroPoisson, mNumXGroups, mNumYGroups, 1);
    for (int i = 0; i < mNumIterations; ++i)
    {
        // Take one step of the Poisson solver.
        solve->Set("poisson", mPoisson0);
        solve->Set("outPoisson", mPoisson1);
        engine->Execute(mSolvePoisson, mNumXGroups, mNumYGroups, 1);

        // Set the boundary to zero.
        xwrite->Set("image", mPoisson1);
        engine->Execute(mWriteXEdge, 1, mNumYGroups, 1);
        ywrite->Set("image", mPoisson1);
        engine->Execute(mWriteYEdge, mNumXGroups, 1, 1);

        std::swap(mPoisson0, mPoisson1);
    }
}


std::string const Fluid2SolvePoisson::msGLSLZeroSource =
"layout(r32f) uniform writeonly image2D poisson;\n"
"\n"
"layout (local_size_x = NUM_X_THREADS, local_size_y = NUM_Y_THREADS, local_size_z = 1) in;\n"
"void main()\n"
"{\n"
"    ivec2 c = ivec2(gl_GlobalInvocationID.xy);\n"
"    imageStore(poisson, c, vec4(0.0f, 0.0f, 0.0f, 0.0f));\n"
"}\n";

std::string const Fluid2SolvePoisson::msGLSLSolveSource =
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
"layout(r32f) uniform readonly image2D divergence;\n"
"layout(r32f) uniform readonly image2D poisson;\n"
"layout(r32f) uniform writeonly image2D outPoisson;\n"
"\n"
"layout (local_size_x = NUM_X_THREADS, local_size_y = NUM_Y_THREADS, local_size_z = 1) in;\n"
"void main()\n"
"{\n"
"    ivec2 c = ivec2(gl_GlobalInvocationID.xy);\n"
"    ivec2 dim = imageSize(divergence);\n"
"\n"
"    int x = int(c.x);\n"
"    int y = int(c.y);\n"
"    int xm = max(x - 1, 0);\n"
"    int xp = min(x + 1, dim.x - 1);\n"
"    int ym = max(y - 1, 0);\n"
"    int yp = min(y + 1, dim.y - 1);\n"
"\n"
"    // Sample the divergence at (x,y).\n"
"    float div = imageLoad(divergence, c).x;\n"
"\n"
"    // Sample Poisson values at (x+dx,y), (x-dx,y), (x,y+dy), (x,y-dy).\n"
"    float poisPZ = imageLoad(poisson, ivec2(xp, y)).x;\n"
"    float poisMZ = imageLoad(poisson, ivec2(xm, y)).x;\n"
"    float poisZP = imageLoad(poisson, ivec2(x, yp)).x;\n"
"    float poisZM = imageLoad(poisson, ivec2(x, ym)).x;\n"
"\n"
"    vec4 temp = vec4(poisPZ + poisMZ, poisZP + poisZM, 0.0f, div);\n"
"    float outPoissonValue = dot(epsilon, temp);\n"
"    imageStore(outPoisson, c, vec4(outPoissonValue, 0.0f, 0.0f, 0.0f));\n"
"}\n";

std::string const Fluid2SolvePoisson::msGLSLEnforceSource =
"#if USE_ZERO_X_EDGE\n"
"layout(r32f) uniform writeonly image2D image;\n"
"\n"
"layout (local_size_x = 1, local_size_y = NUM_Y_THREADS, local_size_z = 1) in;\n"
"void main()\n"
"{\n"
"    ivec2 c = ivec2(gl_GlobalInvocationID.xy);\n"
"    ivec2 dim = imageSize(image);\n"
"    imageStore(image, ivec2(0, c.y), vec4(0.0f, 0.0f, 0.0f, 0.0f));\n"
"    imageStore(image, ivec2(dim.x - 1, c.y), vec4(0.0f, 0.0f, 0.0f, 0.0f));\n"
"}\n"
"#endif\n"
"\n"
"#if USE_ZERO_Y_EDGE\n"
"layout(r32f) uniform writeonly image2D image;\n"
"\n"
"layout (local_size_x = NUM_X_THREADS, local_size_y = 1, local_size_z = 1) in;\n"
"void main()\n"
"{\n"
"    ivec2 c = ivec2(gl_GlobalInvocationID.xy);\n"
"    ivec2 dim = imageSize(image);\n"
"    imageStore(image, ivec2(c.x, 0), vec4(0.0f, 0.0f, 0.0f, 0.0f));\n"
"    imageStore(image, ivec2(c.x, dim.y - 1), vec4(0.0f, 0.0f, 0.0f, 0.0f));\n"
"}\n"
"#endif\n";

std::string const Fluid2SolvePoisson::msHLSLZeroSource =
"RWTexture2D<float> poisson;\n"
"\n"
"[numthreads(NUM_X_THREADS, NUM_Y_THREADS, 1)]\n"
"void CSMain(uint2 c : SV_DispatchThreadID)\n"
"{\n"
"    poisson[c] = 0.0f;\n"
"}\n";

std::string const Fluid2SolvePoisson::msHLSLSolveSource =
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
"Texture2D<float> divergence;\n"
"Texture2D<float> poisson;\n"
"RWTexture2D<float> outPoisson;\n"
"\n"
"[numthreads(NUM_X_THREADS, NUM_Y_THREADS, 1)]\n"
"void CSMain(uint2 c : SV_DispatchThreadID)\n"
"{\n"
"    uint2 dim;\n"
"    divergence.GetDimensions(dim.x, dim.y);\n"
"\n"
"    int x = int(c.x);\n"
"    int y = int(c.y);\n"
"    int xm = max(x - 1, 0);\n"
"    int xp = min(x + 1, dim.x - 1);\n"
"    int ym = max(y - 1, 0);\n"
"    int yp = min(y + 1, dim.y - 1);\n"
"\n"
"    // Sample the divergence at (x,y).\n"
"    float div = divergence[int2(x, y)];\n"
"\n"
"    // Sample Poisson values at (x+dx,y), (x-dx,y), (x,y+dy), (x,y-dy).\n"
"    float poisPZ = poisson[int2(xp, y)];\n"
"    float poisMZ = poisson[int2(xm, y)];\n"
"    float poisZP = poisson[int2(x, yp)];\n"
"    float poisZM = poisson[int2(x, ym)];\n"
"\n"
"    float4 temp = float4(poisPZ + poisMZ, poisZP + poisZM, 0.0f, div);\n"
"    outPoisson[c] = dot(epsilon, temp);\n"
"}\n";

std::string const Fluid2SolvePoisson::msHLSLEnforceSource =
"#if USE_ZERO_X_EDGE\n"
"RWTexture2D<float> image;\n"
"\n"
"[numthreads(1, NUM_Y_THREADS, 1)]\n"
"void CSMain(uint2 c : SV_DispatchThreadID)\n"
"{\n"
"    uint2 dim;\n"
"    image.GetDimensions(dim.x, dim.y);\n"
"    image[uint2(0, c.y)] = 0.0f;\n"
"    image[uint2(dim.x - 1, c.y)] = 0.0f;\n"
"}\n"
"#endif\n"
"\n"
"#if USE_ZERO_Y_EDGE\n"
"RWTexture2D<float> image;\n"
"\n"
"[numthreads(NUM_X_THREADS, 1, 1)]\n"
"void CSMain(uint2 c : SV_DispatchThreadID)\n"
"{\n"
"    uint2 dim;\n"
"    image.GetDimensions(dim.x, dim.y);\n"
"    image[uint2(c.x, 0)] = 0.0f;\n"
"    image[uint2(c.x, dim.y - 1)] = 0.0f;\n"
"}\n"
"#endif\n";

std::string const* Fluid2SolvePoisson::msZeroSource[] =
{
    &msGLSLZeroSource,
    &msHLSLZeroSource
};

std::string const* Fluid2SolvePoisson::msSolveSource[] =
{
    &msGLSLSolveSource,
    &msHLSLSolveSource
};

std::string const* Fluid2SolvePoisson::msEnforceSource[] =
{
    &msGLSLEnforceSource,
    &msHLSLEnforceSource
};
