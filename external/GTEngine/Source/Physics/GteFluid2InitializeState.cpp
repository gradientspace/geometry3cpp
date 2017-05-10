// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#include <GTEnginePCH.h>
#include <Physics/GteFluid2InitializeState.h>
#include <Graphics/GteGraphicsEngine.h>
#include <Graphics/GteProgramFactory.h>
#include <random>
using namespace gte;

Fluid2InitializeState::Fluid2InitializeState(std::shared_ptr<ProgramFactory> const& factory,
    int xSize, int ySize, int numXThreads, int numYThreads)
    :
    mNumXGroups(xSize/numXThreads),
    mNumYGroups(ySize/numYThreads)
{
    // Use a Mersenne twister engine for random numbers.
    std::mt19937 mte;
    std::uniform_real_distribution<float> unirnd(0.0f, 1.0f);

    // Initial density values are randomly generated.
    mDensity = std::make_shared<Texture2>(DF_R32_FLOAT, xSize, ySize);
    float* data = mDensity->Get<float>();
    for (unsigned int i = 0; i < mDensity->GetNumElements(); ++i, ++data)
    {
        *data = unirnd(mte);
    }

    // Initial velocity values are zero.
    mVelocity = std::make_shared<Texture2>(DF_R32G32_FLOAT, xSize, ySize);
    memset(mVelocity->GetData(), 0, mVelocity->GetNumBytes());

    // The states at time 0 and time -dt are initialized by a compute shader.
    mStateTm1 = std::make_shared<Texture2>(DF_R32G32B32A32_FLOAT, xSize, ySize);
    mStateTm1->SetUsage(Resource::SHADER_OUTPUT);

    mStateT = std::make_shared<Texture2>(DF_R32G32B32A32_FLOAT, xSize, ySize);
    mStateT->SetUsage(Resource::SHADER_OUTPUT);

    // Create the shader for initializing velocity and density.
    int i = factory->GetAPI();
    factory->PushDefines();
    factory->defines.Set("NUM_X_THREADS", numXThreads);
    factory->defines.Set("NUM_Y_THREADS", numYThreads);
    mInitializeState = factory->CreateFromSource(*msSource[i]);
    if (mInitializeState)
    {
        std::shared_ptr<ComputeShader> cshader =  mInitializeState->GetCShader();
        cshader->Set("density", mDensity);
        cshader->Set("velocity", mVelocity);
        cshader->Set("stateTm1", mStateTm1);
        cshader->Set("stateT", mStateT);
    }
    factory->PopDefines();
}

void Fluid2InitializeState::Execute(std::shared_ptr<GraphicsEngine> const& engine)
{
    engine->Execute(mInitializeState, mNumXGroups, mNumYGroups, 1);
}


std::string const Fluid2InitializeState::msGLSLSource =
"layout(r32f) uniform readonly image2D density;\n"
"layout(rg32f) uniform readonly image2D velocity;\n"
"layout(rgba32f) uniform writeonly image2D stateTm1;\n"
"layout(rgba32f) uniform writeonly image2D stateT;\n"
"\n"
"layout (local_size_x = NUM_X_THREADS, local_size_y = NUM_Y_THREADS, local_size_z = 1) in;\n"
"void main()\n"
"{\n"
"    ivec2 c = ivec2(gl_GlobalInvocationID.xy);\n"
"    vec4 initial = vec4(imageLoad(velocity, c).xy, 0.0f, imageLoad(density, c).x);\n"
"    imageStore(stateTm1, c, initial);\n"
"    imageStore(stateT, c, initial);\n"
"}\n";

std::string const Fluid2InitializeState::msHLSLSource =
"Texture2D<float> density;\n"
"Texture2D<float2> velocity;\n"
"RWTexture2D<float4> stateTm1;\n"
"RWTexture2D<float4> stateT;\n"
"\n"
"[numthreads(NUM_X_THREADS, NUM_Y_THREADS, 1)]\n"
"void CSMain(uint2 c : SV_DispatchThreadID)\n"
"{\n"
"    float4 initial = float4(velocity[c], 0.0f, density[c]);\n"
"    stateTm1[c] = initial;\n"
"    stateT[c] = initial;\n"
"}\n";

std::string const* Fluid2InitializeState::msSource[] =
{
    &msGLSLSource,
    &msHLSLSource
};
