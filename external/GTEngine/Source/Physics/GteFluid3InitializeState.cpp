// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#include <GTEnginePCH.h>
#include <Physics/GteFluid3InitializeState.h>
#include <Graphics/GteGraphicsEngine.h>
#include <Graphics/GteProgramFactory.h>
#include <random>
using namespace gte;

Fluid3InitializeState::Fluid3InitializeState(std::shared_ptr<ProgramFactory> const& factory,
    int xSize, int ySize, int zSize, int numXThreads, int numYThreads, int numZThreads)
    :
    mNumXGroups(xSize/numXThreads),
    mNumYGroups(ySize/numYThreads),
    mNumZGroups(zSize/numZThreads)
{
    // Use a Mersenne twister engine for random numbers.
    std::mt19937 mte;
    std::uniform_real_distribution<float> unirnd(0.0f, 1.0f);

    // Initial density values are randomly generated.
    mDensity = std::make_shared<Texture3>(DF_R32_FLOAT, xSize, ySize, zSize);
    float* data = mDensity->Get<float>();
    for (unsigned int i = 0; i < mDensity->GetNumElements(); ++i, ++data)
    {
        *data = unirnd(mte);
    }

    // Initial velocity values are zero.
    mVelocity = std::make_shared<Texture3>(DF_R32G32B32A32_FLOAT, xSize, ySize, zSize);
    mVelocity->SetUsage(Resource::SHADER_OUTPUT);
    memset(mVelocity->GetData(), 0, mVelocity->GetNumBytes());

    mStateTm1 = std::make_shared<Texture3>(DF_R32G32B32A32_FLOAT, xSize, ySize, zSize);
    mStateTm1->SetUsage(Resource::SHADER_OUTPUT);

    mStateT = std::make_shared<Texture3>(DF_R32G32B32A32_FLOAT, xSize, ySize, zSize);
    mStateT->SetUsage(Resource::SHADER_OUTPUT);

    // Create the shader for generating velocity from vortices.
    int i = factory->GetAPI();
    factory->PushDefines();
    factory->defines.Set("NUM_X_THREADS", numXThreads);
    factory->defines.Set("NUM_Y_THREADS", numYThreads);
    factory->defines.Set("NUM_Z_THREADS", numZThreads);
    mInitializeState = factory->CreateFromSource(*msSource[i]);
    if (mInitializeState)
    {
        std::shared_ptr<ComputeShader> cshader = mInitializeState->GetCShader();
        cshader->Set("density", mDensity);
        cshader->Set("velocity", mVelocity);
        cshader->Set("stateTm1", mStateTm1);
        cshader->Set("stateT", mStateT);
    }
    factory->PopDefines();
}

void Fluid3InitializeState::Execute(std::shared_ptr<GraphicsEngine> const& engine)
{
    engine->Execute(mInitializeState, mNumXGroups, mNumYGroups, mNumZGroups);
}


std::string const Fluid3InitializeState::msGLSLSource =
"layout(r32f) uniform readonly image3D density;\n"
"layout(rgba32f) uniform readonly image3D velocity;\n"
"layout(rgba32f) uniform writeonly image3D stateTm1;\n"
"layout(rgba32f) uniform writeonly image3D stateT;\n"
"\n"
"layout (local_size_x = NUM_X_THREADS, local_size_y = NUM_Y_THREADS, local_size_z = NUM_Z_THREADS) in;\n"
"void main()\n"
"{\n"
"    ivec3 c = ivec3(gl_GlobalInvocationID.xyz);\n"
"    vec4 initial = vec4(imageLoad(velocity, c).xyz, imageLoad(density, c).x);\n"
"    imageStore(stateTm1, c, initial);\n"
"    imageStore(stateT, c, initial);\n"
"}\n";

std::string const Fluid3InitializeState::msHLSLSource =
"Texture3D<float> density;\n"
"Texture3D<float4> velocity;\n"
"RWTexture3D<float4> stateTm1;\n"
"RWTexture3D<float4> stateT;\n"
"\n"
"[numthreads(NUM_X_THREADS, NUM_Y_THREADS, NUM_Z_THREADS)]\n"
"void CSMain(uint3 c : SV_DispatchThreadID)\n"
"{\n"
"    float4 initial = float4(velocity[c].xyz, density[c]);\n"
"    stateTm1[c] = initial;\n"
"    stateT[c] = initial;\n"
"}\n";

std::string const* Fluid3InitializeState::msSource[] =
{
    &msGLSLSource,
    &msHLSLSource
};
