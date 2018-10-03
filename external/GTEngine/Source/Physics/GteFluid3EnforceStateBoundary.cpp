// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#include <GTEnginePCH.h>
#include <Physics/GteFluid3EnforceStateBoundary.h>
#include <Graphics/GteGraphicsEngine.h>
#include <Graphics/GteProgramFactory.h>
using namespace gte;

Fluid3EnforceStateBoundary::Fluid3EnforceStateBoundary(std::shared_ptr<ProgramFactory> const& factory,
    int xSize, int ySize, int zSize, int numXThreads, int numYThreads, int numZThreads)
    :
    mNumXGroups(xSize/numXThreads),
    mNumYGroups(ySize/numYThreads),
    mNumZGroups(zSize/numZThreads)
{
    mXMin = std::make_shared<Texture2>(DF_R32G32_FLOAT, ySize, zSize);
    mXMin->SetUsage(Resource::SHADER_OUTPUT);
    mXMax = std::make_shared<Texture2>(DF_R32G32_FLOAT, ySize, zSize);
    mXMax->SetUsage(Resource::SHADER_OUTPUT);
    mYMin = std::make_shared<Texture2>(DF_R32G32_FLOAT, xSize, zSize);
    mYMin->SetUsage(Resource::SHADER_OUTPUT);
    mYMax = std::make_shared<Texture2>(DF_R32G32_FLOAT, xSize, zSize);
    mYMax->SetUsage(Resource::SHADER_OUTPUT);
    mZMin = std::make_shared<Texture2>(DF_R32G32_FLOAT, xSize, ySize);
    mZMin->SetUsage(Resource::SHADER_OUTPUT);
    mZMax = std::make_shared<Texture2>(DF_R32G32_FLOAT, xSize, ySize);
    mZMax->SetUsage(Resource::SHADER_OUTPUT);

    int i = factory->GetAPI();
    factory->PushDefines();
    factory->defines.Set("USE_COPY_X_FACE", 1);
    factory->defines.Set("NUM_Y_THREADS", numYThreads);
    factory->defines.Set("NUM_Z_THREADS", numZThreads);
    mCopyXFace = factory->CreateFromSource(*msSource[i]);
    if (mCopyXFace)
    {
        mCopyXFace->GetCShader()->Set("xMin", mXMin);
        mCopyXFace->GetCShader()->Set("xMax", mXMax);
    }

    factory->defines.Clear();
    factory->defines.Set("USE_WRITE_X_FACE", 1);
    factory->defines.Set("NUM_Y_THREADS", numYThreads);
    factory->defines.Set("NUM_Z_THREADS", numZThreads);
    mWriteXFace = factory->CreateFromSource(*msSource[i]);
    if (mWriteXFace)
    {
        mWriteXFace->GetCShader()->Set("xMin", mXMin);
        mWriteXFace->GetCShader()->Set("xMax", mXMax);
    }

    factory->defines.Clear();
    factory->defines.Set("USE_COPY_Y_FACE", 1);
    factory->defines.Set("NUM_X_THREADS", numXThreads);
    factory->defines.Set("NUM_Z_THREADS", numZThreads);
    mCopyYFace = factory->CreateFromSource(*msSource[i]);
    if (mCopyYFace)
    {
        mCopyYFace->GetCShader()->Set("yMin", mYMin);
        mCopyYFace->GetCShader()->Set("yMax", mYMax);
    }

    factory->defines.Clear();
    factory->defines.Set("USE_WRITE_Y_FACE", 1);
    factory->defines.Set("NUM_X_THREADS", numXThreads);
    factory->defines.Set("NUM_Z_THREADS", numZThreads);
    mWriteYFace = factory->CreateFromSource(*msSource[i]);
    if (mWriteYFace)
    {
        mWriteYFace->GetCShader()->Set("yMin", mYMin);
        mWriteYFace->GetCShader()->Set("yMax", mYMax);
    }

    factory->defines.Clear();
    factory->defines.Set("USE_COPY_Z_FACE", 1);
    factory->defines.Set("NUM_X_THREADS", numXThreads);
    factory->defines.Set("NUM_Y_THREADS", numYThreads);
    mCopyZFace = factory->CreateFromSource(*msSource[i]);
    if (mCopyZFace)
    {
        mCopyZFace->GetCShader()->Set("zMin", mZMin);
        mCopyZFace->GetCShader()->Set("zMax", mZMax);
    }

    factory->defines.Clear();
    factory->defines.Set("USE_WRITE_Z_FACE", 1);
    factory->defines.Set("NUM_X_THREADS", numXThreads);
    factory->defines.Set("NUM_Y_THREADS", numYThreads);
    mWriteZFace = factory->CreateFromSource(*msSource[i]);
    if (mWriteZFace)
    {
        mWriteZFace->GetCShader()->Set("zMin", mZMin);
        mWriteZFace->GetCShader()->Set("zMax", mZMax);
    }

    factory->PopDefines();
}

void Fluid3EnforceStateBoundary::Execute(std::shared_ptr<GraphicsEngine> const& engine,
    std::shared_ptr<Texture3> const& state)
{
    // in: state
    // out: mXMin, mXMax
    mCopyXFace->GetCShader()->Set("state", state);
    engine->Execute(mCopyXFace, 1, mNumYGroups, mNumZGroups);

    // in: mXMin, mXMax
    // out: state
    mWriteXFace->GetCShader()->Set("state", state);
    engine->Execute(mWriteXFace, 1, mNumYGroups, mNumZGroups);

    // in: state
    // out: mYMin, mYMax
    mCopyYFace->GetCShader()->Set("state", state);
    engine->Execute(mCopyYFace, mNumXGroups, 1, mNumZGroups);

    // in: mYMin, mYMax
    // out: state
    mWriteYFace->GetCShader()->Set("state", state);
    engine->Execute(mWriteYFace, mNumXGroups, 1, mNumZGroups);

    // in: state
    // out: mZMin, mZMax
    mCopyZFace->GetCShader()->Set("state", state);
    engine->Execute(mCopyZFace, mNumXGroups, mNumYGroups, 1);

    // in: mZMin, mZMax
    // out: state
    mWriteZFace->GetCShader()->Set("state", state);
    engine->Execute(mWriteZFace, mNumXGroups, mNumYGroups, 1);
}


std::string const Fluid3EnforceStateBoundary::msGLSLSource =
"#if USE_COPY_X_FACE\n"
"layout(rgba32f) uniform readonly image3D state;\n"
"layout(rg32f) uniform writeonly image2D xMin;\n"
"layout(rg32f) uniform writeonly image2D xMax;\n"
"\n"
"layout (local_size_x = 1, local_size_y = NUM_Y_THREADS, local_size_z = NUM_Z_THREADS) in;\n"
"void main()\n"
"{\n"
"    ivec3 c = ivec3(gl_GlobalInvocationID.xyz);\n"
"    ivec3 dim = imageSize(state);\n"
"    vec2 xMinValue = imageLoad(state, ivec3(1, c.y, c.z)).yz;\n"
"    vec2 xMaxValue = imageLoad(state, ivec3(dim.x - 2, c.y, c.z)).yz;\n"
"    imageStore(xMin, c.yz, vec4(xMinValue, 0.0f, 0.0f));\n"
"    imageStore(xMax, c.yz, vec4(xMaxValue, 0.0f, 0.0f));\n"
"}\n"
"#endif\n"
"\n"
"#if USE_WRITE_X_FACE\n"
"layout(rg32f) uniform readonly image2D xMin;\n"
"layout(rg32f) uniform readonly image2D xMax;\n"
"layout(rgba32f) uniform writeonly image3D state;\n"
"\n"
"layout (local_size_x = 1, local_size_y = NUM_Y_THREADS, local_size_z = NUM_Z_THREADS) in;\n"
"void main()\n"
"{\n"
"    ivec3 c = ivec3(gl_GlobalInvocationID.xyz);\n"
"    ivec3 dim = imageSize(state);\n"
"    vec2 xMinValue = imageLoad(xMin, c.yz).xy;\n"
"    vec2 xMaxValue = imageLoad(xMax, c.yz).xy;\n"
"    imageStore(state, ivec3(0, c.y, c.z), vec4(0.0f, xMinValue.x, xMinValue.y, 0.0f));\n"
"    imageStore(state, ivec3(dim.x - 1, c.y, c.z), vec4(0.0f, xMaxValue.x, xMaxValue.y, 0.0f));\n"
"}\n"
"#endif\n"
"\n"
"#if USE_COPY_Y_FACE\n"
"layout(rgba32f) uniform readonly image3D state;\n"
"layout(rg32f) uniform writeonly image2D yMin;\n"
"layout(rg32f) uniform writeonly image2D yMax;\n"
"\n"
"layout (local_size_x = NUM_X_THREADS, local_size_y = 1, local_size_z = NUM_Z_THREADS) in;\n"
"void main()\n"
"{\n"
"    ivec3 c = ivec3(gl_GlobalInvocationID.xyz);\n"
"    ivec3 dim = imageSize(state);\n"
"    vec2 yMinValue = imageLoad(state, ivec3(c.x, 1, c.z)).xz;\n"
"    vec2 yMaxValue = imageLoad(state, ivec3(c.x, dim.y - 2, c.z)).xz;\n"
"    imageStore(yMin, c.xz, vec4(yMinValue, 0.0f, 0.0f));\n"
"    imageStore(yMax, c.xz, vec4(yMaxValue, 0.0f, 0.0f));\n"
"}\n"
"#endif\n"
"\n"
"#if USE_WRITE_Y_FACE\n"
"layout(rg32f) uniform readonly image2D yMin;\n"
"layout(rg32f) uniform readonly image2D yMax;\n"
"layout(rgba32f) uniform writeonly image3D state;\n"
"\n"
"layout (local_size_x = NUM_X_THREADS, local_size_y = 1, local_size_z = NUM_Z_THREADS) in;\n"
"void main()\n"
"{\n"
"    ivec3 c = ivec3(gl_GlobalInvocationID.xyz);\n"
"    ivec3 dim = imageSize(state);\n"
"    vec2 yMinValue = imageLoad(yMin, c.xz).xy;\n"
"    vec2 yMaxValue = imageLoad(yMax, c.xz).xy;\n"
"    imageStore(state, ivec3(c.x, 0, c.z), vec4(yMinValue.x, 0.0f, yMinValue.y, 0.0f));\n"
"    imageStore(state, ivec3(c.x, dim.y - 1, c.z), vec4(yMaxValue.x, 0.0f, yMaxValue.y, 0.0f));\n"
"}\n"
"#endif\n"
"\n"
"#if USE_COPY_Z_FACE\n"
"layout(rgba32f) uniform readonly image3D state;\n"
"layout(rg32f) uniform writeonly image2D zMin;\n"
"layout(rg32f) uniform writeonly image2D zMax;\n"
"\n"
"layout (local_size_x = NUM_X_THREADS, local_size_y = NUM_Y_THREADS, local_size_z = 1) in;\n"
"void main()\n"
"{\n"
"    ivec3 c = ivec3(gl_GlobalInvocationID.xyz);\n"
"    ivec3 dim = imageSize(state);\n"
"    vec2 zMinValue = imageLoad(state, ivec3(c.x, c.y, 1)).xy;\n"
"    vec2 zMaxValue = imageLoad(state, ivec3(c.x, c.y, dim.z - 2)).xy;\n"
"    imageStore(zMin, c.xy, vec4(zMinValue, 0.0f, 0.0f));\n"
"    imageStore(zMax, c.xy, vec4(zMaxValue, 0.0f, 0.0f));\n"
"}\n"
"#endif\n"
"\n"
"#if USE_WRITE_Z_FACE\n"
"layout(rg32f) uniform readonly image2D zMin;\n"
"layout(rg32f) uniform readonly image2D zMax;\n"
"layout(rgba32f) uniform writeonly image3D state;\n"
"\n"
"layout (local_size_x = NUM_X_THREADS, local_size_y = NUM_Y_THREADS, local_size_z = 1) in;\n"
"void main()\n"
"{\n"
"    ivec3 c = ivec3(gl_GlobalInvocationID.xyz);\n"
"    ivec3 dim = imageSize(state);\n"
"    vec2 zMinValue = imageLoad(zMin, c.xy).xy;\n"
"    vec2 zMaxValue = imageLoad(zMax, c.xy).xy;\n"
"    imageStore(state, ivec3(c.x, c.y, 0), vec4(zMinValue.x, zMinValue.y, 0.0f, 0.0f));\n"
"    imageStore(state, ivec3(c.x, c.y, dim.z - 1), vec4(zMaxValue.x, zMaxValue.y, 0.0f, 0.0f));\n"
"}\n"
"#endif\n";

std::string const Fluid3EnforceStateBoundary::msHLSLSource =
"#if USE_COPY_X_FACE\n"
"Texture3D<float4> state;\n"
"RWTexture2D<float2> xMin;\n"
"RWTexture2D<float2> xMax;\n"
"\n"
"[numthreads(1, NUM_Y_THREADS, NUM_Z_THREADS)]\n"
"void CSMain(uint3 c : SV_DispatchThreadID)\n"
"{\n"
"    uint3 dim;\n"
"    state.GetDimensions(dim.x, dim.y, dim.z);\n"
"    xMin[c.yz] = state[uint3(1, c.y, c.z)].yz;\n"
"    xMax[c.yz] = state[uint3(dim.x - 2, c.y, c.z)].yz;\n"
"}\n"
"#endif\n"
"\n"
"#if USE_WRITE_X_FACE\n"
"Texture2D<float2> xMin;\n"
"Texture2D<float2> xMax;\n"
"RWTexture3D<float4> state;\n"
"\n"
"[numthreads(1, NUM_Y_THREADS, NUM_Z_THREADS)]\n"
"void CSMain(uint3 c : SV_DispatchThreadID)\n"
"{\n"
"    uint3 dim;\n"
"    state.GetDimensions(dim.x, dim.y, dim.z);\n"
"    state[uint3(0, c.y, c.z)] = float4(0.0f, xMin[c.yz].x, xMin[c.yz].y, 0.0f);\n"
"    state[uint3(dim.x - 1, c.y, c.z)] = float4(0.0f, xMax[c.yz].x, xMax[c.yz].y, 0.0f);\n"
"}\n"
"#endif\n"
"\n"
"#if USE_COPY_Y_FACE\n"
"Texture3D<float4> state;\n"
"RWTexture2D<float2> yMin;\n"
"RWTexture2D<float2> yMax;\n"
"\n"
"[numthreads(NUM_X_THREADS, 1, NUM_Z_THREADS)]\n"
"void CSMain(uint3 c : SV_DispatchThreadID)\n"
"{\n"
"    uint3 dim;\n"
"    state.GetDimensions(dim.x, dim.y, dim.z);\n"
"    yMin[c.xz] = state[uint3(c.x, 1, c.z)].xz;\n"
"    yMax[c.xz] = state[uint3(c.x, dim.y - 2, c.z)].xz;\n"
"}\n"
"#endif\n"
"\n"
"#if USE_WRITE_Y_FACE\n"
"Texture2D<float2> yMin;\n"
"Texture2D<float2> yMax;\n"
"RWTexture3D<float4> state;\n"
"\n"
"[numthreads(NUM_X_THREADS, 1, NUM_Z_THREADS)]\n"
"void CSMain(uint3 c : SV_DispatchThreadID)\n"
"{\n"
"    uint3 dim;\n"
"    state.GetDimensions(dim.x, dim.y, dim.z);\n"
"    state[uint3(c.x, 0, c.z)] = float4(yMin[c.xz].x, 0.0f, yMin[c.xz].y, 0.0f);\n"
"    state[uint3(c.x, dim.y - 1, c.z)] = float4(yMax[c.xz].x, 0.0f, yMax[c.xz].y, 0.0f);\n"
"}\n"
"#endif\n"
"\n"
"#if USE_COPY_Z_FACE\n"
"Texture3D<float4> state;\n"
"RWTexture2D<float2> zMin;\n"
"RWTexture2D<float2> zMax;\n"
"\n"
"[numthreads(NUM_X_THREADS, NUM_Y_THREADS, 1)]\n"
"void CSMain(uint3 c : SV_DispatchThreadID)\n"
"{\n"
"    uint3 dim;\n"
"    state.GetDimensions(dim.x, dim.y, dim.z);\n"
"    zMin[c.xy] = state[uint3(c.x, c.y, 1)].xy;\n"
"    zMax[c.xy] = state[uint3(c.x, c.y, dim.z - 2)].xy;\n"
"}\n"
"#endif\n"
"\n"
"#if USE_WRITE_Z_FACE\n"
"Texture2D<float2> zMin;\n"
"Texture2D<float2> zMax;\n"
"RWTexture3D<float4> state;\n"
"\n"
"[numthreads(NUM_X_THREADS, NUM_Y_THREADS, 1)]\n"
"void CSMain(uint3 c : SV_DispatchThreadID)\n"
"{\n"
"    uint3 dim;\n"
"    state.GetDimensions(dim.x, dim.y, dim.z);\n"
"    state[uint3(c.x, c.y, 0)] = float4(zMin[c.xy].x, zMin[c.xy].y, 0.0f, 0.0f);\n"
"    state[uint3(c.x, c.y, dim.z - 1)] = float4(zMax[c.xy].x, zMax[c.xy].y, 0.0f, 0.0f);\n"
"}\n"
"#endif\n";

std::string const* Fluid3EnforceStateBoundary::msSource[] =
{
    &msGLSLSource,
    &msHLSLSource
};
