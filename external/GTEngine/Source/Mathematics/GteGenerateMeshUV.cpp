// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#include <GTEnginePCH.h>
#include <Mathematics/GteGenerateMeshUV.h>
using namespace gte;

std::string const GenerateMeshUVBase::msGLSLSource =
"uniform Bounds\n"
"{\n"
"    ivec2 bound;\n"
"    int numBoundaryEdges;\n"
"    int numInputs;\n"
"};\n"
"\n"
"struct VertexGraphData\n"
"{\n"
"    int adjacent;\n"
"    Real weight;\n"
"};\n"
"\n"
"buffer vertexGraph { ivec3 data[]; } vertexGraphSB;\n"
"buffer vertexGraphData { VertexGraphData data[]; } vertexGraphDataSB;\n"
"buffer orderedVertices { int data[]; } orderedVerticesSB;\n"
"buffer inTCoords { Real2 data[]; } inTCoordsSB;\n"
"buffer outTCoords { Real2 data[]; } outTCoordsSB;\n"
"\n"
"layout (local_size_x = NUM_X_THREADS, local_size_y = NUM_Y_THREADS, local_size_z = 1) in;\n"
"void main()\n"
"{\n"
"    ivec2 t = ivec2(gl_GlobalInvocationID.xy);\n"
"    int index = t.x + bound.x * t.y;\n"
"    if (step(index, numInputs-1) == 1)\n"
"    {\n"
"        int v = orderedVerticesSB.data[numBoundaryEdges + index];\n"
"        ivec2 range = vertexGraphSB.data[v].yz;\n"
"        Real2 tcoord = Real2(0, 0);\n"
"        Real weightSum = 0;\n"
"        for (int j = 0; j < range.y; ++j)\n"
"        {\n"
"            VertexGraphData vgd = vertexGraphDataSB.data[range.x + j];\n"
"            weightSum += vgd.weight;\n"
"            tcoord += vgd.weight * inTCoordsSB.data[vgd.adjacent];\n"
"        }\n"
"        tcoord /= weightSum;\n"
"        outTCoordsSB.data[v] = tcoord;\n"
"    }\n"
"}\n";

std::string const GenerateMeshUVBase::msHLSLSource =
"cbuffer Bounds\n"
"{\n"
"    int2 bound;\n"
"    int numBoundaryEdges;\n"
"    int numInputs;\n"
"};\n"
"\n"
"struct VertexGraphData\n"
"{\n"
"    int adjacent;\n"
"    Real weight;\n"
"};\n"
"\n"
"StructuredBuffer<int3> vertexGraph;\n"
"StructuredBuffer<VertexGraphData> vertexGraphData;\n"
"StructuredBuffer<int> orderedVertices;\n"
"StructuredBuffer<Real2> inTCoords;\n"
"RWStructuredBuffer<Real2> outTCoords;\n"
"\n"
"[numthreads(NUM_X_THREADS, NUM_Y_THREADS, 1)]\n"
"void CSMain(int2 t : SV_DispatchThreadID)\n"
"{\n"
"    int index = t.x + bound.x * t.y;\n"
"    if (step(index, numInputs-1))\n"
"    {\n"
"        int v = orderedVertices[numBoundaryEdges + index];\n"
"        int2 range = vertexGraph[v].yz;\n"
"        Real2 tcoord = Real2(0, 0);\n"
"        Real weightSum = 0;\n"
"        for (int j = 0; j < range.y; ++j)\n"
"        {\n"
"            VertexGraphData data = vertexGraphData[range.x + j];\n"
"            weightSum += data.weight;\n"
"            tcoord += data.weight * inTCoords[data.adjacent];\n"
"        }\n"
"        tcoord /= weightSum;\n"
"        outTCoords[v] = tcoord;\n"
"    }\n"
"}\n";

std::string const* GenerateMeshUVBase::msSource[] =
{
    &msGLSLSource,
    &msHLSLSource
};
