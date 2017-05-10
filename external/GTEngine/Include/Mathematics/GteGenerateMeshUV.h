// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <LowLevel/GteComputeModel.h>
#include <LowLevel/GteWrapper.h>
#include <Mathematics/GteVector2.h>
#include <Mathematics/GteVector3.h>
#include <Mathematics/GteETManifoldMesh.h>
#include <Mathematics/GteConstants.h>
#if defined(GTE_COMPUTE_MODEL_ALLOW_GPGPU)
#include <Graphics/GteConstantBuffer.h>
#include <Graphics/GteStructuredBuffer.h>
#include <Graphics/GteComputeProgram.h>
#endif
#include <algorithm>
#include <cstring>
#include <functional>
#include <limits>
#include <memory>
#include <set>
#include <thread>
#include <vector>

// This class is an implementation of the barycentric mapping algorithm
// described in Section 5.3 of the book
//     Polygon Mesh Processing
//     Mario Botsch, Leif Kobbelt, Mark Pauly, Pierre Alliez, Bruno Levy
//     AK Peters, Ltd., Natick MA, 2010
// It uses the mean value weights described in Section 5.3.1 to allow the mesh
// geometry to influence the texture coordinate generation, and it uses
// Gauss-Seidel iteration to solve the sparse linear system.  The authors'
// advice is that the Gauss-Seidel approach works well for at most about 5000
// vertices, presumably the convergence rate degrading as the number of
// vertices increases.
//
// The algorithm implemented here has an additional preprocessing step that
// computes a topological distance transform of the vertices.  The boundary
// texture coordinates are propagated inward by updating the vertices in
// topological distance order, leading to fast convergence for large numbers
// of vertices.

namespace gte
{

class GenerateMeshUVBase
{
protected:
    static std::string const msGLSLSource;
    static std::string const msHLSLSource;
    static std::string const* msSource[];
};

template <typename Real>
class GenerateMeshUV : public GenerateMeshUVBase
{
public:
    class UVComputeModel : public ComputeModel
    {
    public:
        virtual ~UVComputeModel();

        UVComputeModel();

        UVComputeModel(unsigned int inNumThreads,
            std::function<void(unsigned int)> const* inProgress);

#if defined(GTE_COMPUTE_MODEL_ALLOW_GPGPU)
        UVComputeModel(unsigned int inNumThreads,
            std::shared_ptr<GraphicsEngine> const& inEngine,
            std::shared_ptr<ProgramFactory> const& inFactory,
            std::function<void(unsigned int)> const* inProgress);
#endif

        std::function<void(unsigned int)> const* progress;
    };

    // Construction.
    GenerateMeshUV(std::shared_ptr<UVComputeModel> const& cmodel);

    // The incoming mesh must be edge-triangle manifold and have rectangle
    // topology (simply connected, closed polyline boundary).  The arrays
    // 'vertices' and 'tcoords' must both have 'numVertices' elements.  Set
    // 'useSquareTopology' to true for the generated coordinates to live
    // in the uv-square [0,1]^2.  Set it to false for the generated
    // coordinates to live in a convex polygon that inscribes the uv-disk
    // of center (1/2,1/2) and radius 1/2.
    void operator()(unsigned int numIterations, bool useSquareTopology,
        int numVertices, Vector3<Real> const* vertices, int numIndices,
        int const* indices, Vector2<Real>* tcoords);

private:
    void TopologicalVertexDistanceTransform();
    void AssignBoundaryTextureCoordinatesSquare();
    void AssignBoundaryTextureCoordinatesDisk();
    void ComputeMeanValueWeights();
    void SolveSystem(unsigned int numIterations);
    void SolveSystemCPUSingle(unsigned int numIterations);
    void SolveSystemCPUMultiple(unsigned int numIterations);

    // Convenience members that store the input parameters to operator().
    int mNumVertices;
    Vector3<Real> const* mVertices;
    Vector2<Real>* mTCoords;

    // The edge-triangle manifold graph, where each edge is shared by at most
    // two triangles.
    ETManifoldMesh mGraph;

    // The mVertexInfo array stores -1 for the interior vertices.  For a
    // boundary edge <v0,v1> that is counterclockwise, mVertexInfo[v0] = v1,
    // which gives us an orded boundary polyline.
    enum { INTERIOR_VERTEX = -1 };
    std::vector<int> mVertexInfo;
    int mNumBoundaryEdges, mBoundaryStart;
    typedef ETManifoldMesh::Edge Edge;
    std::set<std::shared_ptr<Edge>> mInteriorEdges;

    // The vertex graph required to set up a sparse linear system of equations
    // to determine the texture coordinates.
    struct Vertex
    {
        // The topological distance from the boundary of the mesh.
        int distance;

        // The value range[0] is the index into mVertexGraphData for the first
        // adjacent vertex.  The value range[1] is the number of adjacent
        // vertices.
        std::array<int, 2> range;

#if defined(GTE_COMPUTE_MODEL_ALLOW_GPGPU) && defined(GTE_DEV_OPENGL)
        int _padding;  // GLSL will map the ivec3 Vertex to ivec4 in an array
#endif
    };
    std::vector<Vertex> mVertexGraph;
    std::vector<std::pair<int, Real>> mVertexGraphData;

    // The vertices are listed in the order determined by a topological distance
    // transform.  Boundary vertices have 'distance' 0.  Any vertices that are
    // not boundary vertices but are edge-adjacent to boundary vertices have
    // 'distance' 1.  Neighbors of those have distance '2', and so on.  The
    // mOrderedVertices array stores distance-0 vertices first, distance-1
    // vertices second, and so on.
    std::vector<int> mOrderedVertices;

    std::shared_ptr<UVComputeModel> mCModel;

#if defined(GTE_COMPUTE_MODEL_ALLOW_GPGPU)
    // Support for solving the sparse linear system on the GPU.
    void SolveSystemGPU(unsigned int numIterations);
    std::shared_ptr<ComputeProgram> mSolveSystem;
    std::shared_ptr<ConstantBuffer> mBoundBuffer;
    std::shared_ptr<StructuredBuffer> mVGBuffer;
    std::shared_ptr<StructuredBuffer> mVGDBuffer;
    std::shared_ptr<StructuredBuffer> mOVBuffer;
    std::shared_ptr<StructuredBuffer> mTCoordsBuffer[2];
#endif
};


template <typename Real>
GenerateMeshUV<Real>::UVComputeModel::~UVComputeModel()
{
}

template <typename Real>
GenerateMeshUV<Real>::UVComputeModel::UVComputeModel()
    :
    progress(nullptr)
{
}

template <typename Real>
GenerateMeshUV<Real>::UVComputeModel::UVComputeModel(unsigned int inNumThreads,
    std::function<void(unsigned int)> const* inProgress)
    :
    ComputeModel(inNumThreads),
    progress(inProgress)
{
}

#if defined(GTE_COMPUTE_MODEL_ALLOW_GPGPU)

template <typename Real>
GenerateMeshUV<Real>::UVComputeModel::UVComputeModel(unsigned int inNumThreads,
    std::shared_ptr<GraphicsEngine> const& inEngine,
    std::shared_ptr<ProgramFactory> const& inFactory,
    std::function<void(unsigned int)> const* inProgress)
    :
    ComputeModel(inNumThreads, inEngine, inFactory),
    progress(inProgress)
{
}

#endif

template <typename Real>
GenerateMeshUV<Real>::GenerateMeshUV(std::shared_ptr<UVComputeModel> const& cmodel)
    :
    mCModel(cmodel),
    mNumVertices(0),
    mVertices(nullptr),
    mTCoords(nullptr),
    mNumBoundaryEdges(0),
    mBoundaryStart(0)
{
}

template <typename Real>
void GenerateMeshUV<Real>::operator()(unsigned int numIterations,
    bool useSquareTopology, int numVertices, Vector3<Real> const* vertices,
    int numIndices, int const* indices, Vector2<Real>* tcoords)
{
    // Ensure that numIterations is even, which avoids having a memory
    // copy from the temporary ping-pong buffer to 'tcoords'.
    if (numIterations & 1)
    {
        ++numIterations;
    }

    mNumVertices = numVertices;
    mVertices = vertices;
    mTCoords = tcoords;

    // The linear system solver has a first pass to initialize the texture
    // coordinates to ensure the Gauss-Seidel iteration converges rapidly.
    // This requires the texture coordinates all start as (-1,-1).
    for (int i = 0; i < numVertices; ++i)
    {
        mTCoords[i][0] = (Real)-1;
        mTCoords[i][1] = (Real)-1;
    }

    // Create the manifold mesh data structure.
    mGraph.Clear();
    int const numTriangles = numIndices / 3;
    for (int t = 0; t < numTriangles; ++t)
    {
        int v0 = *indices++;
        int v1 = *indices++;
        int v2 = *indices++;
        mGraph.Insert(v0, v1, v2);
    }

    TopologicalVertexDistanceTransform();

    if (useSquareTopology)
    {
        AssignBoundaryTextureCoordinatesSquare();
    }
    else
    {
        AssignBoundaryTextureCoordinatesDisk();
    }

    ComputeMeanValueWeights();
    SolveSystem(numIterations);
}

template <typename Real>
void GenerateMeshUV<Real>::TopologicalVertexDistanceTransform()
{
    // Initialize the graph information.
    mVertexInfo.resize(mNumVertices);
    std::fill(mVertexInfo.begin(), mVertexInfo.end(), INTERIOR_VERTEX);
    mVertexGraph.resize(mNumVertices);
    mVertexGraphData.resize(2 * mGraph.GetEdges().size());
    std::pair<int, Real> initialData = std::make_pair(-1, (Real)-1);
    std::fill(mVertexGraphData.begin(), mVertexGraphData.end(), initialData);
    mOrderedVertices.resize(mNumVertices);
    mInteriorEdges.clear();
    mNumBoundaryEdges = 0;
    mBoundaryStart = std::numeric_limits<int>::max();

    // Count the number of adjacent vertices for each vertex.  For data sets
    // with a large number of vertices, this is a preprocessing step to avoid
    // a dynamic data structure that has a large number of std:map objects
    // that take a very long time to destroy when a debugger is attached to
    // the executable.  Instead, we allocate a single array that stores all
    // the adjacency information.  It is also necessary to bundle the data
    // this way for a GPU version of the algorithm.
    std::vector<int> numAdjacencies(mNumVertices);
    std::fill(numAdjacencies.begin(), numAdjacencies.end(), 0);

    for (auto const& element : mGraph.GetEdges())
    {
        ++numAdjacencies[element.first.V[0]];
        ++numAdjacencies[element.first.V[1]];

        if (element.second->T[1].lock())
        {
            // This is an interior edge.
            mInteriorEdges.insert(element.second);
        }
        else
        {
            // This is a boundary edge.  Determine the ordering of the
            // vertex indices to make the edge counterclockwise.
            ++mNumBoundaryEdges;
            int v0 = element.second->V[0], v1 = element.second->V[1];
            auto tri = element.second->T[0].lock();
            int i;
            for (i = 0; i < 3; ++i)
            {
                int v2 = tri->V[i];
                if (v2 != v0 && v2 != v1)
                {
                    // The vertex is opposite the boundary edge.
                    v0 = tri->V[(i + 1) % 3];
                    v1 = tri->V[(i + 2) % 3];
                    mVertexInfo[v0] = v1;
                    mBoundaryStart = std::min(mBoundaryStart, v0);
                    break;
                }
            }
        }
    }

    // Set the range data for each vertex.
    for (int vIndex = 0, aIndex = 0; vIndex < mNumVertices; ++vIndex)
    {
        int numAdjacent = numAdjacencies[vIndex];
        mVertexGraph[vIndex].range = { { aIndex, numAdjacent } };
        aIndex += numAdjacent;

#if defined(GTE_COMPUTE_MODEL_ALLOW_GPGPU) && defined(GTE_DEV_OPENGL)
        // Initialize the padding, even though it is unused.
        mVertexGraph[vIndex]._padding = 0;
#endif
    }

    // Compute a topological distance transform of the vertices.
    std::set<int> currFront;
    for (auto const& element : mGraph.GetEdges())
    {
        int v0 = element.second->V[0], v1 = element.second->V[1];
        for (int i = 0; i < 2; ++i)
        {
            if (mVertexInfo[v0] == INTERIOR_VERTEX)
            {
                mVertexGraph[v0].distance = -1;
            }
            else
            {
                mVertexGraph[v0].distance = 0;
                currFront.insert(v0);
            }

            // Insert v1 into the first available slot of the adjacency array.
            std::array<int, 2> range = mVertexGraph[v0].range;
            for (int j = 0; j < range[1]; ++j)
            {
                std::pair<int, Real>& data = mVertexGraphData[range[0] + j];
                if (data.second == (Real)-1)
                {
                    data.first = v1;
                    data.second = (Real)0;
                    break;
                }
            }

            std::swap(v0, v1);
        }
    }

    // Use a breadth-first search to propagate the distance information.
    int currDistance = 0, nextDistance = 1;
    size_t numFrontVertices = currFront.size();
    std::copy(currFront.begin(), currFront.end(), mOrderedVertices.begin());
    while (currFront.size() > 0)
    {
        std::set<int> nextFront;
        for (auto v : currFront)
        {
            std::array<int, 2> range = mVertexGraph[v].range;
            auto* current = &mVertexGraphData[range[0]];
            for (int j = 0; j < range[1]; ++j, ++current)
            {
                int a = current->first;
                if (mVertexGraph[a].distance == -1)
                {
                    mVertexGraph[a].distance = nextDistance;
                    nextFront.insert(a);
                }
            }
        }
        std::copy(nextFront.begin(), nextFront.end(), mOrderedVertices.begin() + numFrontVertices);
        numFrontVertices += nextFront.size();
        currFront = std::move(nextFront);
        currDistance = nextDistance++;
    }
}

template <typename Real>
void GenerateMeshUV<Real>::AssignBoundaryTextureCoordinatesSquare()
{
    // Map the boundary of the mesh to the unit square [0,1]^2.  The selection
    // of square vertices is such that the relative distances between boundary
    // vertices and the relative distances between polygon vertices is
    // preserved, except that the four corners of the square are required to
    // have boundary points mapped to them.  The first boundary point has an
    // implied distance of zero.  The value distance[i] is the length of the
    // boundary polyline from vertex 0 to vertex i+1.
    std::vector<Real> distance(mNumBoundaryEdges);
    Real total = (Real)0;
    int v0 = mBoundaryStart, v1, i;
    for (i = 0; i < mNumBoundaryEdges; ++i)
    {
        v1 = mVertexInfo[v0];
        total += Length(mVertices[v1] - mVertices[v0]);
        distance[i] = total;
        v0 = v1;
    }

    Real invTotal = ((Real)1) / total;
    for (auto& d : distance)
    {
        d *= invTotal;
    }

    auto begin = distance.begin(), end = distance.end();
    int endYMin = (int)(std::lower_bound(begin, end, (Real)0.25) - begin);
    int endXMax = (int)(std::lower_bound(begin, end, (Real)0.50) - begin);
    int endYMax = (int)(std::lower_bound(begin, end, (Real)0.75) - begin);
    int endXMin = (int)distance.size() - 1;

    // The first polygon vertex is (0,0).  The remaining vertices are chosen
    // counterclockwise around the square.
    v0 = mBoundaryStart;
    mTCoords[v0][0] = (Real)0;
    mTCoords[v0][1] = (Real)0;
    for (i = 0; i < endYMin; ++i)
    {
        v1 = mVertexInfo[v0];
        mTCoords[v1][0] = distance[i] * (Real)4;
        mTCoords[v1][1] = (Real)0;
        v0 = v1;
    }

    v1 = mVertexInfo[v0];
    mTCoords[v1][0] = (Real)1;
    mTCoords[v1][1] = (Real)0;
    v0 = v1;
    for (++i; i < endXMax; ++i)
    {
        v1 = mVertexInfo[v0];
        mTCoords[v1][0] = (Real)1;
        mTCoords[v1][1] = distance[i] * (Real)4 - (Real)1;
        v0 = v1;
    }

    v1 = mVertexInfo[v0];
    mTCoords[v1][0] = (Real)1;
    mTCoords[v1][1] = (Real)1;
    v0 = v1;
    for (++i; i < endYMax; ++i)
    {
        v1 = mVertexInfo[v0];
        mTCoords[v1][0] = (Real)3 - distance[i] * (Real)4;
        mTCoords[v1][1] = (Real)1;
        v0 = v1;
    }

    v1 = mVertexInfo[v0];
    mTCoords[v1][0] = (Real)0;
    mTCoords[v1][1] = (Real)1;
    v0 = v1;
    for (++i; i < endXMin; ++i)
    {
        v1 = mVertexInfo[v0];
        mTCoords[v1][0] = (Real)0;
        mTCoords[v1][1] = (Real)4 - distance[i] * (Real)4;
        v0 = v1;
    }
}

template <typename Real>
void GenerateMeshUV<Real>::AssignBoundaryTextureCoordinatesDisk()
{
    // Map the boundary of the mesh to a convex polygon.  The selection of
    // convex polygon vertices is such that the relative distances between
    // boundary vertices and the relative distances between polygon vertices
    // is preserved.  The first boundary point has an implied distance of
    // zero.  The value distance[i] is the length of the boundary polyline
    // from vertex 0 to vertex i+1.
    std::vector<Real> distance(mNumBoundaryEdges);
    Real total = (Real)0;
    int v0 = mBoundaryStart;
    for (int i = 0; i < mNumBoundaryEdges; ++i)
    {
        int v1 = mVertexInfo[v0];
        total += Length(mVertices[v1] - mVertices[v0]);
        distance[i] = total;
        v0 = v1;
    }

    // The convex polygon lives in [0,1]^2 and inscribes a circle with center
    // (1/2,1/2) and radius 1/2.  The polygon center is not necessarily the
    // circle center!  This is the case when a boundary edge has length larger
    // than half the total length of the boundary polyline; we do not expect
    // such data for our meshes.  The first polygon vertex is (1/2,0).  The
    // remaining vertices are chosen counterclockwise around the polygon.
    Real multiplier = ((Real)GTE_C_TWO_PI) / total;
    v0 = mBoundaryStart;
    mTCoords[v0][0] = (Real)1;
    mTCoords[v0][1] = (Real)0.5;
    for (int i = 1; i < mNumBoundaryEdges; ++i)
    {
        int v1 = mVertexInfo[v0];
        Real angle = multiplier * distance[i - 1];
        mTCoords[v1][0] = (cos(angle) + (Real)1) * (Real)0.5;
        mTCoords[v1][1] = (sin(angle) + (Real)1) * (Real)0.5;
        v0 = v1;
    }
}

template <typename Real>
void GenerateMeshUV<Real>::ComputeMeanValueWeights()
{
    for (auto const& edge : mInteriorEdges)
    {
        int v0 = edge->V[0], v1 = edge->V[1];
        for (int i = 0; i < 2; ++i)
        {
            // Compute the direction from X0 to X1 and compute the length
            // of the edge (X0,X1).
            Vector3<Real> X0 = mVertices[v0];
            Vector3<Real> X1 = mVertices[v1];
            Vector3<Real> X1mX0 = X1 - X0;
            Real x1mx0length = Normalize(X1mX0);
            Real weight;
            if (x1mx0length >(Real)0)
            {
                // Compute the weight for X0 associated with X1.
                weight = (Real)0;
                for (int j = 0; j < 2; ++j)
                {
                    // Find the vertex of triangle T[j] opposite edge <X0,X1>.
                    auto tri = edge->T[j].lock();
                    int k;
                    for (k = 0; k < 3; ++k)
                    {
                        int v2 = tri->V[k];
                        if (v2 != v0 && v2 != v1)
                        {
                            Vector3<Real> X2 = mVertices[v2];
                            Vector3<Real> X2mX0 = X2 - X0;
                            Real x2mx0Length = Normalize(X2mX0);
                            if (x2mx0Length >(Real)0)
                            {
                                Real dot = Dot(X2mX0, X1mX0);
                                Real cs = std::min(std::max(dot, (Real)-1), (Real)1);
                                Real angle = acos(cs);
                                weight += tan(angle * (Real)0.5);
                            }
                            else
                            {
                                weight += (Real)1;
                            }
                            break;
                        }
                    }
                }
                weight /= x1mx0length;
            }
            else
            {
                weight = (Real)1;
            }

            std::array<int, 2> range = mVertexGraph[v0].range;
            for (int j = 0; j < range[1]; ++j)
            {
                std::pair<int, Real>& data = mVertexGraphData[range[0] + j];
                if (data.first == v1)
                {
                    data.second = weight;
                }
            }

            std::swap(v0, v1);
        }
    }
}

template <typename Real>
void GenerateMeshUV<Real>::SolveSystem(unsigned int numIterations)
{
    // On the first pass, average only neighbors whose texture coordinates
    // have been computed.  This is a good initial guess for the linear system
    // and leads to relatively fast convergence of the Gauss-Seidel iterates.
    Real zero = (Real)0;
    for (int i = mNumBoundaryEdges; i < mNumVertices; ++i)
    {
        int v0 = mOrderedVertices[i];
        std::array<int, 2> range = mVertexGraph[v0].range;
        auto const* current = &mVertexGraphData[range[0]];
        Vector2<Real> tcoord{ zero, zero };
        Real weight, weightSum = zero;
        for (int j = 0; j < range[1]; ++j, ++current)
        {
            int v1 = current->first;
            if (mTCoords[v1][0] != -1.0f)
            {
                weight = current->second;
                weightSum += weight;
                tcoord += weight * mTCoords[v1];
            }
        }
        tcoord /= weightSum;
        mTCoords[v0] = tcoord;
    }

#if defined(GTE_COMPUTE_MODEL_ALLOW_GPGPU)
    if (mCModel->engine)
    {
        SolveSystemGPU(numIterations);
    }
    else
#endif
    {
        if (mCModel->numThreads > 1)
        {
            SolveSystemCPUMultiple(numIterations);
        }
        else
        {
            SolveSystemCPUSingle(numIterations);
        }
    }
}

template <typename Real>
void GenerateMeshUV<Real>::SolveSystemCPUSingle(unsigned int numIterations)
{
    // Use ping-pong buffers for the texture coordinates.
    std::vector<Vector2<Real>> tcoords(mNumVertices);
    size_t numBytes = mNumVertices * sizeof(Vector2<Real>);
    Memcpy(&tcoords[0], mTCoords, numBytes);
    Vector2<Real>* inTCoords = mTCoords;
    Vector2<Real>* outTCoords = &tcoords[0];

    // The value numIterations is even, so we always swap an even number
    // of times.  This ensures that on exit from the loop, outTCoords is
    // tcoords.
    for (unsigned int i = 1; i <= numIterations; ++i)
    {
        if (mCModel->progress)
        {
            (*mCModel->progress)(i);
        }

        for (int j = mNumBoundaryEdges; j < mNumVertices; ++j)
        {
            int v0 = mOrderedVertices[j];
            std::array<int, 2> range = mVertexGraph[v0].range;
            auto const* current = &mVertexGraphData[range[0]];
            Vector2<Real> tcoord{ (Real)0, (Real)0 };
            Real weight, weightSum = (Real)0;
            for (int k = 0; k < range[1]; ++k, ++current)
            {
                int v1 = current->first;
                weight = current->second;
                weightSum += weight;
                tcoord += weight * inTCoords[v1];
            }
            tcoord /= weightSum;
            outTCoords[v0] = tcoord;
        }

        std::swap(inTCoords, outTCoords);
    }
}

template <typename Real>
void GenerateMeshUV<Real>::SolveSystemCPUMultiple(unsigned int numIterations)
{
    // Use ping-pong buffers for the texture coordinates.
    std::vector<Vector2<Real>> tcoords(mNumVertices);
    size_t numBytes = mNumVertices * sizeof(Vector2<Real>);
    Memcpy(&tcoords[0], mTCoords, numBytes);
    Vector2<Real>* inTCoords = mTCoords;
    Vector2<Real>* outTCoords = &tcoords[0];

    // Partition the data for multiple threads.
    int numV = mNumVertices - mNumBoundaryEdges;
    int numVPerThread = numV / mCModel->numThreads;
    std::vector<int> vmin(mCModel->numThreads), vmax(mCModel->numThreads);
    for (unsigned int t = 0; t < mCModel->numThreads; ++t)
    {
        vmin[t] = mNumBoundaryEdges + t * numVPerThread;
        vmax[t] = vmin[t] + numVPerThread - 1;
    }
    vmax[mCModel->numThreads - 1] = mNumVertices - 1;

    // The value numIterations is even, so we always swap an even number
    // of times.  This ensures that on exit from the loop, outTCoords is
    // tcoords.
    for (unsigned int i = 1; i <= numIterations; ++i)
    {
        if (mCModel->progress)
        {
            (*mCModel->progress)(i);
        }

        // Execute Gauss-Seidel iterations in multiple threads.
        std::vector<std::thread> process(mCModel->numThreads);
        for (unsigned int t = 0; t < mCModel->numThreads; ++t)
        {
            process[t] = std::thread([this, t, &vmin, &vmax, inTCoords,
                outTCoords]()
            {
                for (int j = vmin[t]; j <= vmax[t]; ++j)
                {
                    int v0 = mOrderedVertices[j];
                    std::array<int, 2> range = mVertexGraph[v0].range;
                    auto const* current = &mVertexGraphData[range[0]];
                    Vector2<Real> tcoord{ (Real)0, (Real)0 };
                    Real weight, weightSum = (Real)0;
                    for (int k = 0; k < range[1]; ++k, ++current)
                    {
                        int v1 = current->first;
                        weight = current->second;
                        weightSum += weight;
                        tcoord += weight * inTCoords[v1];
                    }
                    tcoord /= weightSum;
                    outTCoords[v0] = tcoord;
                }
            });
        }

        // Wait for all threads to finish.
        for (unsigned int t = 0; t < mCModel->numThreads; ++t)
        {
            process[t].join();
        }

        std::swap(inTCoords, outTCoords);
    }
}

#if defined(GTE_COMPUTE_MODEL_ALLOW_GPGPU)

template <typename Real>
void GenerateMeshUV<Real>::SolveSystemGPU(unsigned int numIterations)
{
    mCModel->factory->defines.Set("NUM_X_THREADS", 8);
    mCModel->factory->defines.Set("NUM_Y_THREADS", 8);
    if (std::numeric_limits<Real>::max() == std::numeric_limits<float>::max())
    {
        mCModel->factory->defines.Set("Real", "float");
#if defined(GTE_DEV_OPENGL)
        mCModel->factory->defines.Set("Real2", "vec2");
#else
        mCModel->factory->defines.Set("Real2", "float2");
#endif
    }
    else
    {
        mCModel->factory->defines.Set("Real", "double");
#if defined(GTE_DEV_OPENGL)
        mCModel->factory->defines.Set("Real2", "dvec2");
#else
        mCModel->factory->defines.Set("Real2", "double2");
#endif
    }

    // TODO: Test mSolveSystem for null and respond accordingly.
    int api = mCModel->factory->GetAPI();
    mSolveSystem = mCModel->factory->CreateFromSource(*msSource[api]);
    std::shared_ptr<ComputeShader> cshader = mSolveSystem->GetCShader();

    // Compute the number of thread groups.
    int numInputs = mNumVertices - mNumBoundaryEdges;
    Real factor0 = ceil(sqrt((Real)numInputs));
    Real factor1 = ceil((Real)numInputs / factor0);
    int xElements = static_cast<int>(factor0);
    int yElements = static_cast<int>(factor1);
    int xRem = (xElements % 8);
    if (xRem > 0)
    {
        xElements += 8 - xRem;
    }
    int yRem = (yElements % 8);
    if (yRem > 0)
    {
        yElements += 8 - yRem;
    }
    unsigned int numXGroups = xElements / 8;
    unsigned int numYGroups = yElements / 8;

    mBoundBuffer = std::make_shared<ConstantBuffer>(4 * sizeof(int), false);
    int* data = mBoundBuffer->Get<int>();
    data[0] = xElements;
    data[1] = yElements;
    data[2] = mNumBoundaryEdges;
    data[3] = numInputs;
    cshader->Set("Bounds", mBoundBuffer);

    unsigned int const vgSize = static_cast<unsigned int>(mVertexGraph.size());
    mVGBuffer = std::make_shared<StructuredBuffer>(vgSize, sizeof(Vertex));
    Memcpy(mVGBuffer->GetData(), &mVertexGraph[0], mVGBuffer->GetNumBytes());
    cshader->Set("vertexGraph", mVGBuffer);

    unsigned int const vgdSize = static_cast<unsigned int>(mVertexGraphData.size());
    mVGDBuffer = std::make_shared<StructuredBuffer>(vgdSize, sizeof(std::pair<int, Real>));
    Memcpy(mVGDBuffer->GetData(), &mVertexGraphData[0], mVGDBuffer->GetNumBytes());
    cshader->Set("vertexGraphData", mVGDBuffer);

    unsigned int const ovSize = static_cast<unsigned int>(mOrderedVertices.size());
    mOVBuffer = std::make_shared<StructuredBuffer>(ovSize, sizeof(int));
    Memcpy(mOVBuffer->GetData(), &mOrderedVertices[0], mOVBuffer->GetNumBytes());
    cshader->Set("orderedVertices", mOVBuffer);

    for (int j = 0; j < 2; ++j)
    {
        mTCoordsBuffer[j] = std::make_shared<StructuredBuffer>(mNumVertices, sizeof(Vector2<Real>));
        mTCoordsBuffer[j]->SetUsage(Resource::SHADER_OUTPUT);
        Memcpy(mTCoordsBuffer[j]->GetData(), mTCoords, mTCoordsBuffer[j]->GetNumBytes());
    }
    mTCoordsBuffer[0]->SetCopyType(Resource::COPY_STAGING_TO_CPU);

    // The value numIterations is even, so we always swap an even number
    // of times.  This ensures that on exit from the loop,
    // mTCoordsBuffer[0] has the final output.
    for (unsigned int i = 1; i <= numIterations; ++i)
    {
        if (mCModel->progress)
        {
            (*mCModel->progress)(i);
        }

        cshader->Set("inTCoords", mTCoordsBuffer[0]);
        cshader->Set("outTCoords", mTCoordsBuffer[1]);
        mCModel->engine->Execute(mSolveSystem, numXGroups, numYGroups, 1);
        std::swap(mTCoordsBuffer[0], mTCoordsBuffer[1]);
    }

    mCModel->engine->CopyGpuToCpu(mTCoordsBuffer[0]);
    Memcpy(mTCoords, mTCoordsBuffer[0]->GetData(), mTCoordsBuffer[0]->GetNumBytes());
}

#endif

}
