// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.8 (2015/11/21)

#include "Wm5MathematicsPCH.h"
#include "Wm5Delaunay2.h"
#include "Wm5Query2Filtered.h"
#include "Wm5Query2Int64.h"
#include "Wm5Query2Integer.h"
#include "Wm5Query2Rational.h"

namespace Wm5
{

template <typename Real>
const int Delaunay2<Real>::msIndex[3][2] = {{0,1},{1,2},{2,0}};

//----------------------------------------------------------------------------
template <typename Real>
Delaunay2<Real>::Delaunay2 (int numVertices, Vector2<Real>* vertices,
    Real epsilon, bool owner, Query::Type queryType)
    :
    Delaunay<Real>(numVertices, epsilon, owner, queryType),
    mVertices(vertices),
    mNumUniqueVertices(0),
    mSVertices(0),
    mQuery(0),
    mLineOrigin(Vector2<Real>::ZERO),
    mLineDirection(Vector2<Real>::ZERO),
    mPathLast(-1),
    mPath(0),
    mLastEdgeV0(-1),
    mLastEdgeV1(-1),
    mLastEdgeOpposite(-1),
    mLastEdgeOppositeIndex(-1)
{
    typename Vector2<Real>::Information info;
    Vector2<Real>::GetInformation(mNumVertices, mVertices, mEpsilon, info);
    if (info.mDimension == 0)
    {
        // The values of mDimension, mIndices, and mAdjacencies were
        // already initialized by the Delaunay base class.
        return;
    }

    if (info.mDimension == 1)
    {
        // The set is (nearly) collinear.  The caller is responsible for
        // creating a Delaunay1 object.
        mDimension = 1;
        mLineOrigin = info.mOrigin;
        mLineDirection = info.mDirection[0];
        return;
    }

    mDimension = 2;

    // Allocate storage for the input vertices and the supertriangle
    // vertices.
    mSVertices = new1<Vector2<Real> >(mNumVertices);
    int i;

    if (queryType != Query::QT_RATIONAL && queryType != Query::QT_FILTERED)
    {
        // Transform the vertices to the square [0,1]^2.
        mMin = Vector2<Real>(info.mMin[0], info.mMin[1]);
        mScale = ((Real)1)/info.mMaxRange;
        for (i = 0; i < mNumVertices; ++i)
        {
            mSVertices[i] = (mVertices[i] - mMin)*mScale;
        }

        Real expand;
        if (queryType == Query::QT_INT64)
        {
            // Scale the vertices to the square [0,2^{16}]^2 to allow use of
            // 64-bit integers for triangulation.
            expand = (Real)(1 << 16);
            mQuery = new0 Query2Int64<Real>(mNumVertices, mSVertices);
        }
        else if (queryType == Query::QT_INTEGER)
        {
            // Scale the vertices to the square [0,2^{20}]^2 to get more
            // precision for TInteger than for 64-bit integers for
            // triangulation.
            expand = (Real)(1 << 20);
            mQuery = new0 Query2Integer<Real>(mNumVertices, mSVertices);
        }
        else // queryType == Query::QT_REAL
        {
            // No scaling for floating point.
            expand = (Real)1;
            mQuery = new0 Query2<Real>(mNumVertices, mSVertices);
        }

        mScale *= expand;
        for (i = 0; i < mNumVertices; ++i)
        {
            mSVertices[i] *= expand;
        }
    }
    else
    {
        // No transformation needed for exact rational arithmetic.
        mMin = Vector2<Real>::ZERO;
        mScale = (Real)1;
        memcpy(mSVertices, mVertices, mNumVertices*sizeof(Vector2<Real>));

        if (queryType == Query::QT_RATIONAL)
        {
            mQuery = new0 Query2Rational<Real>(mNumVertices, mSVertices);
        }
        else // queryType == Query::QT_FILTERED
        {
            mQuery = new0 Query2Filtered<Real>(mNumVertices, mSVertices,
                mEpsilon);
        }
    }

    // Insert the (nondegenerate) triangle constructed by the call to
    // GetInformation. This is necessary for the circumcircle-visibility
    // algorithm to work correctly.
    if (!info.mExtremeCCW)
    {
        std::swap(info.mExtreme[1], info.mExtreme[2]);
    }
    mTriMesh.InsertTriangle(info.mExtreme[0], info.mExtreme[1],
        info.mExtreme[2]);

    // Incrementally update the triangulation.  The set of processed points
    // is maintained to eliminate duplicates, either in the original input
    // points or in the points obtained by snap rounding.
    std::set<Vector2<Real> > processed;
    for (i = 0; i < 3; ++i)
    {
        processed.insert(mSVertices[info.mExtreme[i]]);
    }
    for (i = 0; i < mNumVertices; ++i)
    {
        if (processed.find(mSVertices[i]) == processed.end())
        {
            Update(i);
            processed.insert(mSVertices[i]);
        }
    }
    mNumUniqueVertices = (int)processed.size();

    // Assign integer values to the triangles for use by the caller.
    std::map<Triangle*,int> permute;
    i = -1;
    permute[(Triangle*)0] = i++;
    const ETManifoldMesh::TMap& tmap = mTriMesh.GetTriangles();
    ETManifoldMesh::TMapCIterator element;
    for (element = tmap.begin(); element != tmap.end(); ++element)
    {
        permute[element->second] = i++;
    }

    // Put Delaunay triangles into an array (vertices and adjacency info).
    mNumSimplices = (int)mTriMesh.GetTriangles().size();
    if (mNumSimplices > 0)
    {
        mIndices = new1<int>(3*mNumSimplices);
        mAdjacencies = new1<int>(3*mNumSimplices);
        i = 0;
        for (element = tmap.begin(); element != tmap.end(); ++element)
        {
            const ETManifoldMesh::Triangle* tri = element->second;
            for (int j = 0; j < 3; ++j, ++i)
            {
                mIndices[i] = tri->V[j];
                mAdjacencies[i] = permute[tri->T[j]];
            }
        }
        assertion(i == 3*mNumSimplices, "Unexpected mismatch\n");

        mPathLast = -1;
        mPath = new1<int>(mNumSimplices + 1);
        memset(mPath, 0, (mNumSimplices + 1)*sizeof(int));
    }
}
//----------------------------------------------------------------------------
template <typename Real>
Delaunay2<Real>::~Delaunay2 ()
{
    delete0(mQuery);
    delete1(mSVertices);
    delete1(mPath);
    if (mOwner)
    {
        delete1(mVertices);
    }
}
//----------------------------------------------------------------------------
template <typename Real>
const Vector2<Real>* Delaunay2<Real>::GetVertices () const
{
    return mVertices;
}
//----------------------------------------------------------------------------
template <typename Real>
int Delaunay2<Real>::GetNumUniqueVertices () const
{
    return mNumUniqueVertices;
}
//----------------------------------------------------------------------------
template <typename Real>
const Vector2<Real>& Delaunay2<Real>::GetLineOrigin () const
{
    return mLineOrigin;
}
//----------------------------------------------------------------------------
template <typename Real>
const Vector2<Real>& Delaunay2<Real>::GetLineDirection () const
{
    return mLineDirection;
}
//----------------------------------------------------------------------------
template <typename Real>
Delaunay1<Real>* Delaunay2<Real>::GetDelaunay1 () const
{
    assertion(mDimension == 1, "The dimension must be 1\n");
    if (mDimension != 1)
    {
        return 0;
    }

    Real* projection = new1<Real>(mNumVertices);
    for (int i = 0; i < mNumVertices; ++i)
    {
        Vector2<Real> diff = mVertices[i] - mLineOrigin;
        projection[i] = mLineDirection.Dot(diff);
    }

    return new0 Delaunay1<Real>(mNumVertices, projection, mEpsilon, true,
        mQueryType);
}
//----------------------------------------------------------------------------
template <typename Real>
bool Delaunay2<Real>::GetHull (int& numEdges, int*& indices)
{
    assertion(mDimension == 2, "The dimension must be 2\n");
    if (mDimension != 2)
    {
        return false;
    }

    numEdges = 0;
    indices = 0;

    // Count the number of edges that are not shared by two triangles.
    int i, numAdjacent = 3*mNumSimplices;
    for (i = 0; i < numAdjacent; ++i)
    {
        if (mAdjacencies[i] == -1)
        {
            numEdges++;
        }
    }
    assertion(numEdges > 0, "There must be at least one triangle\n");
    if (numEdges == 0)
    {
        return false;
    }

    // Enumerate the edges.
    indices = new1<int>(2*numEdges);
    int* currentIndex = indices;
    for (i = 0; i < numAdjacent; ++i)
    {
        if (mAdjacencies[i] == -1)
        {
            int tri = i/3, j = i%3;
            *currentIndex++ = mIndices[3*tri + j];
            *currentIndex++ = mIndices[3*tri + ((j+1)%3)];
        }
    }

    return true;
}
//----------------------------------------------------------------------------
template <typename Real>
int Delaunay2<Real>::GetContainingTriangle (const Vector2<Real>& p) const
{
    assertion(mDimension == 2, "The dimension must be 2\n");
    if (mDimension != 2)
    {
        return -1;
    }

    // Convert to scaled coordinates.
    Vector2<Real> scP = (p - mMin)*mScale;

    // Start at first triangle in mesh.
    int index = (mPathLast >= 0 ? mPath[mPathLast] : 0);
    mPathLast = -1;
    mLastEdgeV0 = -1;
    mLastEdgeV1 = -1;
    mLastEdgeOpposite = -1;
    mLastEdgeOppositeIndex = -1;

    // Use triangle edges as binary separating lines.
    for (int i = 0; i < mNumSimplices; ++i)
    {
        mPath[++mPathLast] = index;

        int* vertices = &mIndices[3*index];

        if (mQuery->ToLine(scP, vertices[0], vertices[1]) > 0)
        {
            index = mAdjacencies[3*index];
            if (index == -1)
            {
                mLastEdgeV0 = vertices[0];
                mLastEdgeV1 = vertices[1];
                mLastEdgeOpposite = vertices[2];
                mLastEdgeOppositeIndex = 2;
                return -1;
            }
            continue;
        }

        if (mQuery->ToLine(scP, vertices[1], vertices[2]) > 0)
        {
            index = mAdjacencies[3*index + 1];
            if (index == -1)
            {
                mLastEdgeV0 = vertices[1];
                mLastEdgeV1 = vertices[2];
                mLastEdgeOpposite = vertices[0];
                mLastEdgeOppositeIndex = 0;
                return -1;
            }
            continue;
        }

        if (mQuery->ToLine(scP, vertices[2], vertices[0]) > 0)
        {
            index = mAdjacencies[3*index + 2];
            if (index == -1)
            {
                mLastEdgeV0 = vertices[2];
                mLastEdgeV1 = vertices[0];
                mLastEdgeOpposite = vertices[1];
                mLastEdgeOppositeIndex = 1;
                return -1;
            }
            continue;
        }

        mLastEdgeV0 = -1;
        mLastEdgeV1 = -1;
        mLastEdgeOpposite = -1;
        mLastEdgeOppositeIndex = -1;
        return index;
    }

    return -1;
}
//----------------------------------------------------------------------------
template <typename Real>
int Delaunay2<Real>::GetPathLast () const
{
    return mPathLast;
}
//----------------------------------------------------------------------------
template <typename Real>
const int* Delaunay2<Real>::GetPath () const
{
    return mPath;
}
//----------------------------------------------------------------------------
template <typename Real>
int Delaunay2<Real>::GetLastEdge (int& v0, int& v1, int& v2) const
{
    v0 = mLastEdgeV0;
    v1 = mLastEdgeV1;
    v2 = mLastEdgeOpposite;
    return mLastEdgeOppositeIndex;
}
//----------------------------------------------------------------------------
template <typename Real>
bool Delaunay2<Real>::GetVertexSet (int i, Vector2<Real> vertices[3]) const
{
    assertion(mDimension == 2, "The dimension must be 2\n");
    if (mDimension != 2)
    {
        return false;
    }

    if (0 <= i && i < mNumSimplices)
    {
        vertices[0] = mVertices[mIndices[3*i  ]];
        vertices[1] = mVertices[mIndices[3*i + 1]];
        vertices[2] = mVertices[mIndices[3*i + 2]];
        return true;
    }

    return false;
}
//----------------------------------------------------------------------------
template <typename Real>
bool Delaunay2<Real>::GetIndexSet (int i, int indices[3]) const
{
    assertion(mDimension == 2, "The dimension must be 2\n");
    if (mDimension != 2)
    {
        return false;
    }

    if (0 <= i && i < mNumSimplices)
    {
        indices[0] = mIndices[3*i  ];
        indices[1] = mIndices[3*i + 1];
        indices[2] = mIndices[3*i + 2];
        return true;
    }

    return false;
}
//----------------------------------------------------------------------------
template <typename Real>
bool Delaunay2<Real>::GetAdjacentSet (int i, int adjacencies[3]) const
{
    assertion(mDimension == 2, "The dimension must be 2\n");
    if (mDimension != 2)
    {
        return false;
    }

    if (0 <= i && i < mNumSimplices)
    {
        adjacencies[0] = mAdjacencies[3*i  ];
        adjacencies[1] = mAdjacencies[3*i + 1];
        adjacencies[2] = mAdjacencies[3*i + 2];
        return true;
    }

    return false;
}
//----------------------------------------------------------------------------
template <typename Real>
bool Delaunay2<Real>::GetBarycentricSet (int i, const Vector2<Real>& p,
    Real bary[3]) const
{
    assertion(mDimension == 2, "The dimension must be 2\n");
    if (mDimension != 2)
    {
        return false;
    }

    if (0 <= i && i < mNumSimplices)
    {
        Vector2<Real> v0 = mVertices[mIndices[3*i  ]];
        Vector2<Real> v1 = mVertices[mIndices[3*i + 1]];
        Vector2<Real> v2 = mVertices[mIndices[3*i + 2]];
        p.GetBarycentrics(v0, v1, v2, bary);
        return true;
    }

    return false;
}
//----------------------------------------------------------------------------
template <typename Real>
Delaunay2<Real>::Delaunay2 (const char* filename, int mode)
    :
    Delaunay<Real>(0, (Real)0, false, Query::QT_REAL),
    mVertices(0),
    mSVertices(0),
    mQuery(0),
    mPath(0)
{
    bool loaded = Load(filename, mode);
    assertion(loaded, "Cannot open file %s\n", filename);
    WM5_UNUSED(loaded);
}
//----------------------------------------------------------------------------
template <typename Real>
bool Delaunay2<Real>::Load (const char* filename, int mode)
{
    FileIO inFile(filename, mode);
    if (!inFile)
    {
        return false;
    }

    Delaunay<Real>::Load(inFile);

    delete0(mQuery);
    delete1(mSVertices);
    delete1(mPath);
    if (mOwner)
    {
        delete1(mVertices);
    }

    mOwner = true;
    mVertices = new1<Vector2<Real> >(mNumVertices);
    mSVertices = new1<Vector2<Real> >(mNumVertices);
    mPath = new1<int>(mNumSimplices + 1);

    inFile.Read(sizeof(int), &mNumUniqueVertices);
    inFile.Read(sizeof(int), &mPathLast);
    inFile.Read(sizeof(int), &mLastEdgeV0);
    inFile.Read(sizeof(int), &mLastEdgeV1);
    inFile.Read(sizeof(int), &mLastEdgeOpposite);
    inFile.Read(sizeof(int), &mLastEdgeOppositeIndex);
    inFile.Read(sizeof(int), mNumSimplices + 1, mPath);

    inFile.Read(sizeof(Real), 2*mNumVertices, mVertices);
    inFile.Read(sizeof(Real), 2*mNumVertices, mSVertices);
    inFile.Read(sizeof(Real), 2, &mMin);
    inFile.Read(sizeof(Real), 2, &mScale);
    inFile.Read(sizeof(Real), 2, &mLineOrigin);
    inFile.Read(sizeof(Real), 2, &mLineDirection);

    inFile.Close();

    switch (mQueryType)
    {
    case Query::QT_INT64:
    {
        mQuery = new0 Query2Int64<Real>(mNumVertices, mSVertices);
        break;
    }
    case Query::QT_INTEGER:
    {
        mQuery = new0 Query2Integer<Real>(mNumVertices, mSVertices);
        break;
    }
    case Query::QT_RATIONAL:
    {
        mQuery = new0 Query2Rational<Real>(mNumVertices, mSVertices);
        break;
    }
    case Query::QT_REAL:
    {
        mQuery = new0 Query2<Real>(mNumVertices, mSVertices);
        break;
    }
    case Query::QT_FILTERED:
    {
        mQuery = new0 Query2Filtered<Real>(mNumVertices, mSVertices,
            mEpsilon);
        break;
    }
    }

    return true;
}
//----------------------------------------------------------------------------
template <typename Real>
bool Delaunay2<Real>::Save (const char* filename, int mode) const
{
    FileIO outFile(filename, mode);
    if (!outFile)
    {
        return false;
    }

    Delaunay<Real>::Save(outFile);

    outFile.Write(sizeof(int), &mNumUniqueVertices);
    outFile.Write(sizeof(int), &mPathLast);
    outFile.Write(sizeof(int), &mLastEdgeV0);
    outFile.Write(sizeof(int), &mLastEdgeV1);
    outFile.Write(sizeof(int), &mLastEdgeOpposite);
    outFile.Write(sizeof(int), &mLastEdgeOppositeIndex);
    outFile.Write(sizeof(int), mNumSimplices + 1, mPath);

    outFile.Write(sizeof(Real), 2*mNumVertices, mVertices);
    outFile.Write(sizeof(Real), 2*mNumVertices, mVertices);
    outFile.Write(sizeof(Real), 2, &mMin);
    outFile.Write(sizeof(Real), 2, &mScale);
    outFile.Write(sizeof(Real), 2, &mLineOrigin);
    outFile.Write(sizeof(Real), 2, &mLineDirection);

    outFile.Close();
    return true;
}
//----------------------------------------------------------------------------
template <typename Real>
bool Delaunay2<Real>::GetContainingTriangle (int i, Triangle*& tri) const
{
    size_t const numTriangles = mTriMesh.GetTriangles().size();
    for (size_t t = 0; t < numTriangles; ++t)
    {
        int j;
        for (j = 0; j < 3; ++j)
        {
            const int edge[2] = { msIndex[j][0], msIndex[j][1] };
            if (mQuery->ToLine(i, tri->V[edge[0]], tri->V[edge[1]]) > 0)
            {
                // Point i sees edge <v0,v1> from outside the triangle.
                if (tri->T[j])
                {
                    // Traverse to the triangle sharing the face.
                    tri = tri->T[j];
                    break;
                }
                else
                {
                    // We reached a hull edge, so the point is outside the
                    // hull.  TODO (for WM6):  Once a hull data structure is
                    // in place, return tri->T[j] as the candidate for
                    // starting a search for visible hull edges.
                    return false;
                }
            }

        }

        if (j == 3)
        {
            // The point is inside all four edges, so the point is inside
            // a triangle.
            return true;
        }
    }

    assertion(false, "Unexpected termination of GetContainingTriangle\n");
    return false;
}
//----------------------------------------------------------------------------
template <typename Real>
void Delaunay2<Real>::GetAndRemoveInsertionPolygon (int i,
    std::set<Triangle*>& candidates, std::set<OrderedEdgeKey>& boundary)
{
    // Locate the triangles that make up the insertion polygon.
    ETManifoldMesh polygon;
    while (candidates.size() > 0)
    {
        Triangle* tri = *candidates.begin();
        candidates.erase(candidates.begin());

        for (int j = 0; j < 3; ++j)
        {
            Triangle* adj = tri->T[j];
            if (adj && candidates.find(adj) == candidates.end())
            {
                if (mQuery->ToCircumcircle(i, adj->V[0], adj->V[1],
                    adj->V[2]) <= 0)
                {
                    // Point i is in the circumcircle.
                    candidates.insert(adj);
                }
            }
        }

        polygon.InsertTriangle(tri->V[0], tri->V[1], tri->V[2]);
        mTriMesh.RemoveTriangle(tri->V[0], tri->V[1], tri->V[2]);
    }

    // Get the boundary edges of the insertion polygon.
    const ETManifoldMesh::TMap& tmap = polygon.GetTriangles();
    ETManifoldMesh::TMapCIterator element;
    for (element = tmap.begin(); element != tmap.end(); ++element)
    {
        const ETManifoldMesh::Triangle* tri = element->second;
        for (int j = 0; j < 3; ++j)
        {
            if (!tri->T[j])
            {
                const int edge[2] = { msIndex[j][0], msIndex[j][1] };
                boundary.insert(OrderedEdgeKey(tri->V[edge[0]],
                    tri->V[edge[1]]));
            }
        }
    }
}
//----------------------------------------------------------------------------
template <typename Real>
void Delaunay2<Real>::Update (int i)
{
    const ETManifoldMesh::TMap& tmap = mTriMesh.GetTriangles();
    ETManifoldMesh::Triangle* tri = tmap.begin()->second;
    if (GetContainingTriangle(i, tri))
    {
        // The point is inside the convex hull.  The insertion polygon
        // contains only triangles in the current triangulation; the
        // hull does not change.

        // Use a depth-first search for those triangles whose circumcircles
        // contain point i.
        std::set<Triangle*> candidates;
        candidates.insert(tri);

        // Get the boundary of the insertion polygon C that contains the
        // triangles whose circumcircles contain point i.  C contains the
        // point i.
        std::set<OrderedEdgeKey> boundary;
        GetAndRemoveInsertionPolygon(i, candidates, boundary);

        // The insertion polygon consists of the triangles formed by
        // point i and the faces of C.
        std::set<OrderedEdgeKey>::const_iterator key = boundary.begin();
        for (key = boundary.begin(); key != boundary.end(); ++key)
        {
            if (mQuery->ToLine(i, key->V[0], key->V[1]) < 0)
            {
                mTriMesh.InsertTriangle(i, key->V[0], key->V[1]);
            }
            // else:  Point i is on an edge of 'tri', so the
            // subdivision has degenerate triangles.  Ignore these.
        }
    }
    else
    {
        // The point is outside the convex hull.  The insertion polygon
        // is formed by point i and any triangles in the current
        // triangulation whose circumcircles contain point i.

        // Locate the convex hull of the triangles.  TODO:  In WM6, maintain
        // a hull data structure that is updated incrementally.
        std::set<OrderedEdgeKey> hull;
        const ETManifoldMesh::TMap& ttmap = mTriMesh.GetTriangles();
        ETManifoldMesh::TMapCIterator element;
        for (element = ttmap.begin(); element != ttmap.end(); ++element)
        {
            const ETManifoldMesh::Triangle* ttri = element->second;
            for (int j = 0; j < 3; ++j)
            {
                if (!ttri->T[j])
                {
                    const int edge[2] = { msIndex[j][0], msIndex[j][1] };
                    hull.insert(OrderedEdgeKey(ttri->V[edge[0]],
                        ttri->V[edge[1]]));
                }
            }
        }

        // TODO:  Until the hull change in WM6, for now just iterate over all
        // the hull edges and use the ones visible to point i to locate the
        // insertion polygon.
        const ETManifoldMesh::EMap& edgemap = mTriMesh.GetEdges();
        std::set<Triangle*> candidates;
        std::set<OrderedEdgeKey> visible;
        std::set<OrderedEdgeKey>::const_iterator key;
        for (key = hull.begin(); key != hull.end(); ++key)
        {
            if (mQuery->ToLine(i, key->V[0], key->V[1]) > 0)
            {
                ETManifoldMesh::EMapCIterator iter =
                    edgemap.find(EdgeKey(key->V[0], key->V[1]));
                assertion(iter != edgemap.end(), "Unexpected condition\n");
                assertion(iter->second->T[1] == 0, "Unexpected condition\n");
                Triangle* adj = iter->second->T[0];
                if (adj && candidates.find(adj) == candidates.end())
                {
                    if (mQuery->ToCircumcircle(i, adj->V[0], adj->V[1],
                        adj->V[2]) <= 0)
                    {
                        // Point i is in the circumcircle.
                        candidates.insert(adj);
                    }
                    else
                    {
                        // Point i is not in the circumcircle but the hull edge
                        // is visible.
                        visible.insert(*key);
                    }
                }
            }
        }

        // Get the boundary of the insertion subpolygon C that contains the
        // triangles whose circumcircles contain point i.
        std::set<OrderedEdgeKey> boundary;
        GetAndRemoveInsertionPolygon(i, candidates, boundary);

        // The insertion polygon P consists of the triangles formed by
        // point i and the back edges of C *and* the visible edges of
        // mTriMesh-C.
        for (key = boundary.begin(); key != boundary.end(); ++key)
        {
            if (mQuery->ToLine(i, key->V[0], key->V[1]) < 0)
            {
                // This is a back edge of the boundary.
                mTriMesh.InsertTriangle(i, key->V[0], key->V[1]);
            }
        }
        for (key = visible.begin(); key != visible.end(); ++key)
        {
            mTriMesh.InsertTriangle(i, key->V[1], key->V[0]);
        }
    }
}
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Explicit instantiation.
//----------------------------------------------------------------------------
template WM5_MATHEMATICS_ITEM
class Delaunay2<float>;

template WM5_MATHEMATICS_ITEM
class Delaunay2<double>;
//----------------------------------------------------------------------------
}
