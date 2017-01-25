// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.9 (2015/11/21)

#include "Wm5MathematicsPCH.h"
#include "Wm5Delaunay3.h"
#include "Wm5Delaunay2.h"
#include "Wm5Query3Filtered.h"
#include "Wm5Query3Int64.h"
#include "Wm5Query3Integer.h"
#include "Wm5Query3Rational.h"

namespace Wm5
{
    
//----------------------------------------------------------------------------
template <typename Real>
Delaunay3<Real>::Delaunay3 (int numVertices, Vector3<Real>* vertices,
    Real epsilon, bool owner, Query::Type queryType)
    :
    Delaunay<Real>(numVertices, epsilon, owner, queryType),
    mVertices(vertices),
    mNumUniqueVertices(0),
    mSVertices(0),
    mQuery(0),
    mLineOrigin(Vector3<Real>::ZERO),
    mLineDirection(Vector3<Real>::ZERO),
    mPlaneOrigin(Vector3<Real>::ZERO),
    mPathLast(-1),
    mPath(0),
    mLastFaceV0(-1),
    mLastFaceV1(-1),
    mLastFaceV2(-1),
    mLastFaceOpposite(-1),
    mLastFaceOppositeIndex(-1)
{
    mPlaneDirection[0] = Vector3<Real>::ZERO;
    mPlaneDirection[1] = Vector3<Real>::ZERO;

    typename Vector3<Real>::Information info;
    Vector3<Real>::GetInformation(mNumVertices, mVertices, mEpsilon, info);
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

    if (info.mDimension == 2)
    {
        // The set is (nearly) coplanar.  The caller is responsible for
        // creating a Delaunay2 object.
        mDimension = 2;
        mPlaneOrigin = info.mOrigin;
        mPlaneDirection[0] = info.mDirection[0];
        mPlaneDirection[1] = info.mDirection[1];
        return;
    }

    mDimension = 3;

    // Allocate storage for the input vertices.
    mSVertices = new1<Vector3<Real> >(mNumVertices);
    int i;

    if (queryType != Query::QT_RATIONAL && queryType != Query::QT_FILTERED)
    {
        // Transform the vertices to the cube [0,1]^3.
        mMin = Vector3<Real>(info.mMin[0], info.mMin[1], info.mMin[2]);
        mScale = ((Real)1)/info.mMaxRange;
        for (i = 0; i < mNumVertices; ++i)
        {
            mSVertices[i] = (mVertices[i] - mMin)*mScale;
        }

        Real expand;
        if (queryType == Query::QT_INT64)
        {
            // Scale the vertices to the cube [0,2^{10}]^3 to allow use of
            // 64-bit integers for tetrahedralization.
            expand = (Real)(1 << 10);
            mQuery = new0 Query3Int64<Real>(mNumVertices, mSVertices);
        }
        else if (queryType == Query::QT_INTEGER)
        {
            // Scale the vertices to the cube [0,2^{20}]^3 to get more
            // precision for TInteger than for 64-bit integers for
            // tetrahedralization.
            expand = (Real)(1 << 20);
            mQuery = new0 Query3Integer<Real>(mNumVertices, mSVertices);
        }
        else // queryType == Query::QT_REAL
        {
            // No scaling for floating point.
            expand = (Real)1;
            mQuery = new0 Query3<Real>(mNumVertices, mSVertices);
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
        mMin = Vector3<Real>::ZERO;
        mScale = (Real)1;
        memcpy(mSVertices, mVertices, mNumVertices*sizeof(Vector3<Real>));

        if (queryType == Query::QT_RATIONAL)
        {
            mQuery = new0 Query3Rational<Real>(mNumVertices, mSVertices);
        }
        else // queryType == Query::QT_FILTERED
        {
            mQuery = new0 Query3Filtered<Real>(mNumVertices, mSVertices,
                mEpsilon);
        }
    }

    // Insert the (nondegenerate) tetrahedron constructed by the call to
    // GetInformation. This is necessary for the circumsphere-visibility
    // algorithm to work correctly.
    if (!info.mExtremeCCW)
    {
        std::swap(info.mExtreme[2], info.mExtreme[3]);
    }
    mTetraMesh.Insert(info.mExtreme[0], info.mExtreme[1], info.mExtreme[2],
        info.mExtreme[3]);

    // Incrementally update the tetrahedralization.  The set of processed
    // points is maintained to eliminate duplicates, either in the original
    // input points or in the points obtained by snap rounding.
    std::set<Vector3<Real> > processed;
    for (i = 0; i < 4; ++i)
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

    // Assign integer values to the tetrahedra for use by the caller.
    std::map<Tetrahedron*,int> permute;
    i = -1;
    permute[(Tetrahedron*)0] = i++;
    const TSManifoldMesh::SMap& smap = mTetraMesh.GetTetrahedra();
    TSManifoldMesh::SMapCIterator element;
    for (element = smap.begin(); element != smap.end(); ++element)
    {
        permute[element->second] = i++;
    }

    // Put Delaunay tetrahedra into an array (vertices and adjacency info).
    mNumSimplices = (int)mTetraMesh.GetTetrahedra().size();
    if (mNumSimplices > 0)
    {
        mIndices = new1<int>(4*mNumSimplices);
        mAdjacencies = new1<int>(4*mNumSimplices);
        i = 0;
        for (element = smap.begin(); element != smap.end(); ++element)
        {
            const TSManifoldMesh::Tetrahedron* tetra = element->second;
            for (int j = 0; j < 4; ++j, ++i)
            {
                mIndices[i] = tetra->V[j];
                mAdjacencies[i] = permute[tetra->S[j]];
            }
        }
        assertion(i == 4*mNumSimplices, "Unexpected mismatch\n");

        mPathLast = -1;
        mPath = new1<int>(mNumSimplices + 1);
        memset(mPath, 0, (mNumSimplices + 1)*sizeof(int));
    }
}
//----------------------------------------------------------------------------
template <typename Real>
Delaunay3<Real>::~Delaunay3 ()
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
const Vector3<Real>* Delaunay3<Real>::GetVertices () const
{
    return mVertices;
}
//----------------------------------------------------------------------------
template <typename Real>
int Delaunay3<Real>::GetNumUniqueVertices () const
{
    return mNumUniqueVertices;
}
//----------------------------------------------------------------------------
template <typename Real>
const Vector3<Real>& Delaunay3<Real>::GetLineOrigin () const
{
    return mLineOrigin;
}
//----------------------------------------------------------------------------
template <typename Real>
const Vector3<Real>& Delaunay3<Real>::GetLineDirection () const
{
    return mLineDirection;
}
//----------------------------------------------------------------------------
template <typename Real>
Delaunay1<Real>* Delaunay3<Real>::GetDelaunay1 () const
{
    assertion(mDimension == 1, "The dimension must be 1\n");
    if (mDimension != 1)
    {
        return 0;
    }

    Real* projection = new1<Real>(mNumVertices);
    for (int i = 0; i < mNumVertices; ++i)
    {
        Vector3<Real> diff = mVertices[i] - mLineOrigin;
        projection[i] = mLineDirection.Dot(diff);
    }

    return new0 Delaunay1<Real>(mNumVertices, projection, mEpsilon, true,
        mQueryType);
}
//----------------------------------------------------------------------------
template <typename Real>
const Vector3<Real>& Delaunay3<Real>::GetPlaneOrigin () const
{
    return mPlaneOrigin;
}
//----------------------------------------------------------------------------
template <typename Real>
const Vector3<Real>& Delaunay3<Real>::GetPlaneDirection (int i) const
{
    return mPlaneDirection[i];
}
//----------------------------------------------------------------------------
template <typename Real>
Delaunay2<Real>* Delaunay3<Real>::GetDelaunay2 () const
{
    assertion(mDimension == 2, "The dimension must be 2\n");
    if (mDimension != 2)
    {
        return 0;
    }

    Vector2<Real>* projection = new1<Vector2<Real> >(mNumVertices);
    for (int i = 0; i < mNumVertices; ++i)
    {
        Vector3<Real> diff = mVertices[i] - mPlaneOrigin;
        projection[i][0] = mPlaneDirection[0].Dot(diff);
        projection[i][1] = mPlaneDirection[1].Dot(diff);
    }

    return new0 Delaunay2<Real>(mNumVertices, projection, mEpsilon, true,
        mQueryType);
}
//----------------------------------------------------------------------------
template <typename Real>
bool Delaunay3<Real>::GetHull (int& numTriangles, int*& indices) const
{
    assertion(mDimension == 3, "The dimension must be 3\n");
    if (mDimension != 3)
    {
        return false;
    }

    numTriangles = 0;
    indices = 0;

    // Count the number of triangles that are not shared by two tetrahedra.
    int i, numAdjacent = 4*mNumSimplices;
    for (i = 0; i < numAdjacent; ++i)
    {
        if (mAdjacencies[i] == -1)
        {
            numTriangles++;
        }
    }
    assertion(numTriangles > 0, "There must be at least one tetrahedron\n");
    if (numTriangles == 0)
    {
        return false;
    }

    // Enumerate the triangles.  The prototypical case is the single
    // tetrahedron V[0] = (0,0,0), V[1] = (1,0,0), V[2] = (0,1,0), and
    // V[3] = (0,0,1) with no adjacent tetrahedra.  The mIndices[] array
    // is <0,1,2,3>.
    //   i = 0, face = 0:
    //     skip index 0, <x,1,2,3>, no swap, triangle = <1,2,3>
    //   i = 1, face = 1:
    //     skip index 1, <0,x,2,3>, swap,    triangle = <0,3,2>
    //   i = 2, face = 2:
    //     skip index 2, <0,1,x,3>, no swap, triangle = <0,1,3>
    //   i = 3, face = 3:
    //     skip index 3, <0,1,2,x>, swap,    triangle = <0,2,1>
    // To guarantee counterclockwise order of triangles when viewed outside
    // the tetrahedron, the swap of the last two indices occurs when
    // iFace is an odd number:  (iFace % 2) != 0
    indices = new1<int>(3*numTriangles);
    int* currentIndex = indices;
    for (i = 0; i < numAdjacent; ++i)
    {
        if (mAdjacencies[i] == -1)
        {
            int tetra = i/4, face = i%4;
            for (int j = 0; j < 4; ++j)
            {
                if (j != face)
                {
                    *currentIndex++ = mIndices[4*tetra + j];
                }
            }
            if ((face % 2) != 0)
            {
                int save = *(currentIndex-1);
                *(currentIndex-1) = *(currentIndex-2);
                *(currentIndex-2) = save;
            }
        }
    }

    return true;
}
//----------------------------------------------------------------------------
template <typename Real>
int Delaunay3<Real>::GetContainingTetrahedron (const Vector3<Real>& p) const
{
    assertion(mDimension == 3, "The dimension must be 3\n");
    if (mDimension != 3)
    {
        return -1;
    }

    // Convert to scaled coordinates.
    Vector3<Real> scP = (p - mMin)*mScale;

    // Start at first tetrahedron in mesh.
    int index = (mPathLast >= 0 ? mPath[mPathLast] : 0);
    mPathLast = -1;
    mLastFaceV0 = -1;
    mLastFaceV1 = -1;
    mLastFaceV2 = -1;
    mLastFaceOpposite = -1;
    mLastFaceOppositeIndex = -1;

    // Use tetrahedron faces as binary separating planes.
    for (int i = 0; i < mNumSimplices; ++i)
    {
        mPath[++mPathLast] = index;

        int* vertices = &mIndices[4*index];

        // <V1,V2,V3> counterclockwise when viewed outside tetrahedron.
        if (mQuery->ToPlane(scP, vertices[1], vertices[2], vertices[3]) > 0)
        {
            index = mAdjacencies[4*index];
            if (index == -1)
            {
                mLastFaceV0 = vertices[1];
                mLastFaceV1 = vertices[2];
                mLastFaceV2 = vertices[3];
                mLastFaceOpposite = vertices[0];
                mLastFaceOppositeIndex = 0;
                return -1;
            }
            continue;
        }

        // <V0,V3,V2> counterclockwise when viewed outside tetrahedron.
        if (mQuery->ToPlane(scP, vertices[0], vertices[2], vertices[3]) < 0)
        {
            index = mAdjacencies[4*index + 1];
            if (index == -1)
            {
                mLastFaceV0 = vertices[0];
                mLastFaceV1 = vertices[2];
                mLastFaceV2 = vertices[3];
                mLastFaceOpposite = vertices[1];
                mLastFaceOppositeIndex = 1;
                return -1;
            }
            continue;
        }

        // <V0,V1,V3> counterclockwise when viewed outside tetrahedron.
        if (mQuery->ToPlane(scP, vertices[0], vertices[1], vertices[3]) > 0)
        {
            index = mAdjacencies[4*index + 2];
            if (index == -1)
            {
                mLastFaceV0 = vertices[0];
                mLastFaceV1 = vertices[1];
                mLastFaceV2 = vertices[3];
                mLastFaceOpposite = vertices[2];
                mLastFaceOppositeIndex = 2;
                return -1;
            }
            continue;
        }

        // <V0,V2,V1> counterclockwise when viewed outside tetrahedron.
        if (mQuery->ToPlane(scP, vertices[0], vertices[1], vertices[2]) < 0)
        {
            index = mAdjacencies[4*index + 3];
            if (index == -1)
            {
                mLastFaceV0 = vertices[0];
                mLastFaceV1 = vertices[1];
                mLastFaceV2 = vertices[2];
                mLastFaceOpposite = vertices[3];
                mLastFaceOppositeIndex = 3;
                return -1;
            }
            continue;
        }

        mLastFaceV0 = -1;
        mLastFaceV1 = -1;
        mLastFaceV2 = -1;
        mLastFaceOppositeIndex = -1;
        return index;
    }

    return -1;
}
//----------------------------------------------------------------------------
template <typename Real>
int Delaunay3<Real>::GetPathLast () const
{
    return mPathLast;
}
//----------------------------------------------------------------------------
template <typename Real>
const int* Delaunay3<Real>::GetPath () const
{
    return mPath;
}
//----------------------------------------------------------------------------
template <typename Real>
int Delaunay3<Real>::GetLastFace (int& v0, int& v1, int& v2, int& v3) const
{
    v0 = mLastFaceV0;
    v1 = mLastFaceV1;
    v2 = mLastFaceV2;
    v3 = mLastFaceOpposite;
    return mLastFaceOppositeIndex;
}
//----------------------------------------------------------------------------
template <typename Real>
bool Delaunay3<Real>::GetVertexSet (int i, Vector3<Real> vertices[4]) const
{
    assertion(mDimension == 3, "The dimension must be 3\n");
    if (mDimension != 3)
    {
        return false;
    }

    if (0 <= i && i < mNumSimplices)
    {
        vertices[0] = mVertices[mIndices[4*i  ]];
        vertices[1] = mVertices[mIndices[4*i + 1]];
        vertices[2] = mVertices[mIndices[4*i + 2]];
        vertices[3] = mVertices[mIndices[4*i + 3]];
        return true;
    }

    return false;
}
//----------------------------------------------------------------------------
template <typename Real>
bool Delaunay3<Real>::GetIndexSet (int i, int indices[4]) const
{
    assertion(mDimension == 3, "The dimension must be 3\n");
    if (mDimension != 3)
    {
        return false;
    }

    if (0 <= i && i < mNumSimplices)
    {
        indices[0] = mIndices[4*i  ];
        indices[1] = mIndices[4*i + 1];
        indices[2] = mIndices[4*i + 2];
        indices[3] = mIndices[4*i + 3];
        return true;
    }

    return false;
}
//----------------------------------------------------------------------------
template <typename Real>
bool Delaunay3<Real>::GetAdjacentSet (int i, int adjacencies[4]) const
{
    assertion(mDimension == 3, "The dimension must be 3\n");
    if (mDimension != 3)
    {
        return false;
    }

    if (0 <= i && i < mNumSimplices)
    {
        adjacencies[0] = mAdjacencies[4*i  ];
        adjacencies[1] = mAdjacencies[4*i + 1];
        adjacencies[2] = mAdjacencies[4*i + 2];
        adjacencies[3] = mAdjacencies[4*i + 3];
        return true;
    }

    return false;
}
//----------------------------------------------------------------------------
template <typename Real>
bool Delaunay3<Real>::GetBarycentricSet (int i, const Vector3<Real>& p,
    Real bary[4]) const
{
    assertion(mDimension == 3, "The dimension must be 3\n");
    if (mDimension != 3)
    {
        return false;
    }

    if (0 <= i && i < mNumSimplices)
    {
        Vector3<Real> v0 = mVertices[mIndices[4*i  ]];
        Vector3<Real> v1 = mVertices[mIndices[4*i + 1]];
        Vector3<Real> v2 = mVertices[mIndices[4*i + 2]];
        Vector3<Real> v3 = mVertices[mIndices[4*i + 3]];
        p.GetBarycentrics(v0, v1, v2, v3, bary);
        return true;
    }

    return false;
}
//----------------------------------------------------------------------------
template <typename Real>
Delaunay3<Real>::Delaunay3 (const char* filename, int mode)
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
bool Delaunay3<Real>::Load (const char* filename, int mode)
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
    mVertices = new1<Vector3<Real> >(mNumVertices);
    mSVertices = new1<Vector3<Real> >(mNumVertices);
    mPath = new1<int>(mNumSimplices + 1);

    inFile.Read(sizeof(int), &mNumUniqueVertices);
    inFile.Read(sizeof(int), &mPathLast);
    inFile.Read(sizeof(int), &mLastFaceV0);
    inFile.Read(sizeof(int), &mLastFaceV1);
    inFile.Read(sizeof(int), &mLastFaceV2);
    inFile.Read(sizeof(int), &mLastFaceOpposite);
    inFile.Read(sizeof(int), &mLastFaceOppositeIndex);
    inFile.Read(sizeof(int), mNumSimplices + 1, mPath);

    inFile.Read(sizeof(Real), 3*mNumVertices, mVertices);
    inFile.Read(sizeof(Real), 3*mNumVertices, mSVertices);
    inFile.Read(sizeof(Real), 3, &mMin);
    inFile.Read(sizeof(Real), 3, &mScale);
    inFile.Read(sizeof(Real), 3, &mLineOrigin);
    inFile.Read(sizeof(Real), 3, &mLineDirection);
    inFile.Read(sizeof(Real), 3, &mPlaneOrigin);
    inFile.Read(sizeof(Real), 6, mPlaneDirection);

    inFile.Close();

    switch (mQueryType)
    {
    case Query::QT_INT64:
    {
        mQuery = new0 Query3Int64<Real>(mNumVertices, mSVertices);
        break;
    }
    case Query::QT_INTEGER:
    {
        mQuery = new0 Query3Integer<Real>(mNumVertices, mSVertices);
        break;
    }
    case Query::QT_RATIONAL:
    {
        mQuery = new0 Query3Rational<Real>(mNumVertices, mSVertices);
        break;
    }
    case Query::QT_REAL:
    {
        mQuery = new0 Query3<Real>(mNumVertices, mSVertices);
        break;
    }
    case Query::QT_FILTERED:
    {
        mQuery = new0 Query3Filtered<Real>(mNumVertices, mSVertices,
            mEpsilon);
        break;
    }
    }

    return true;
}
//----------------------------------------------------------------------------
template <typename Real>
bool Delaunay3<Real>::Save (const char* filename, int mode) const
{
    FileIO outFile(filename, mode);
    if (!outFile)
    {
        return false;
    }

    Delaunay<Real>::Save(outFile);

    outFile.Write(sizeof(int), &mNumUniqueVertices);
    outFile.Write(sizeof(int), &mPathLast);
    outFile.Write(sizeof(int), &mLastFaceV0);
    outFile.Write(sizeof(int), &mLastFaceV1);
    outFile.Write(sizeof(int), &mLastFaceV2);
    outFile.Write(sizeof(int), &mLastFaceOpposite);
    outFile.Write(sizeof(int), &mLastFaceOpposite);
    outFile.Write(sizeof(int), mNumSimplices + 1 ,mPath);

    outFile.Write(sizeof(Real), 3*mNumVertices, mVertices);
    outFile.Write(sizeof(Real), 3*mNumVertices, mSVertices);
    outFile.Write(sizeof(Real), 3, &mMin);
    outFile.Write(sizeof(Real), 3, &mScale);
    outFile.Write(sizeof(Real), 3, &mLineOrigin);
    outFile.Write(sizeof(Real), 3, &mLineDirection);
    outFile.Write(sizeof(Real), 3, &mPlaneOrigin);
    outFile.Write(sizeof(Real), 6, mPlaneDirection);

    outFile.Close();
    return true;
}
//----------------------------------------------------------------------------
template <typename Real>
bool Delaunay3<Real>::GetContainingTetrahedron (int i, Tetrahedron*& tetra)
    const
{
    size_t const numTetra = mTetraMesh.GetTetrahedra().size();
    for (size_t t = 0; t < numTetra; ++t)
    {
        int j;
        for (j = 0; j < 4; ++j)
        {
            const int face[3] = {
                TetrahedronKey::oppositeFace[j][0],
                TetrahedronKey::oppositeFace[j][1],
                TetrahedronKey::oppositeFace[j][2] };
            if (mQuery->ToPlane(i, tetra->V[face[0]], tetra->V[face[1]],
                tetra->V[face[2]]) > 0)
            {
                // Point i sees face <v0,v1,v2> from outside the tetrahedron.
                if (tetra->S[j])
                {
                    // Traverse to the tetrahedron sharing the face.
                    tetra = tetra->S[j];
                    break;
                }
                else
                {
                    // We reached a hull face, so the point is outside the
                    // hull.  TODO (for WM6):  Once a hull data structure is
                    // in place, return tetra->S[j] as the candidate for
                    // starting a search for visible hull faces.
                    return false;
                }
            }

        }

        if (j == 4)
        {
            // The point is inside all four faces, so the point is inside
            // a tetrahedron.
            return true;
        }
    }

    assertion(false, "Unexpected termination of GetContainingTetrahedron\n");
    return false;
}
//----------------------------------------------------------------------------
template <typename Real>
void Delaunay3<Real>::GetAndRemoveInsertionPolyhedron (int i,
    std::set<Tetrahedron*>& candidates, std::set<TriangleKey>& boundary)
{
    // Locate the tetrahedra that make up the insertion polyhedron.
    TSManifoldMesh polyhedron;
    while (candidates.size() > 0)
    {
        Tetrahedron* tetra = *candidates.begin();
        candidates.erase(candidates.begin());

        for (int j = 0; j < 4; ++j)
        {
            Tetrahedron* adj = tetra->S[j];
            if (adj && candidates.find(adj) == candidates.end())
            {
                if (mQuery->ToCircumsphere(i, adj->V[0], adj->V[1], adj->V[2],
                    adj->V[3]) <= 0)
                {
                    // Point i is in the circumsphere.
                    candidates.insert(adj);
                }
            }
        }

        polyhedron.Insert(tetra->V[0], tetra->V[1], tetra->V[2], tetra->V[3]);
        mTetraMesh.Remove(tetra->V[0], tetra->V[1], tetra->V[2], tetra->V[3]);
    }

    // Get the boundary triangles of the insertion polyhedron.
	const TSManifoldMesh::SMap& smap = polyhedron.GetTetrahedra();
    TSManifoldMesh::SMapCIterator element;
	for (element = smap.begin(); element != smap.end(); ++element)
    {
        const TSManifoldMesh::Tetrahedron* tetra = element->second;
        for (int j = 0; j < 4; ++j)
        {
            if (!tetra->S[j])
            {
                const int face[3] = {
                    TetrahedronKey::oppositeFace[j][0],
                    TetrahedronKey::oppositeFace[j][1],
                    TetrahedronKey::oppositeFace[j][2] };
                boundary.insert(TriangleKey(tetra->V[face[0]],
                    tetra->V[face[1]], tetra->V[face[2]]));
            }
        }
    }
}
//----------------------------------------------------------------------------
template <typename Real>
void Delaunay3<Real>::Update (int i)
{
    const TSManifoldMesh::SMap& smap = mTetraMesh.GetTetrahedra();
    TSManifoldMesh::Tetrahedron* tetra = smap.begin()->second;
    if (GetContainingTetrahedron(i, tetra))
    {
        // The point is inside the convex hull.  The insertion polyhedron
        // contains only tetrahedra in the current tetrahedralization; the
        // hull does not change.

        // Use a depth-first search for those tetrahedra whose circumspheres
        // contain point i.
        std::set<Tetrahedron*> candidates;
        candidates.insert(tetra);

        // Get the boundary of the insertion polyhedron C that contains the
        // tetrahedra whose circumspheres contain point i.  C contains the
        // point i.
        std::set<TriangleKey> boundary;
        GetAndRemoveInsertionPolyhedron(i, candidates, boundary);

        // The insertion polyhedron consists of the tetrahedra formed by
        // point i and the faces of C.
        std::set<TriangleKey>::const_iterator key;
        for (key = boundary.begin(); key != boundary.end(); ++key)
        {
            if (mQuery->ToPlane(i, key->V[0], key->V[1], key->V[2]) < 0)
            {
                mTetraMesh.Insert(i, key->V[0], key->V[1], key->V[2]);
            }
            // else:  Point i is on an edge or face of 'tetra', so the
            // subdivision has degenerate tetrahedra.  Ignore these.
        }
    }
    else
    {
        // The point is outside the convex hull.  The insertion polyhedron
        // is formed by point i and any tetrahedra in the current
        // tetrahedralization whose circumspheres contain point i.

        // Locate the convex hull of the tetrahedra.  TODO:  In WM6, maintain
        // a hull data structure that is updated incrementally.
        std::set<TriangleKey> hull;
        const TSManifoldMesh::SMap& ssmap = mTetraMesh.GetTetrahedra();
        TSManifoldMesh::SMapCIterator element;
        for (element = ssmap.begin(); element != ssmap.end(); ++element)
        {
            const TSManifoldMesh::Tetrahedron* ttetra = element->second;
            for (int j = 0; j < 4; ++j)
            {
                if (!ttetra->S[j])
                {
                    const int face[3] = {
                        TetrahedronKey::oppositeFace[j][0],
                        TetrahedronKey::oppositeFace[j][1],
                        TetrahedronKey::oppositeFace[j][2] };
                    hull.insert(TriangleKey(ttetra->V[face[0]],
                        ttetra->V[face[1]], ttetra->V[face[2]]));
                }
            }
        }

        // TODO:  Until the hull change in WM6, for now just iterate over all
        // the hull faces and use the ones visible to point i to locate the
        // insertion polyhedron.
        const TSManifoldMesh::TMap& trimap = mTetraMesh.GetTriangles();
        std::set<Tetrahedron*> candidates;
        std::set<TriangleKey> visible;
        std::set<TriangleKey>::const_iterator key;
        for (key = hull.begin(); key != hull.end(); ++key)
        {
            if (mQuery->ToPlane(i, key->V[0], key->V[1], key->V[2]) > 0)
            {
                TSManifoldMesh::TMapCIterator iter = trimap.find(
                    UnorderedTriangleKey(key->V[0], key->V[1], key->V[2]));
                assertion(iter != trimap.end(), "Unexpected condition\n");
                assertion(iter->second->T[1] == 0, "Unexpected condition\n");
                Tetrahedron* adj = iter->second->T[0];
                if (adj && candidates.find(adj) == candidates.end())
                {
                    if (mQuery->ToCircumsphere(i, adj->V[0], adj->V[1],
                        adj->V[2], adj->V[3]) <= 0)
                    {
                        // Point i is in the circumsphere.
                        candidates.insert(adj);
                    }
                    else
                    {
                        // Point i is not in the circumsphere but the hull face
                        // is visible.
                        visible.insert(*key);
                    }
                }
            }
        }

        // Get the boundary of the insertion subpolyhedron C that contains the
        // tetrahedra whose circumspheres contain point i.
        std::set<TriangleKey> boundary;
        GetAndRemoveInsertionPolyhedron(i, candidates, boundary);

        // The insertion polyhedron P consists of the tetrahedra formed by
        // point i and the back faces of C *and* the visible faces of
        // mTetraMesh-C.
        for (key = boundary.begin(); key != boundary.end(); ++key)
        {
            if (mQuery->ToPlane(i, key->V[0], key->V[1], key->V[2]) < 0)
            {
                // This is a back face of the boundary.
                mTetraMesh.Insert(i, key->V[0], key->V[1], key->V[2]);
            }
        }
        for (key = visible.begin(); key != visible.end(); ++key)
        {
            mTetraMesh.Insert(i, key->V[0], key->V[2], key->V[1]);
        }
    }
}
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Explicit instantiation.
//----------------------------------------------------------------------------
template WM5_MATHEMATICS_ITEM
class Delaunay3<float>;

template WM5_MATHEMATICS_ITEM
class Delaunay3<double>;
//----------------------------------------------------------------------------
}
