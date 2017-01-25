// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <Mathematics/GteVector3.h>
#include <set>
#include <vector>

// The Polyhedron3 object represents a simple polyhedron.  The 'vertexPool'
// array can contain more points than needed to define the polyhedron, which
// allows the vertex pool to have multiple polyhedra associated with it.
// Thus, the programmer must ensure that the vertex pool persists as long as
// any Polyhedron3 objects exist that depend on the pool.  The number of
// polyhedron indices is 'numIndices' and must be 6 or larger  The 'indices'
// array refers to the points in 'vertexPool' that form the triangle faces,
// so 'numIndices' must be a multiple of 3.  The number of vertices is
// the number of unique elements in 'indices' and is determined during
// construction.  The programmer should ensure the polyhedron is simple.  The
// geometric queries are valid regardless of whether the polyhedron triangles
// are oriented clockwise or counterclockwise.
//
// NOTE:  Comparison operators are not provided.  The semantics of equal
// polyhedra is complicated and (at the moment) not useful.  The vertex pools
// can be different and indices do not match, but the vertices they reference
// can match.  Even with a shared vertex pool, the indices can be permuted,
// leading to the same polyhedron abstractly but the data structures do not
// match.

namespace gte
{

template <typename Real>
class Polyhedron3
{
public:
    // Construction and destruction.  The constructor succeeds when
    // numIndices >= 6, numVertices >= 4 (determined from the indices), and
    // 'vertices' and 'indices' are not null; we cannot test whether you have
    // a valid number of elements in the input arrays.  A copy is made of
    // 'indices', but the 'vertexPool' is not copied.  If the constructor
    // fails, the internal vertex pointer is set to null, the number of
    // vertices is set to zero, the index array has no elements, and the
    // triangle face orientation is set to clockwise.
    ~Polyhedron3();
    Polyhedron3(Vector3<Real> const* vertices, int numIndices,
        int const* indices, bool counterClockwise);

    // To validate construction, create an object as shown:
    //     Polyhedron3<Real> polyhedron(parameters);
    //     if (!polyhedron) { <constructor failed, handle accordingly>; }
    inline operator bool() const;

    // Member access.
    inline Vector3<Real> const* GetVertexPool() const;
    inline std::set<int> const& GetVertices() const;
    inline std::vector<int> const& GetIndices() const;
    inline bool CounterClockwise() const;

    // Geometric queries.
    Vector3<Real> ComputeVertexAverage() const;
    Real ComputeSurfaceArea() const;
    Real ComputeVolume() const;

private:
    Vector3<Real> const* mVertexPool;
    std::set<int> mVertices;
    std::vector<int> mIndices;
    bool mCounterClockwise;
};


template <typename Real>
Polyhedron3<Real>::~Polyhedron3()
{

}

template <typename Real>
Polyhedron3<Real>::Polyhedron3(Vector3<Real> const* vertexPool,
    int numIndices, int const* indices, bool counterClockwise)
    :
    mVertexPool(vertexPool),
    mCounterClockwise(counterClockwise)
{
    if (numIndices >= 4 && vertexPool && indices)
    {
        for (int i = 0; i < numIndices; ++i)
        {
            mVertices.insert(indices[i]);
        }

        if (numIndices == static_cast<int>(mVertices.size()))
        {
            mIndices.resize(numIndices);
            std::copy(indices, indices + numIndices, mIndices.begin());
            return;
        }

        mVertices.clear();
    }

    // Invalid input to the Polyhedron3 constructor.
    mVertexPool = nullptr;
    mCounterClockwise = false;
}

template <typename Real> inline
Polyhedron3<Real>::operator bool() const
{
    return mVertexPool != nullptr;
}

template <typename Real> inline
Vector3<Real> const* Polyhedron3<Real>::GetVertexPool() const
{
    return mVertexPool;
}

template <typename Real> inline
std::set<int> const& Polyhedron3<Real>::GetVertices() const
{
    return mVertices;
}

template <typename Real> inline
std::vector<int> const& Polyhedron3<Real>::GetIndices() const
{
    return mIndices;
}

template <typename Real> inline
bool Polyhedron3<Real>::CounterClockwise() const
{
    return mCounterClockwise;
}

template <typename Real>
Vector3<Real> Polyhedron3<Real>::ComputeVertexAverage() const
{
    Vector3<Real> average = Vector3<Real>::Zero();
    if (mVertexPool)
    {
        for (int index : mVertices)
        {
            average += mVertexPool[index];
        }
        average /= static_cast<Real>(mVertices.size());
    }
    return average;
}

template <typename Real>
Real Polyhedron3<Real>::ComputeSurfaceArea() const
{
    Real surfaceArea = (Real)0;
    if (mVertexPool)
    {
        int const numTriangles = static_cast<int>(mIndices.size()) / 3;
        int const* indices = &mIndices[0];
        for (int t = 0; t < numTriangles; ++t)
        {
            int v0 = *indices++;
            int v1 = *indices++;
            int v2 = *indices++;
            Vector3<Real> edge0 = mVertexPool[v1] - mVertexPool[v0];
            Vector3<Real> edge1 = mVertexPool[v2] - mVertexPool[v0];
            Vector3<Real> cross = Cross(edge0, edge1);
            surfaceArea += Length(cross);
        }
        surfaceArea *= (Real)0.5;
    }
    return surfaceArea;
}

template <typename Real>
Real Polyhedron3<Real>::ComputeVolume() const
{
    Real volume = (Real)0;
    if (mVertexPool)
    {
        int const numTriangles = static_cast<int>(mIndices.size()) / 3;
        int const* indices = &mIndices[0];
        for (int t = 0; t < numTriangles; ++t)
        {
            int v0 = *indices++;
            int v1 = *indices++;
            int v2 = *indices++;
            volume +=
                DotCross(mVertexPool[v0], mVertexPool[v1], mVertexPool[v2]);
        }
        volume /= (Real)6;
    }
    return std::abs(volume);
}


}
