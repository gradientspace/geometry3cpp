// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <Mathematics/GteVector2.h>
#include <set>
#include <vector>

// The Polygon2 object represents a simple polygon:  No duplicate vertices,
// closed (each vertex is shared by exactly two edges), and no
// self-intersections at interior edge points.  The 'vertexPool' array can
// contain more points than needed to define the polygon, which allows the
// vertex pool to have multiple polygons associated with it.  Thus, the
// programmer must ensure that the vertex pool persists as long as any
// Polygon2 objects exist that depend on the pool.  The number of polygon
// vertices is 'numIndices' and must be 3 or larger.  The 'indices' array
// refers to the points in 'vertexPool' that are part of the polygon and must
// have 'numIndices' unique elements.  The edges of the polygon are pairs of
// indices into 'vertexPool',
//   edge[0] = (indices[0], indices[1])
//   :
//   edge[numIndices-2] = (indices[numIndices-2], indices[numIndices-1])
//   edge[numIndices-1] = (indices[numIndices-1], indices[0])
// The programmer should ensure the polygon is simple.  The geometric
// queries are valid regardless of whether the polygon is oriented clockwise
// or counterclockwise.
//
// NOTE: Comparison operators are not provided.  The semantics of equal
// polygons is complicated and (at the moment) not useful.  The vertex pools
// can be different and indices do not match, but the vertices they reference
// can match.  Even with a shared vertex pool, the indices can be "rotated",
// leading to the same polygon abstractly but the data structures do not
// match.

namespace gte
{

template <typename Real>
class Polygon2
{
public:
    // Construction and destruction.  The constructor succeeds when
    // numVertices >= 3 and 'vertices' and 'indices' are not null; we
    // cannot test whether you have a valid number of elements in the input
    // arrays.  A copy is made of 'indices', but the 'vertexPool' is not
    // copied.  If the constructor fails, the internal vertex pointer is
    // set to null, the  index array has no elements, and the orientation
    // is set to clockwise.
    ~Polygon2();
    Polygon2(Vector2<Real> const* vertices, int numIndices,
        int const* indices, bool counterClockwise);

    // To validate construction, create an object as shown:
    //     Polygon2<Real> polygon(parameters);
    //     if (!polygon) { <constructor failed, handle accordingly>; }
    inline operator bool() const;

    // Member access.
    inline Vector2<Real> const* GetVertexPool() const;
    inline std::set<int> const& GetVertices() const;
    inline std::vector<int> const& GetIndices() const;
    inline bool CounterClockwise() const;

    // Geometric queries.
    Vector2<Real> ComputeVertexAverage() const;
    Real ComputePerimeterLength() const;
    Real ComputeArea() const;

private:
    Vector2<Real> const* mVertexPool;
    std::set<int> mVertices;
    std::vector<int> mIndices;
    bool mCounterClockwise;
};


template <typename Real>
Polygon2<Real>::~Polygon2()
{

}

template <typename Real>
Polygon2<Real>::Polygon2(Vector2<Real> const* vertexPool, int numIndices,
    int const* indices, bool counterClockwise)
    :
    mVertexPool(vertexPool),
    mCounterClockwise(counterClockwise)
{
    if (numIndices >= 3 && vertexPool && indices)
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

    // Invalid input to the Polygon2 constructor.
    mVertexPool = nullptr;
    mCounterClockwise = false;
}

template <typename Real> inline
Polygon2<Real>::operator bool() const
{
    return mVertexPool != nullptr;
}

template <typename Real> inline
Vector2<Real> const* Polygon2<Real>::GetVertexPool() const
{
    return mVertexPool;
}

template <typename Real> inline
std::set<int> const& Polygon2<Real>::GetVertices() const
{
    return mVertices;
}

template <typename Real> inline
std::vector<int> const& Polygon2<Real>::GetIndices() const
{
    return mIndices;
}

template <typename Real> inline
bool Polygon2<Real>::CounterClockwise() const
{
    return mCounterClockwise;
}

template <typename Real>
Vector2<Real> Polygon2<Real>::ComputeVertexAverage() const
{
    Vector2<Real> average = Vector2<Real>::Zero();
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
Real Polygon2<Real>::ComputePerimeterLength() const
{
    Real length = (Real)0;
    if (mVertexPool)
    {
        Vector2<Real> v0 = mVertexPool[mIndices.back()];
        for (int index : mIndices)
        {
            Vector2<Real> v1 = mVertexPool[index];
            length += Length(v1 - v0);
            v0 = v1;
        }
    }
    return length;
}

template <typename Real>
Real Polygon2<Real>::ComputeArea() const
{
    Real area = (Real)0;
    if (mVertexPool)
    {
        int const numIndices = static_cast<int>(mIndices.size());
        Vector2<Real> v0 = mVertexPool[mIndices[numIndices - 2]];
        Vector2<Real> v1 = mVertexPool[mIndices[numIndices - 1]];
        for (int index : mIndices)
        {
            Vector2<Real> v2 = mVertexPool[index];
            area += v1[0] * (v2[1] - v0[1]);
            v0 = v1;
            v1 = v2;
        }
        area *= (Real)0.5;
    }
    return std::abs(area);
}


}
