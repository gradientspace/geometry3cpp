// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <Mathematics/GteConstrainedDelaunay2.h>
#include <Mathematics/GteContPointInPolygon2.h>
#include <algorithm>
#include <queue>

// The triangulation is based on constrained Delaunay triangulations (CDT),
// which does not use divisions, so ComputeType may be chosen using BSNumber.
// The input constraints are relaxed compared to TriangulateEC; specifically,
// the inner polygons are allowed to share vertices with the outer polygons.
// The CDT produces a triangulation of the convex hull of the input, which
// includes triangles outside the top-level outer polygon and inside the
// inner polygons.  Only the triangles relevant to the input are returned
// via the 'std::vector<int>& triangles', but the other triangles are
// retained so that you can perform linear walks in search of points inside
// the original polygon (nested polygon, tree of nested polygons).  This is
// useful, for example, when subsampling the polygon triangles for
// interpolation of function data specified at the vertices.  A linear walk
// does not work for a mesh consisting only of the polygon triangles, but
// with the additional triangles, the walk can navigate through holes in
// the polygon to find the containing triangle of the specified point.

namespace gte
{

template <typename InputType, typename ComputeType>
class TriangulateCDT
{
public:
    // The class is a functor to support triangulating multiple polygons that
    // share vertices in a collection of points.  The interpretation of
    // 'numPoints' and 'points' is described before each operator() function.
    // Preconditions are numPoints >= 3 and points is a nonnull pointer to an
    // array of at least numPoints elements.  If the preconditions are
    // satisfied, then operator() functions will return 'true'; otherwise,
    // they return 'false'.
    TriangulateCDT(int numPoints, Vector2<InputType> const* points);

    // The triangles of the polygon triangulation.
    inline std::vector<std::array<int, 3>> const& GetTriangles() const;

    // The triangles inside the convex hull of the points but outside the
    // triangulation.
    inline std::vector<std::array<int, 3>> const& GetOutsideTriangles() const;

    // The triangles of the convex hull of the inputs to the constructor.
    inline std::vector<std::array<int, 3>> const& GetAllTriangles() const;

    // The classification of whether a triangle is part of the triangulation
    // or outside the triangulation.  These may be used in conjunction with
    // the array returned by GetAllTriangles().
    inline std::vector<bool> const& GetIsInside() const;
    inline bool IsInside(int triIndex) const;
    inline bool IsOutside(int triIndex) const;

    // A simple wrapper for indexed polygons.  Although the polygon inputs
    // could be std::vector<int> objects, the wrapper supports users whose
    // input is not of std::vector<int> form, thus avoiding copying of
    // their input arrays to std::vector<int>.  The outer polygons have
    // counterclockwise ordered vertices.  The inner polygons have clockwise
    // ordered vertices.
    struct Polygon
    {
        int numIndices;
        int* indices;
    };

    // The input 'points' represents an array of vertices for a simple
    // polygon. The vertices are points[0] through points[numPoints-1] and
    // are listed in counterclockwise order.
    bool operator()();

    // The input 'points' represents an array of vertices that contains the
    // vertices of a simple polygon.
    bool operator()(Polygon const& polygon);

    // The input 'points' is a shared array of vertices that contains the
    // vertices for two simple polygons, an outer polygon and an inner
    // polygon.  The inner polygon must be strictly inside the outer polygon.
    bool operator()(Polygon const& outer, Polygon const& inner);

    // The input 'points' is a shared array of vertices that contains the
    // vertices for multiple simple polygons, an outer polygon and one or more
    // inner polygons.  The inner polygons must be nonoverlapping and strictly
    // inside the outer polygon.
    bool operator()(Polygon const& outer, std::vector<Polygon> const& inners);

    // A tree of nested polygons.  The root node corresponds to an outer
    // polygon.  The children of the root correspond to inner polygons, which
    // are nonoverlapping polygons strictly contained in the outer polygon.
    // Each inner polygon may itself contain an outer polygon, thus leading
    // to a hierarchy of polygons.  The outer polygons have vertices listed
    // in counterclockwise order.  The inner polygons have vertices listed in
    // clockwise order.
    class Tree
    {
    public:
        Polygon polygon;
        std::vector<Tree*> child;
    };

    // The input 'positions' is a shared array of vertices that contains the
    // vertices for multiple simple polygons in a tree of polygons.
    bool operator()(Tree const& tree);

    // Helper function to delete a tree that was dynamically allocated.
    // The function deletes the child pointers.  You may choose to have it
    // also delete the polygon.indices pointers, assuming you dynamically
    // allocated them.
    static void Delete(Tree*& tree, bool deletePolygons);

private:
    // Triangulate the points referenced by an operator(...) query.  The
    // mAllTriangles and mIsInside are populated by this function, but the
    // indices of mAllTriangles are relative to the packed 'points'.  After
    // the call, the indices need to be mapped back to the original set
    // provided by the input arrays to operator(...).  The mTriangles and
    // mOutsideTriangles are generated after the call by examining
    // mAllTriangles and mIsInside.
    bool TriangulatePacked(int numPoints, Vector2<InputType> const* points,
        Tree const& tree, std::map<Tree const*, int> const& offsets);

    int GetNumPointsAndOffsets(Tree const& tree,
        std::map<Tree const*, int>& offsets) const;

    void PackPoints(Tree const& tree,
        std::vector<Vector2<InputType>>& points);

    bool InsertEdges(Tree const& tree);

    void LookupIndex(Tree const& tree, int& v,
        std::map<Tree const*, int> const& offsets) const;

    bool IsInside(Tree const& tree, Vector2<ComputeType> const* points,
        Vector2<ComputeType> const& test,
        std::map<Tree const*, int> const& offsets) const;

    // The input polygon.
    int mNumPoints;
    Vector2<InputType> const* mPoints;

    // The output triangulation and those triangle inside the hull of the
    // points but outside the triangulation.
    std::vector<std::array<int, 3>> mTriangles;
    std::vector<std::array<int, 3>> mOutsideTriangles;
    std::vector<std::array<int, 3>> mAllTriangles;
    std::vector<bool> mIsInside;

    ConstrainedDelaunay2<InputType, ComputeType> mConstrainedDelaunay;
};


template <typename InputType, typename ComputeType>
TriangulateCDT<InputType, ComputeType>::TriangulateCDT(int numPoints,
    Vector2<InputType> const* points)
    :
    mNumPoints(numPoints),
    mPoints(points)
{
    if (mNumPoints < 3 || !mPoints)
    {
        LogError("Invalid input.");
        mNumPoints = 0;
        mPoints = nullptr;
        // The operator() functions will triangulate only when mPoints is
        // not null.  The test mNumPoints >= 3 is redundant because of the
        // logic in the constructor.
    }
}

template <typename InputType, typename ComputeType> inline
std::vector<std::array<int, 3>> const&
TriangulateCDT<InputType, ComputeType>::GetTriangles() const
{
    return mTriangles;
}

template <typename InputType, typename ComputeType> inline
std::vector<std::array<int, 3>> const&
TriangulateCDT<InputType, ComputeType>::GetOutsideTriangles() const
{
    return mOutsideTriangles;
}

template <typename InputType, typename ComputeType> inline
std::vector<std::array<int, 3>> const&
TriangulateCDT<InputType, ComputeType>::GetAllTriangles() const
{
    return mAllTriangles;
}

template <typename InputType, typename ComputeType> inline
std::vector<bool> const&
TriangulateCDT<InputType, ComputeType>::GetIsInside() const
{
    return mIsInside;
}

template <typename InputType, typename ComputeType> inline
bool TriangulateCDT<InputType, ComputeType>::IsInside(int triIndex) const
{
    if (0 <= triIndex && triIndex < static_cast<int>(mIsInside.size()))
    {
        return mIsInside[triIndex];
    }
    else
    {
        return false;
    }
}

template <typename InputType, typename ComputeType> inline
bool TriangulateCDT<InputType, ComputeType>::IsOutside(int triIndex) const
{
    if (0 <= triIndex && triIndex < static_cast<int>(mIsInside.size()))
    {
        return !mIsInside[triIndex];
    }
    else
    {
        return false;
    }
}

template <typename InputType, typename ComputeType>
bool TriangulateCDT<InputType, ComputeType>::operator()()
{
    if (mPoints)
    {
        std::vector<int> indices(mNumPoints);
        for (int i = 0; i < mNumPoints; ++i)
        {
            indices[i] = i;
        }

        Tree tree;
        tree.polygon.numIndices = mNumPoints;
        tree.polygon.indices = &indices[0];

        return operator()(tree);
    }
    return false;
}

template <typename InputType, typename ComputeType>
bool TriangulateCDT<InputType, ComputeType>::operator()(
    Polygon const& polygon)
{
    if (mPoints)
    {
        Tree tree;
        tree.polygon = polygon;

        return operator()(tree);
    }
    return false;
}

template <typename InputType, typename ComputeType>
bool TriangulateCDT<InputType, ComputeType>::operator()(Polygon const& outer,
    Polygon const& inner)
{
    if (mPoints)
    {
        Tree otree;
        otree.polygon = outer;
        otree.child.resize(1);

        Tree itree;
        itree.polygon = inner;
        otree.child[0] = &itree;

        return operator()(otree);
    }
    return false;
}

template <typename InputType, typename ComputeType>
bool TriangulateCDT<InputType, ComputeType>::operator()(
    Polygon const& outer, std::vector<Polygon> const& inners)
{
    if (mPoints)
    {
        Tree otree;
        otree.polygon = outer;
        otree.child.resize(inners.size());

        std::vector<Tree> itree(inners.size());
        for (size_t i = 0; i < inners.size(); ++i)
        {
            itree[i].polygon = inners[i];
            otree.child[i] = &itree[i];
        }

        return operator()(otree);
    }
    return false;
}

template <typename InputType, typename ComputeType>
bool TriangulateCDT<InputType, ComputeType>::operator()(Tree const& tree)
{
    if (mPoints)
    {
        std::map<Tree const*, int> offsets;
        int numPoints = GetNumPointsAndOffsets(tree, offsets);
        std::vector<Vector2<InputType>> points(numPoints);
        PackPoints(tree, points);

        if (TriangulatePacked(numPoints, &points[0], tree, offsets))
        {
            int numTriangles = static_cast<int>(mAllTriangles.size());
            for (int t = 0; t < numTriangles; ++t)
            {
                auto& tri = mAllTriangles[t];
                for (int j = 0; j < 3; ++j)
                {
                    LookupIndex(tree, tri[j], offsets);
                }

                if (mIsInside[t])
                {
                    mTriangles.push_back(tri);
                }
                else
                {
                    mOutsideTriangles.push_back(tri);
                }
            }
            return true;
        }
    }

    return false;
}

template <typename InputType, typename ComputeType>
void TriangulateCDT<InputType, ComputeType>::Delete(Tree*& tree,
    bool deletePolygons)
{
    if (tree)
    {
        std::queue<Tree*> treeQueue;
        treeQueue.push(tree);
        while (treeQueue.size() > 0)
        {
            Tree* front = treeQueue.front();
            treeQueue.pop();

            if (deletePolygons)
            {
                delete[] front->polygon.indices;
            }

            int numChildren = static_cast<int>(front->child.size());
            for (int c = 0; c < numChildren; ++c)
            {
                treeQueue.push(front->child[c]);
            }
            delete front;
        }
        tree = nullptr;
    }
}

template <typename InputType, typename ComputeType>
bool TriangulateCDT<InputType, ComputeType>::TriangulatePacked(int numPoints,
    Vector2<InputType> const* points, Tree const& tree,
    std::map<Tree const*, int> const& offsets)
{
    mTriangles.clear();
    mOutsideTriangles.clear();
    mAllTriangles.clear();
    mIsInside.clear();

    mConstrainedDelaunay(numPoints, points, (InputType)0);
    InsertEdges(tree);

    ComputeType half = (ComputeType)0.5;
    ComputeType fourth = (ComputeType)0.25;
    auto const& query = mConstrainedDelaunay.GetQuery();
    auto const* ctPoints = query.GetVertices();
    int numTriangles = mConstrainedDelaunay.GetNumTriangles();
    int const* indices = &mConstrainedDelaunay.GetIndices()[0];
    mIsInside.resize(numTriangles);
    for (int t = 0; t < numTriangles; ++t)
    {
        int v0 = *indices++;
        int v1 = *indices++;
        int v2 = *indices++;
        auto ctInside = fourth * ctPoints[v0] + half * ctPoints[v1] +
            fourth *ctPoints[v2];
        mIsInside[t] = IsInside(tree, ctPoints, ctInside, offsets);
        mAllTriangles.push_back({ { v0, v1, v2 } });
    }
    return true;
}

template <typename InputType, typename ComputeType>
int TriangulateCDT<InputType, ComputeType>::GetNumPointsAndOffsets(
    Tree const& tree, std::map<Tree const*, int>& offsets) const
{
    int numPoints = 0;
    std::queue<const Tree*> treeQueue;
    treeQueue.push(&tree);
    while (treeQueue.size() > 0)
    {
        Tree const* outer = treeQueue.front();
        treeQueue.pop();
        offsets.insert(std::make_pair(outer, numPoints));
        numPoints += outer->polygon.numIndices;

        int numChildren = static_cast<int>(outer->child.size());
        for (int c = 0; c < numChildren; ++c)
        {
            Tree const* inner = outer->child[c];
            offsets.insert(std::make_pair(inner, numPoints));
            numPoints += inner->polygon.numIndices;

            int numGrandChildren = static_cast<int>(inner->child.size());
            for (int g = 0; g < numGrandChildren; ++g)
            {
                treeQueue.push(inner->child[g]);
            }
        }
    }
    return numPoints;
}

template <typename InputType, typename ComputeType>
void TriangulateCDT<InputType, ComputeType>::PackPoints(Tree const& tree,
    std::vector<Vector2<InputType>>& points)
{
    int numPoints = 0;
    std::queue<const Tree*> treeQueue;
    treeQueue.push(&tree);
    while (treeQueue.size() > 0)
    {
        Tree const* outer = treeQueue.front();
        treeQueue.pop();
        for (int i = 0; i < outer->polygon.numIndices; ++i)
        {
            points[numPoints++] = mPoints[outer->polygon.indices[i]];
        }

        int numChildren = static_cast<int>(outer->child.size());
        for (int c = 0; c < numChildren; ++c)
        {
            Tree const* inner = outer->child[c];
            for (int i = 0; i < inner->polygon.numIndices; ++i)
            {
                points[numPoints++] = mPoints[inner->polygon.indices[i]];
            }

            int numGrandChildren = static_cast<int>(inner->child.size());
            for (int g = 0; g < numGrandChildren; ++g)
            {
                treeQueue.push(inner->child[g]);
            }
        }
    }
}

template <typename InputType, typename ComputeType>
bool TriangulateCDT<InputType, ComputeType>::InsertEdges(Tree const& tree)
{
    int numPoints = 0;
    std::array<int, 2> edge;
    std::vector<int> ignoreOutEdge;
    std::queue<const Tree*> treeQueue;
    treeQueue.push(&tree);
    while (treeQueue.size() > 0)
    {
        Tree const* outer = treeQueue.front();
        treeQueue.pop();
        int numOuter = outer->polygon.numIndices;
        for (int i0 = numOuter - 1, i1 = 0; i1 < numOuter; i0 = i1++)
        {
            edge[0] = numPoints + i0;
            edge[1] = numPoints + i1;
            if (!mConstrainedDelaunay.Insert(edge, ignoreOutEdge))
            {
                return false;
            }
        }
        numPoints += numOuter;

        int numChildren = static_cast<int>(outer->child.size());
        for (int c = 0; c < numChildren; ++c)
        {
            Tree const* inner = outer->child[c];
            int numInner = inner->polygon.numIndices;
            for (int i0 = numInner - 1, i1 = 0; i1 < numInner; i0 = i1++)
            {
                edge[0] = numPoints + i0;
                edge[1] = numPoints + i1;
                if (!mConstrainedDelaunay.Insert(edge, ignoreOutEdge))
                {
                    return false;
                }
            }
            numPoints += numInner;

            int numGrandChildren = static_cast<int>(inner->child.size());
            for (int g = 0; g < numGrandChildren; ++g)
            {
                treeQueue.push(inner->child[g]);
            }
        }
    }
    return true;
}

template <typename InputType, typename ComputeType>
void TriangulateCDT<InputType, ComputeType>::LookupIndex(Tree const& tree,
    int& v, std::map<Tree const*, int> const& offsets) const
{
    std::queue<const Tree*> treeQueue;
    treeQueue.push(&tree);
    while (treeQueue.size() > 0)
    {
        Tree const* outer = treeQueue.front();
        treeQueue.pop();
        int offset = offsets.find(outer)->second;
        if (v < offset + outer->polygon.numIndices)
        {
            v = outer->polygon.indices[v - offset];
            return;
        }

        int numChildren = static_cast<int>(outer->child.size());
        for (int c = 0; c < numChildren; ++c)
        {
            Tree const* inner = outer->child[c];
            offset = offsets.find(inner)->second;
            if (v < offset + inner->polygon.numIndices)
            {
                v = inner->polygon.indices[v - offset];
                return;
            }

            int numGrandChildren = static_cast<int>(inner->child.size());
            for (int g = 0; g < numGrandChildren; ++g)
            {
                treeQueue.push(inner->child[g]);
            }
        }
    }
}

template <typename InputType, typename ComputeType>
bool TriangulateCDT<InputType, ComputeType>::IsInside(Tree const& tree,
    Vector2<ComputeType> const* points, Vector2<ComputeType> const& test,
    std::map<Tree const*, int> const& offsets) const
{
    std::queue<const Tree*> treeQueue;
    treeQueue.push(&tree);
    while (treeQueue.size() > 0)
    {
        Tree const* outer = treeQueue.front();
        treeQueue.pop();
        int offset = offsets.find(outer)->second;
        PointInPolygon2<ComputeType> piOuter(outer->polygon.numIndices,
            points + offset);
        if (piOuter.Contains(test))
        {
            int numChildren = static_cast<int>(outer->child.size());
            int c;
            for (c = 0; c < numChildren; ++c)
            {
                Tree const* inner = outer->child[c];
                offset = offsets.find(inner)->second;
                PointInPolygon2<ComputeType> piInner(
                    inner->polygon.numIndices, points + offset);
                if (piInner.Contains(test))
                {
                    int numGrandChildren =
                        static_cast<int>(inner->child.size());
                    for (int g = 0; g < numGrandChildren; ++g)
                    {
                        treeQueue.push(inner->child[g]);
                    }
                    break;
                }
            }
            if (c == numChildren)
            {
                return true;
            }
        }
    }
    return false;
}


}
