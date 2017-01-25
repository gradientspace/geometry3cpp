// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#include <GTEnginePCH.h>
#include <LowLevel/GteLogger.h>
#include <Mathematics/GteETManifoldMesh.h>
#include <fstream>
using namespace gte;


ETManifoldMesh::~ETManifoldMesh()
{
    for (auto& element : mEMap)
    {
        delete element.second;
    }

    for (auto& element : mTMap)
    {
        delete element.second;
    }
}

ETManifoldMesh::ETManifoldMesh(ECreator eCreator, TCreator tCreator)
    :
    mECreator(eCreator ? eCreator : CreateEdge),
    mTCreator(tCreator ? tCreator : CreateTriangle),
    mAssertOnNonmanifoldInsertion(true)
{
}

ETManifoldMesh::ETManifoldMesh(ETManifoldMesh const& mesh)
{
    *this = mesh;
}

ETManifoldMesh& ETManifoldMesh::operator=(ETManifoldMesh const& mesh)
{
    Clear();

    mECreator = mesh.mECreator;
    mTCreator = mesh.mTCreator;
    mAssertOnNonmanifoldInsertion = mesh.mAssertOnNonmanifoldInsertion;
    for (auto const& element : mesh.mTMap)
    {
        Insert(element.first.V[0], element.first.V[1], element.first.V[2]);
    }

    return *this;
}

ETManifoldMesh::EMap const& ETManifoldMesh::GetEdges() const
{
    return mEMap;
}

ETManifoldMesh::TMap const& ETManifoldMesh::GetTriangles() const
{
    return mTMap;
}

ETManifoldMesh::Edge* ETManifoldMesh::CreateEdge(int v0, int v1)
{
    return new Edge(v0, v1);
}

ETManifoldMesh::Triangle* ETManifoldMesh::CreateTriangle(int v0, int v1,
    int v2)
{
    return new Triangle(v0, v1, v2);
}

void ETManifoldMesh::AssertOnNonmanifoldInsertion(bool doAssert)
{
    mAssertOnNonmanifoldInsertion = doAssert;
}

ETManifoldMesh::Triangle* ETManifoldMesh::Insert(int v0, int v1, int v2)
{
    TriangleKey<true> tkey(v0, v1, v2);
    if (mTMap.find(tkey) != mTMap.end())
    {
        // The triangle already exists.  Return a null pointer as a signal to
        // the caller that the insertion failed.
        return nullptr;
    }

    // Add the new triangle.
    Triangle* tri = mTCreator(v0, v1, v2);
    mTMap[tkey] = tri;

    // Add the edges to the mesh if they do not already exist.
    for (int i0 = 2, i1 = 0; i1 < 3; i0 = i1++)
    {
        EdgeKey<false> ekey(tri->V[i0], tri->V[i1]);
        Edge* edge;
        auto eiter = mEMap.find(ekey);
        if (eiter == mEMap.end())
        {
            // This is the first time the edge is encountered.
            edge = mECreator(tri->V[i0], tri->V[i1]);
            mEMap[ekey] = edge;

            // Update the edge and triangle.
            edge->T[0] = tri;
            tri->E[i0] = edge;
        }
        else
        {
            // This is the second time the edge is encountered.
            edge = eiter->second;
            if (!edge)
            {
                LogError("Unexpected condition.");
                return nullptr;
            }

            // Update the edge.
            if (edge->T[1])
            {
                if (mAssertOnNonmanifoldInsertion)
                {
                    LogInformation("The mesh must be manifold.");
                }
                return nullptr;
            }
            edge->T[1] = tri;

            // Update the adjacent triangles.
            Triangle* adjacent = edge->T[0];
            if (!adjacent)
            {
                LogError("Unexpected condition.");
                return nullptr;
            }
            for (int j = 0; j < 3; ++j)
            {
                if (adjacent->E[j] == edge)
                {
                    adjacent->T[j] = tri;
                    break;
                }
            }

            // Update the triangle.
            tri->E[i0] = edge;
            tri->T[i0] = adjacent;
        }
    }

    return tri;
}

bool ETManifoldMesh::Remove(int v0, int v1, int v2)
{
    TriangleKey<true> tkey(v0, v1, v2);
    auto titer = mTMap.find(tkey);
    if (titer == mTMap.end())
    {
        // The triangle does not exist.
        return false;
    }

    // Get the triangle.
    Triangle* tri = titer->second;

    // Remove the edges and update adjacent triangles if necessary.
    for (int i = 0; i < 3; ++i)
    {
        // Inform the edges the triangle is being deleted.
        Edge* edge = tri->E[i];
        if (!edge)
        {
            // The triangle edge should be nonnull.
            LogError("Unexpected condition.");
            return false;
        }

        if (edge->T[0] == tri)
        {
            // One-triangle edges always have pointer at index zero.
            edge->T[0] = edge->T[1];
            edge->T[1] = nullptr;
        }
        else if (edge->T[1] == tri)
        {
            edge->T[1] = nullptr;
        }
        else
        {
            LogError("Unexpected condition.");
            return false;
        }

        // Remove the edge if you have the last reference to it.
        if (!edge->T[0] && !edge->T[1])
        {
            EdgeKey<false> ekey(edge->V[0], edge->V[1]);
            mEMap.erase(ekey);
            delete edge;
        }

        // Inform adjacent triangles the triangle is being deleted.
        Triangle* adjacent = tri->T[i];
        if (adjacent)
        {
            for (int j = 0; j < 3; ++j)
            {
                if (adjacent->T[j] == tri)
                {
                    adjacent->T[j] = nullptr;
                    break;
                }
            }
        }
    }

    mTMap.erase(tkey);
    delete tri;
    return true;
}

void ETManifoldMesh::Clear()
{
    for (auto& element : mEMap)
    {
        delete element.second;
    }

    for (auto& element : mTMap)
    {
        delete element.second;
    }

    mEMap.clear();
    mTMap.clear();
}

bool ETManifoldMesh::IsClosed() const
{
    for (auto const& element : mEMap)
    {
        Edge const* edge = element.second;
        if (!edge->T[0] || !edge->T[1])
        {
            return false;
        }
    }
    return true;
}

bool ETManifoldMesh::IsOriented() const
{
    for (auto const& element : mEMap)
    {
        Edge const* edge = element.second;
        if (edge->T[0] && edge->T[1])
        {
            // In each triangle, find the ordered edge that corresponds to the
            // unordered edge element.first.  Also find the vertex opposite
            // that edge.
            bool edgePositive[2] = { false, false };
            int vOpposite[2] = { -1, -1 };
            for (int j = 0; j < 2; ++j)
            {
                for (int i = 0; i < 3; ++i)
                {
                    if (edge->T[j]->V[i] == element.first.V[0])
                    {
                        int vNext = edge->T[j]->V[(i + 1) % 3];
                        if (vNext == element.first.V[1])
                        {
                            edgePositive[j] = true;
                            vOpposite[j] = edge->T[j]->V[(i + 2) % 3];
                        }
                        else
                        {
                            edgePositive[j] = false;
                            vOpposite[j] = vNext;
                        }
                        break;
                    }
                }
            }

            // To be oriented consistently, the edges must have reversed
            // ordering and the oppositive vertices cannot match.
            if (edgePositive[0] == edgePositive[1]
                || vOpposite[0] == vOpposite[1])
            {
                return false;
            }
        }
    }
    return true;
}

void ETManifoldMesh::GetComponents(
    std::vector<std::vector<Triangle const*>>& components) const
{
    // visited: 0 (unvisited), 1 (discovered), 2 (finished)
    std::map<Triangle const*, int> visited;
    for (auto const& element : mTMap)
    {
        visited.insert(std::make_pair(element.second, 0));
    }

    for (auto& element : mTMap)
    {
        Triangle const* tri = element.second;
        if (visited[tri] == 0)
        {
            std::vector<Triangle const*> component;
            DepthFirstSearch(tri, visited, component);
            components.push_back(component);
        }
    }
}

void ETManifoldMesh::GetComponents(
    std::vector<std::vector<TriangleKey<true>>>& components) const
{
    // visited: 0 (unvisited), 1 (discovered), 2 (finished)
    std::map<Triangle const*, int> visited;
    for (auto const& element : mTMap)
    {
        visited.insert(std::make_pair(element.second, 0));
    }

    for (auto& element : mTMap)
    {
        Triangle const* tri = element.second;
        if (visited[tri] == 0)
        {
            std::vector<Triangle const*> component;
            DepthFirstSearch(tri, visited, component);

            std::vector<TriangleKey<true>> keyComponent;
            keyComponent.reserve(component.size());
            for (auto const* t : component)
            {
                keyComponent.push_back(
                    TriangleKey<true>(t->V[0], t->V[1], t->V[2]));
            }
            components.push_back(keyComponent);
        }
    }
}

bool ETManifoldMesh::Print(std::string const& filename)
{
    std::ofstream outFile(filename);
    if (!outFile)
    {
        return false;
    }

    // Assign unique indices to the edges.
    std::map<Edge*,int> edgeIndex;
    edgeIndex[nullptr] = 0;
    int i = 1;
    for (auto const& element : mEMap)
    {
        if (element.second)
        {
            edgeIndex[element.second] = i++;
        }
    }

    // Assign unique indices to the triangles.
    std::map<Triangle*,int> triIndex;
    triIndex[nullptr] = 0;
    i = 1;
    for (auto const& element : mTMap)
    {
        if (element.second)
        {
            triIndex[element.second] = i++;
        }
    }

    // Print the edges.
    outFile << "edge quantity = " << mEMap.size() << std::endl;
    for (auto const& element : mEMap)
    {
        Edge const& edge = *element.second;
        outFile << 'e' << edgeIndex[element.second] << " <"
              << 'v' << edge.V[0] << ",v" << edge.V[1] << "; ";
        for (int j = 0; j < 2; ++j)
        {
            if (edge.T[j])
            {
                outFile << 't' << triIndex[edge.T[0]];
            }
            else
            {
                outFile << '*';
            }
            outFile << (j == 0 ? ',' : '>');
        }
        outFile << std::endl;
    }
    outFile << std::endl;

    // Print the triangles.
    outFile << "triangle quantity = " << mTMap.size() << std::endl;
    for (auto const& element : mTMap)
    {
        Triangle const& tri = *element.second;
        outFile << 't' << triIndex[element.second] << " <"
              << 'v' << tri.V[0] << ",v" << tri.V[1] << ",v"
              << tri.V[2] << "; ";
        for (int j = 0; j < 3; ++j)
        {
            if (tri.E[j])
            {
                outFile << 'e' << edgeIndex[tri.E[j]];
            }
            else
            {
                outFile << '*';
            }
            outFile << (j < 2 ? "," : "; ");
        }

        for (int j = 0; j < 3; ++j)
        {
            if (tri.T[j])
            {
                outFile << 't' << triIndex[tri.T[j]];
            }
            else
            {
                outFile << '*';
            }
            outFile << (j < 2 ? ',' : '>');
        }
        outFile << std::endl;
    }
    outFile << std::endl;
    return true;
}

void ETManifoldMesh::DepthFirstSearch(Triangle const* tInitial,
    std::map<Triangle const*, int>& visited,
    std::vector<Triangle const*>& component) const
{
    // Allocate the maximum-size stack that can occur in the depth-first
    // search.  The stack is empty when the index top is -1.
    std::vector<Triangle const*> tStack(mTMap.size());
    int top = -1;
    tStack[++top] = tInitial;
    while (top >= 0)
    {
        Triangle const* tri = tStack[top];
        visited[tri] = 1;
        int i;
        for (i = 0; i < 3; ++i)
        {
            Triangle const* adj = tri->T[i];
            if (adj && visited[adj] == 0)
            {
                tStack[++top] = adj;
                break;
            }
        }
        if (i == 3)
        {
            visited[tri] = 2;
            component.push_back(tri);
            --top;
        }
    }
}

ETManifoldMesh::Edge::~Edge()
{
}

ETManifoldMesh::Edge::Edge(int v0, int v1)
{
    V[0] = v0;
    V[1] = v1;
    T[0] = nullptr;
    T[1] = nullptr;
}

ETManifoldMesh::Triangle::~Triangle()
{
}

ETManifoldMesh::Triangle::Triangle(int v0, int v1, int v2)
{
    V[0] = v0;
    V[1] = v1;
    V[2] = v2;
    for (int i = 0; i < 3; ++i)
    {
        E[i] = nullptr;
        T[i] = nullptr;
    }
}

