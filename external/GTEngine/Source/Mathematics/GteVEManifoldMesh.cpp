// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#include <GTEnginePCH.h>
#include <LowLevel/GteLogger.h>
#include <Mathematics/GteVEManifoldMesh.h>
#include <fstream>
using namespace gte;


VEManifoldMesh::~VEManifoldMesh()
{
    for (auto& element : mVMap)
    {
        delete element.second;
    }

    for (auto& element : mEMap)
    {
        delete element.second;
    }
}

VEManifoldMesh::VEManifoldMesh(VCreator vCreator, ECreator eCreator)
    :
    mVCreator(vCreator ? vCreator : CreateVertex),
    mECreator(eCreator ? eCreator : CreateEdge),
    mAssertOnNonmanifoldInsertion(true)
{
}

VEManifoldMesh::VMap const& VEManifoldMesh::GetVertices() const
{
    return mVMap;
}

VEManifoldMesh::EMap const& VEManifoldMesh::GetEdges() const
{
    return mEMap;
}

VEManifoldMesh::Vertex* VEManifoldMesh::CreateVertex(int v)
{
    return new Vertex(v);
}

VEManifoldMesh::Edge* VEManifoldMesh::CreateEdge(int v0, int v1)
{
    return new Edge(v0, v1);
}

void VEManifoldMesh::AssertOnNonmanifoldInsertion(bool doAssert)
{
    mAssertOnNonmanifoldInsertion = doAssert;
}

VEManifoldMesh::Edge* VEManifoldMesh::Insert(int v0, int v1)
{
    std::pair<int,int> ekey(v0, v1);
    if (mEMap.find(ekey) != mEMap.end())
    {
        // The edge already exists.  Return a null pointer as a signal to
        // the caller that the insertion failed.
        return nullptr;
    }

    // Add the new edge.
    Edge* edge = mECreator(v0,v1);
    mEMap[ekey] = edge;

    // Add the vertices if they do not already exist.
    for (int i = 0; i < 2; ++i)
    {
        int v = edge->V[i];
        Vertex* vertex;
        auto viter = mVMap.find(v);
        if (viter == mVMap.end())
        {
            // This is the first time the vertex is encountered.
            vertex = mVCreator(v);
            mVMap[v] = vertex;

            // Update the vertex.
            vertex->E[0] = edge;
        }
        else
        {
            // This is the second time the vertex is encountered.
            vertex = viter->second;
            if (!vertex)
            {
                LogError("Unexpected condition.");
                return nullptr;
            }

            // Update the vertex.
            if (vertex->E[1])
            {
                if (mAssertOnNonmanifoldInsertion)
                {
                    LogInformation("The mesh must be manifold.");
                }
                return nullptr;
            }
            vertex->E[1] = edge;

            // Update the adjacent edge.
            Edge* adjacent = vertex->E[0];
            if (!adjacent)
            {
                LogError("Unexpected condition.");
                return nullptr;
            }
            for (int j = 0; j < 2; ++j)
            {
                if (adjacent->V[j] == v)
                {
                    adjacent->E[j] = edge;
                    break;
                }
            }

            // Update the edge.
            edge->E[i] = adjacent;
        }
    }

    return edge;
}

bool VEManifoldMesh::Remove(int v0, int v1)
{
    std::pair<int,int> ekey(v0, v1);
    auto eiter = mEMap.find(ekey);
    if (eiter == mEMap.end())
    {
        // The edge does not exist.
        return false;
    }

    // Get the edge.
    Edge* edge = eiter->second;

    // Remove the vertices if necessary (when they are not shared).
    for (int i = 0; i < 2; ++i)
    {
        // Inform the vertices the edge is being deleted.
        auto viter = mVMap.find(edge->V[i]);
        if (viter == mVMap.end())
        {
            // The edge vertices should be in the map.
            LogError("Unexpected condition.");
            return false;
        }

        Vertex* vertex = viter->second;
        if (!vertex)
        {
            // Any <vkey,vertex> in the map must have a dynamically allocated
            // Vertex object.
            LogError("Unexpected condition.");
            return false;
        }
        if (vertex->E[0] == edge)
        {
            // One-edge vertices always have pointer at index zero.
            vertex->E[0] = vertex->E[1];
            vertex->E[1] = nullptr;
        }
        else if (vertex->E[1] == edge)
        {
            vertex->E[1] = nullptr;
        }
        else
        {
            LogError("Unexpected condition.");
            return false;
        }

        // Remove the vertex if you have the last reference to it.
        if (!vertex->E[0] && !vertex->E[1])
        {
            mVMap.erase(vertex->V);
            delete vertex;
        }

        // Inform adjacent edges the edge is being deleted.
        Edge* adjacent = edge->E[i];
        if (adjacent)
        {
            for (int j = 0; j < 2; ++j)
            {
                if (adjacent->E[j] == edge)
                {
                    adjacent->E[j] = nullptr;
                    break;
                }
            }
        }
    }

    mEMap.erase(ekey);
    delete edge;
    return true;
}

bool VEManifoldMesh::IsClosed() const
{
    for (auto const& element : mVMap)
    {
        Vertex const* vertex = element.second;
        if (!vertex->E[0] || !vertex->E[1])
        {
            return false;
        }
    }
    return true;
}

bool VEManifoldMesh::Print(std::string const& filename) const
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

    // Print the vertices.
    outFile << "vertex quantity = " << mVMap.size() << std::endl;
    for (auto const& element : mVMap)
    {
        Vertex const& vertex = *element.second;
        outFile << 'v' << vertex.V << " <";
        if (vertex.E[0])
        {
            outFile << 'e' << edgeIndex[vertex.E[0]];
        }
        else
        {
            outFile << '*';
        }
        outFile << ',';
        if (vertex.E[1])
        {
            outFile << 'e' << edgeIndex[vertex.E[1]];
        }
        else
        {
            outFile << '*';
        }
        outFile << '>' << std::endl;
    }

    // Print the edges.
    outFile << "edge quantity = " << mEMap.size() << std::endl;
    for (auto const& element : mEMap)
    {
        Edge const& edge = *element.second;
        outFile << 'e' << edgeIndex[element.second] << " <"
            << 'v' << edge.V[0] << ",v" << edge.V[1] << "; ";
        if (edge.E[0])
        {
            outFile << 'e' << edgeIndex[edge.E[0]];
        }
        else
        {
            outFile << '*';
        }
        outFile << ',';
        if (edge.E[1])
        {
            outFile << 'e' << edgeIndex[edge.E[1]];
        }
        else
        {
            outFile << '*';
        }
        outFile << '>' << std::endl;
    }
    outFile << std::endl;
    return true;
}

VEManifoldMesh::Vertex::~Vertex()
{
}

VEManifoldMesh::Vertex::Vertex(int v)
{
    V = v;
    E[0] = nullptr;
    E[1] = nullptr;
}

VEManifoldMesh::Edge::~Edge()
{
}

VEManifoldMesh::Edge::Edge(int v0, int v1)
{
    V[0] = v0;
    V[1] = v1;
    E[0] = nullptr;
    E[1] = nullptr;
}

