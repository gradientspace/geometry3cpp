// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#include <GTEnginePCH.h>
#include <LowLevel/GteLogger.h>
#include <Mathematics/GteTSManifoldMesh.h>
#include <fstream>
using namespace gte;


TSManifoldMesh::~TSManifoldMesh()
{
    for (auto& element : mTMap)
    {
        delete element.second;
    }

    for (auto& element : mSMap)
    {
        delete element.second;
    }
}

TSManifoldMesh::TSManifoldMesh(TCreator tCreator, SCreator sCreator)
    :
    mTCreator(tCreator ? tCreator : CreateTriangle),
    mSCreator(sCreator ? sCreator : CreateTetrahedron),
    mAssertOnNonmanifoldInsertion(true)
{
}

TSManifoldMesh::TMap const& TSManifoldMesh::GetTriangles() const
{
    return mTMap;
}

TSManifoldMesh::SMap const& TSManifoldMesh::GetTetrahedra() const
{
    return mSMap;
}

TSManifoldMesh::Triangle* TSManifoldMesh::CreateTriangle(int v0, int v1,
    int v2)
{
    return new Triangle(v0, v1, v2);
}

TSManifoldMesh::Tetrahedron* TSManifoldMesh::CreateTetrahedron(int v0,
    int v1, int v2, int v3)
{
    return new Tetrahedron(v0, v1, v2, v3);
}

void TSManifoldMesh::AssertOnNonmanifoldInsertion(bool doAssert)
{
    mAssertOnNonmanifoldInsertion = doAssert;
}

TSManifoldMesh::Tetrahedron* TSManifoldMesh::Insert(int v0, int v1, int v2,
    int v3)
{
    TetrahedronKey<true> skey(v0, v1, v2, v3);
    if (mSMap.find(skey) != mSMap.end())
    {
        // The tetrahedron already exists.  Return a null pointer as a signal
        // to the caller that the insertion failed.
        return nullptr;
    }

    // Add the new tetrahedron.
    Tetrahedron* tetra = mSCreator(v0, v1, v2, v3);
    mSMap[skey] = tetra;

    // Add the faces to the mesh if they do not already exist.
    for (int i = 0; i < 4; ++i)
    {
        auto opposite = TetrahedronKey<true>::oppositeFace[i];
        TriangleKey<false> tkey(
            tetra->V[opposite[0]],
            tetra->V[opposite[1]],
            tetra->V[opposite[2]]);
        Triangle* face;
        auto titer = mTMap.find(tkey);
        if (titer == mTMap.end())
        {
            // This is the first time the face is encountered.
            face = mTCreator(
                tetra->V[opposite[0]],
                tetra->V[opposite[1]],
                tetra->V[opposite[2]]);
            mTMap[tkey] = face;

            // Update the face and tetrahedron.
            face->T[0] = tetra;
            tetra->T[i] = face;
        }
        else
        {
            // This is the second time the face is encountered.
            face = titer->second;
            if (!face)
            {
                LogError("Unexpected condition.");
                return nullptr;
            }

            // Update the face.
            if (face->T[1])
            {
                if (mAssertOnNonmanifoldInsertion)
                {
                    LogInformation("The mesh must be manifold.");
                }
                return nullptr;
            }
            face->T[1] = tetra;

            // Update the adjacent triangles.
            Tetrahedron* adjacent = face->T[0];
            if (!adjacent)
            {
                LogError("Unexpected condition.");
                return nullptr;
            }
            for (int j = 0; j < 4; ++j)
            {
                if (adjacent->T[j] == face)
                {
                    adjacent->S[j] = tetra;
                    break;
                }
            }

            // Update the tetrahedron.
            tetra->T[i] = face;
            tetra->S[i] = adjacent;
        }
    }

    return tetra;
}

bool TSManifoldMesh::Remove(int v0, int v1, int v2, int v3)
{
    TetrahedronKey<true> skey(v0, v1, v2, v3);
    auto siter = mSMap.find(skey);
    if (siter == mSMap.end())
    {
        // The tetrahedron does not exist.
        return false;
    }

    // Get the tetrahedron.
    Tetrahedron* tetra = siter->second;

    // Remove the faces and update adjacent tetrahedra if necessary.
    for (int i = 0; i < 4; ++i)
    {
        // Inform the faces the tetrahedron is being deleted.
        Triangle* face = tetra->T[i];
        if (!face)
        {
            // The triangle edge should be nonnull.
            LogError("Unexpected condition.");
            return false;
        }

        if (face->T[0] == tetra)
        {
            // One-tetrahedron faces always have pointer at index zero.
            face->T[0] = face->T[1];
            face->T[1] = nullptr;
        }
        else if (face->T[1] == tetra)
        {
            face->T[1] = nullptr;
        }
        else
        {
            LogError("Unexpected condition.");
            return false;
        }

        // Remove the face if you have the last reference to it.
        if (!face->T[0] && !face->T[1])
        {
            TriangleKey<false> tkey(face->V[0], face->V[1], face->V[2]);
            mTMap.erase(tkey);
            delete face;
        }

        // Inform adjacent tetrahedra the tetrahedron is being deleted.
        Tetrahedron* adjacent = tetra->S[i];
        if (adjacent)
        {
            for (int j = 0; j < 4; ++j)
            {
                if (adjacent->S[j] == tetra)
                {
                    adjacent->S[j] = nullptr;
                    break;
                }
            }
        }
    }

    mSMap.erase(skey);
    delete tetra;
    return true;
}

bool TSManifoldMesh::IsClosed() const
{
    for (auto const& element : mTMap)
    {
        Triangle const* tri = element.second;
        if (!tri->T[0] || !tri->T[1])
        {
            return false;
        }
    }
    return true;
}

bool TSManifoldMesh::Print(std::string const& filename)
{
    std::ofstream outFile(filename);
    if (!outFile)
    {
        return false;
    }

    // Assign unique indices to the triangles.
    std::map<Triangle*,int> triIndex;
    triIndex[nullptr] = 0;
    int i = 1;
    for (auto const& element : mTMap)
    {
        if (element.second)
        {
            triIndex[element.second] = i++;
        }
    }

    // Assign unique indices to the tetrahedra.
    std::map<Tetrahedron*, int> tetraIndex;
    tetraIndex[nullptr] = 0;
    i = 1;
    for (auto const& element : mSMap)
    {
        if (element.second)
        {
            tetraIndex[element.second] = i++;
        }
    }

    // Print the triangles.
    outFile << "triangle quantity = " << mTMap.size() << std::endl;
    for (auto const& element : mTMap)
    {
        Triangle const& tri = *element.second;
        outFile << 't' << triIndex[element.second] << " <"
            << 'v' << tri.V[0] << ",v" << tri.V[1] << ",v"
            << tri.V[2] << "; ";
        for (int j = 0; j < 2; ++j)
        {
            if (tri.T[j])
            {
                outFile << 's' << tetraIndex[tri.T[j]];
            }
            else
            {
                outFile << '*';
            }
            outFile << (j == 0 ? '*' : '>');
        }
        outFile << std::endl;
    }
    outFile << std::endl;

    // Print the tetrahedra.
    outFile << "tetrahedron quantity = " << mSMap.size() << std::endl;
    for (auto const& element : mSMap)
    {
        Tetrahedron const& tetra = *element.second;
        outFile << 's' << tetraIndex[element.second] << " <"
            << 'v' << tetra.V[0] << ",v" << tetra.V[1] << ",v"
            << tetra.V[2] << ",v" << tetra.V[3] << "; ";
        for (int j = 0; j < 4; ++j)
        {
            if (tetra.T[j])
            {
                outFile << 't' << triIndex[tetra.T[j]];
            }
            else
            {
                outFile << '*';
            }
            outFile << (j < 3 ? "," : "; ");
        }

        for (int j = 0; j < 4; ++j)
        {
            if (tetra.S[j])
            {
                outFile << 't' << tetraIndex[tetra.S[j]];
            }
            else
            {
                outFile << '*';
            }
            outFile << (j < 3 ? ',' : '>');
        }

        outFile << std::endl;
    }
    outFile << std::endl;
    return true;
}

TSManifoldMesh::Triangle::~Triangle()
{
}

TSManifoldMesh::Triangle::Triangle(int v0, int v1, int v2)
{
    V[0] = v0;
    V[1] = v1;
    V[2] = v2;
    T[0] = nullptr;
    T[1] = nullptr;
}

TSManifoldMesh::Tetrahedron::~Tetrahedron()
{
}

TSManifoldMesh::Tetrahedron::Tetrahedron(int v0, int v1, int v2, int v3)
{
    V[0] = v0;
    V[1] = v1;
    V[2] = v2;
    V[3] = v3;
    for (int i = 0; i < 4; ++i)
    {
        T[i] = nullptr;
        S[i] = nullptr;
    }
}

