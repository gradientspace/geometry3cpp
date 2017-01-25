// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.10.1 (2014/01/21)

#include "Wm5MathematicsPCH.h"
#include "Wm5TSManifoldMesh.h"
#include "Wm5Memory.h"
using namespace Wm5;

//----------------------------------------------------------------------------
TSManifoldMesh::~TSManifoldMesh()
{
    TMapIterator telement;
    for (telement = mTMap.begin(); telement != mTMap.end(); ++telement)
    {
        delete0(telement->second);
    }

    SMapIterator selement;
    for (selement = mSMap.begin(); selement != mSMap.end(); ++selement)
    {
        delete0(selement->second);
    }
}
//----------------------------------------------------------------------------
TSManifoldMesh::TSManifoldMesh(TCreator tCreator, SCreator sCreator)
{
    mTCreator = (tCreator ? tCreator : CreateTriangle);
    mSCreator = (sCreator ? sCreator : CreateTetrahedron);
}
//----------------------------------------------------------------------------
const TSManifoldMesh::TMap& TSManifoldMesh::GetTriangles() const
{
    return mTMap;
}
//----------------------------------------------------------------------------
const TSManifoldMesh::SMap& TSManifoldMesh::GetTetrahedra() const
{
    return mSMap;
}
//----------------------------------------------------------------------------
TSManifoldMesh::Triangle* TSManifoldMesh::CreateTriangle(int v0, int v1,
    int v2)
{
    return new0 Triangle(v0, v1, v2);
}
//----------------------------------------------------------------------------
TSManifoldMesh::Tetrahedron* TSManifoldMesh::CreateTetrahedron(int v0,
    int v1, int v2, int v3)
{
    return new0 Tetrahedron(v0, v1, v2, v3);
}
//----------------------------------------------------------------------------
TSManifoldMesh::Tetrahedron* TSManifoldMesh::Insert(int v0, int v1, int v2,
    int v3)
{
    TetrahedronKey skey(v0, v1, v2, v3);
    if (mSMap.find(skey) != mSMap.end())
    {
        // The tetrahedron already exists.  Return a null pointer as a signal
        // to the caller that the insertion failed.
        return 0;
    }

    // Add the new tetrahedron.
    Tetrahedron* tetra = mSCreator(v0, v1, v2, v3);
    mSMap[skey] = tetra;

    // Add the faces to the mesh if they do not already exist.
    for (int i = 0; i < 4; ++i)
    {
        int opposite[3] = {
            TetrahedronKey::oppositeFace[i][0],
            TetrahedronKey::oppositeFace[i][1],
            TetrahedronKey::oppositeFace[i][2]};
        UnorderedTriangleKey tkey(
            tetra->V[opposite[0]],
            tetra->V[opposite[1]],
            tetra->V[opposite[2]]);
        Triangle* face;
        TMapIterator titer = mTMap.find(tkey);
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
                assertion(false, "Unexpected condition.");
                return 0;
            }

            // Update the face.
            if (face->T[1])
            {
                assertion(false, "The mesh must be manifold.");
                return 0;
            }
            face->T[1] = tetra;

            // Update the adjacent triangles.
            Tetrahedron* adjacent = face->T[0];
            if (!adjacent)
            {
                assertion(false, "Unexpected condition.");
                return 0;
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
//----------------------------------------------------------------------------
bool TSManifoldMesh::Remove(int v0, int v1, int v2, int v3)
{
    TetrahedronKey skey(v0, v1, v2, v3);
    SMapIterator siter = mSMap.find(skey);
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
            assertion(false, "Unexpected condition.");
            return false;
        }

        if (face->T[0] == tetra)
        {
            // One-tetrahedron faces always have pointer at index zero.
            face->T[0] = face->T[1];
            face->T[1] = 0;
        }
        else if (face->T[1] == tetra)
        {
            face->T[1] = 0;
        }
        else
        {
            assertion(false, "Unexpected condition.");
            return false;
        }

        // Remove the face if you have the last reference to it.
        if (!face->T[0] && !face->T[1])
        {
            UnorderedTriangleKey tkey(face->V[0], face->V[1], face->V[2]);
            mTMap.erase(tkey);
            delete0(face);
        }

        // Inform adjacent tetrahedra the tetrahedron is being deleted.
        Tetrahedron* adjacent = tetra->S[i];
        if (adjacent)
        {
            for (int j = 0; j < 4; ++j)
            {
                if (adjacent->S[j] == tetra)
                {
                    adjacent->S[j] = 0;
                    break;
                }
            }
        }
    }

    mSMap.erase(skey);
    delete0(tetra);
    return true;
}
//----------------------------------------------------------------------------
bool TSManifoldMesh::IsClosed() const
{
    TMapCIterator element;
    for (element = mTMap.begin(); element != mTMap.end(); ++element)
    {
        Triangle const* tri = element->second;
        if (!tri->T[0] || !tri->T[1])
        {
            return false;
        }
    }
    return true;
}
//----------------------------------------------------------------------------
bool TSManifoldMesh::Print(const char* filename)
{
    std::ofstream outFile(filename);
    if (!outFile)
    {
        return false;
    }

    // Assign unique indices to the triangles.
    std::map<Triangle*,int> triIndex;
    triIndex[(Triangle*)0] = 0;
    int i = 1;
    TMapCIterator telement;
    for (telement = mTMap.begin(); telement != mTMap.end(); ++telement)
    {
        if (telement->second)
        {
            triIndex[telement->second] = i++;
        }
    }

    // Assign unique indices to the tetrahedra.
    std::map<Tetrahedron*, int> tetraIndex;
    tetraIndex[(Tetrahedron*)0] = 0;
    i = 1;
    SMapCIterator selement;
    for (selement = mSMap.begin(); selement != mSMap.end(); ++selement)
    {
        if (selement->second)
        {
            tetraIndex[selement->second] = i++;
        }
    }

    // Print the triangles.
    outFile << "triangle quantity = " << mTMap.size() << std::endl;
    for (telement = mTMap.begin(); telement != mTMap.end(); ++telement)
    {
        Triangle const& tri = *telement->second;
        outFile << 't' << triIndex[telement->second] << " <"
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
    for (selement = mSMap.begin(); selement != mSMap.end(); ++selement)
    {
        Tetrahedron const& tetra = *selement->second;
        outFile << 's' << tetraIndex[selement->second] << " <"
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
//----------------------------------------------------------------------------
TSManifoldMesh::Triangle::~Triangle()
{
}
//----------------------------------------------------------------------------
TSManifoldMesh::Triangle::Triangle(int v0, int v1, int v2)
{
    V[0] = v0;
    V[1] = v1;
    V[2] = v2;
    T[0] = 0;
    T[1] = 0;
}
//----------------------------------------------------------------------------
TSManifoldMesh::Tetrahedron::~Tetrahedron()
{
}
//----------------------------------------------------------------------------
TSManifoldMesh::Tetrahedron::Tetrahedron(int v0, int v1, int v2, int v3)
{
    V[0] = v0;
    V[1] = v1;
    V[2] = v2;
    V[3] = v3;
    for (int i = 0; i < 4; ++i)
    {
        T[i] = 0;
        S[i] = 0;
    }
}
//----------------------------------------------------------------------------
