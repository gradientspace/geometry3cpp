// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <Mathematics/GteTetrahedronKey.h>
#include <Mathematics/GteTriangleKey.h>
#include <map>

namespace gte
{

class GTE_IMPEXP TSManifoldMesh
{
public:
    // Triangle data types.
    class Triangle;
    typedef Triangle* (*TCreator)(int, int, int);
    typedef std::map<TriangleKey<false>, Triangle*> TMap;

    // Tetrahedron data types.
    class Tetrahedron;
    typedef Tetrahedron* (*SCreator)(int, int, int, int);
    typedef std::map<TetrahedronKey<true>, Tetrahedron*> SMap;

    // Triangle object.
    class GTE_IMPEXP Triangle
    {
    public:
        virtual ~Triangle();
        Triangle(int v0, int v1, int v2);

        // Vertices of the face.
        int V[3];

        // Tetrahedra sharing the face.
        Tetrahedron* T[2];
    };

    // Tetrahedron object.
    class GTE_IMPEXP Tetrahedron
    {
    public:
        virtual ~Tetrahedron();
        Tetrahedron(int v0, int v1, int v2, int v3);

        // Vertices, listed in an order so that each face vertices in
        // counterclockwise order when viewed from outside the tetrahedron.
        int V[4];

        // Adjacent faces.  T[i] points to the triangle face opposite V[i].
        //   T[0] points to face (V[1],V[2],V[3])
        //   T[1] points to face (V[0],V[3],V[2])
        //   T[2] points to face (V[0],V[1],V[3])
        //   T[3] points to face (V[0],V[2],V[1])
        Triangle* T[4];

        // Adjacent tetrahedra.  S[i] points to the adjacent tetrahedron
        // sharing face T[i].
        Tetrahedron* S[4];
    };


    // Construction and destruction.
    virtual ~TSManifoldMesh();
    TSManifoldMesh(TCreator tCreator = nullptr, SCreator sCreator = nullptr);

    // Member access.
    TMap const& GetTriangles() const;
    SMap const& GetTetrahedra() const;

    // If the insertion of a tetrahedron fails because the mesh would become
    // nonmanifold, the default behavior is to trigger a LogError message.
    // You can disable this behavior in situations where you want the Logger
    // system on but you want to continue gracefully without an assertion.
    void AssertOnNonmanifoldInsertion(bool doAssert);

    // If <v0,v1,v2,v3> is not in the mesh, a Tetrahedron object is created
    // and returned; otherwise, <v0,v1,v2,v3> is in the mesh and nullptr is
    // returned.  If the insertion leads to a nonmanifold mesh, the call
    // fails with a nullptr returned.
    Tetrahedron* Insert(int v0, int v1, int v2, int v3);

    // If <v0,v1,v2,v3> is in the mesh, it is removed and 'true' is returned;
    // otherwise, <v0,v1,v2,v3> is not in the mesh and 'false' is returned.
    bool Remove(int v0, int v1, int v2, int v3);

    // A manifold mesh is closed if each face is shared twice.
    bool IsClosed() const;

    // For debugging.  The function returns 'true' iff the text file has been
    // created and saved.
    bool Print(std::string const& filename);

protected:
    // The triangle data and default triangle creation.
    static Triangle* CreateTriangle(int v0, int v1, int v2);
    TCreator mTCreator;
    TMap mTMap;

    // The tetrahedron data and default tetrahedron creation.
    static Tetrahedron* CreateTetrahedron(int v0, int v1, int v2, int v3);
    SCreator mSCreator;
    SMap mSMap;
    bool mAssertOnNonmanifoldInsertion;  // default: true
};

}
