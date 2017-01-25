// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <Imagics/GteMarchingCubes.h>
#include <cstring>

namespace gte
{

template <typename Real>
class SurfaceExtractor : public MarchingCubes
{
public:
    // Construction and destruction.
    virtual ~SurfaceExtractor();
    SurfaceExtractor();

    // The input function values must be stored as
    //  F[0] = function(0,0,0), F[4] = function(0,0,1),
    //  F[1] = function(1,0,0), F[5] = function(1,0,1),
    //  F[2] = function(0,1,0), F[6] = function(0,1,1),
    //  F[3] = function(1,1,0), F[7] = function(1,1,1).
    // Thus, F[k] = function(k & 1, (k & 2) >> 1, (k & 4) >> 2).
    // The return value is 'true' iff the F[] values are all nonzero.
    // If they are not, the returned 'mesh' has no vertices and no
    // triangles--as if F[] had all positive or all negative values.

    struct Mesh
    {
        int numVertices;
        int numTriangles;
        int vpair[2*MAX_VERTICES];
        int itriple[3*MAX_TRIANGLES];
        Real vertices[MAX_VERTICES][3];
    };

    bool Extract(Real const F[8], Mesh& mesh) const;
};


template <typename Real>
SurfaceExtractor<Real>::~SurfaceExtractor()
{
}

template <typename Real>
SurfaceExtractor<Real>::SurfaceExtractor()
{
}

template <typename Real>
bool SurfaceExtractor<Real>::Extract(Real const F[8], Mesh& mesh) const
{
    mesh.numVertices = 0;
    mesh.numTriangles = 0;
    memset(mesh.vpair, 0, 2 * MAX_VERTICES*sizeof(int));
    memset(mesh.itriple, 0, 3 * MAX_TRIANGLES*sizeof(int));
    memset(&mesh.vertices[0][0], 0, 3 * MAX_VERTICES*sizeof(Real));

    int entry = 0;
    for (int i = 0, mask = 1; i < 8; ++i, mask <<= 1)
    {
        if (F[i] < (Real)0)
        {
            entry |= mask;
        }
        else if (F[i] == (Real)0)
        {
            return false;
        }
    }

    Unpack(entry, mesh.numVertices, mesh.vpair, mesh.numTriangles,
        mesh.itriple);

    for (int i = 0; i < mesh.numVertices; ++i)
    {
        int j0 = mesh.vpair[2 * i + 0];
        int j1 = mesh.vpair[2 * i + 1];

        Real corner0[3];
        corner0[0] = static_cast<Real>(j0 & 1);
        corner0[1] = static_cast<Real>((j0 & 2) >> 1);
        corner0[2] = static_cast<Real>((j0 & 4) >> 2);

        Real corner1[3];
        corner1[0] = static_cast<Real>(j1 & 1);
        corner1[1] = static_cast<Real>((j1 & 2) >> 1);
        corner1[2] = static_cast<Real>((j1 & 4) >> 2);

        Real invDenom = ((Real)1) / (F[j0] - F[j1]);
        for (int k = 0; k < 3; ++k)
        {
            Real numer = F[j0] * corner1[k] - F[j1] * corner0[k];
            mesh.vertices[i][k] = numer*invDenom;
        }
    }
    return true;
}


}
