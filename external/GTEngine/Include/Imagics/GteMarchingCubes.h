// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <GTEngineDEF.h>
#include <array>
#include <string>

namespace gte
{

// Create the lookup table for the Marching Cubes algorithm that is used to
// extract a triangular mesh that represents a level surface of a 3D image
// sampled on a regular lattice.  The assumption is that no data sample is
// zero, which allows us to have a table with 256 entries: 2 signs per
// sample, 8 samples per volume element (voxel).  Each entry corresponds to
// the pattern of 8 signs at the corners of a voxel.  The signs are stored as
// bits (b7,b6,b5,b4,b3,b2,b1,b0).  The bit assignments to voxel corners is
//   b0 = (x,y,z),   b1 = (x+1,y,z),   b2 = (x,y+1,z),   b3 = (x+1,y+1,z)
//   b4 = (x,y,z+1), b5 = (x+1,y,z+1), b6 = (x,y+1,z+1), b7 = (x+1,y+1,z+1)
// If a bit is zero, then the voxel value at the corresponding corner is
// positive; otherwise, the bit is one and the value is negative.  The
// triangles are counterclockwise ordered according to an observer viewing
// the triangle from the negative side of the level surface.

class GTE_IMPEXP MarchingCubes
{
public:
    // Construction and destruction.
    virtual ~MarchingCubes();
    MarchingCubes();

    enum GTE_IMPEXP
    {
        MAX_VERTICES = 12,
        MAX_TRIANGLES = 5,
    };

    struct GTE_IMPEXP Topology
    {
        // All members are set to zeros.
        Topology();

        int numVertices;
        int numTriangles;
        std::array<std::array<int, 2>, MAX_VERTICES> vpair;
        std::array<std::array<int, 3>, MAX_TRIANGLES> itriple;
    };

    // The entry must be in {0..255}.
    inline Topology const& GetTable(int entry) const;

    // The table has 256 entries, each 41 integers, stored as table[256][41].
    // The return value is a pointer to the table via &table[0][0].
    inline int const* GetTable() const;

    // Get the configuration type for the voxel, which is one of the string
    // names of the 'void Bits* (int[8])' functions.
    static std::string GetConfigurationType(int entry);

protected:
    // Support for lookup construction and access.
    // mTable[i][0] = numVertices
    // mTable[i][1] = numTriangles
    // mTable[i][2..25] = pairs of corner indices (maximum of 12 pairs)
    // mTable[i][26..40] = triples of indices (maximum of 5 triples)
    Topology mTable[256];

    // The constructor iterates mEntry from 0 to 255 and calls configuration
    // functions, each calling SetTable(...).  The mEntry value is the table
    // index to be used.
    int mEntry;

    void SetTable(int numV, int const* vpair, int numT, int const* itriple);

    // The precomputed information about the 256 configurations for voxels.
    void Bits0(int index[8]);
    void Bits1(int index[8]);
    void Bits7(int index[8]);
    void Bits2Edge(int index[8]);
    void Bits6Edge(int index[8]);
    void Bits2FaceDiag(int index[8]);
    void Bits6FaceDiag(int index[8]);
    void Bits2BoxDiag(int index[8]);
    void Bits6BoxDiag(int index[8]);
    void Bits3SameFace(int index[8]);
    void Bits5SameFace(int index[8]);
    void Bits3EdgeFaceDiag(int index[8]);
    void Bits5EdgeFaceDiag(int index[8]);
    void Bits3FaceDiagFaceDiag(int index[8]);
    void Bits5FaceDiagFaceDiag(int index[8]);
    void Bits4SameFace(int index[8]);
    void Bits4FaceEdge(int index[8]);
    void Bits4FaceFaceDiagL(int index[8]);
    void Bits4FaceFaceDiagR(int index[8]);
    void Bits4FaceBoxDiag(int index[8]);
    void Bits4EdgeEdgePara(int index[8]);
    void Bits4EdgeEdgePerp(int index[8]);

    enum GTE_IMPEXP ConfigurationType
    {
        CT_Bits0,
        CT_Bits1,
        CT_Bits7,
        CT_Bits2Edge,
        CT_Bits6Edge,
        CT_Bits2FaceDiag,
        CT_Bits6FaceDiag,
        CT_Bits2BoxDiag,
        CT_Bits6BoxDiag,
        CT_Bits3SameFace,
        CT_Bits5SameFace,
        CT_Bits3EdgeFaceDiag,
        CT_Bits5EdgeFaceDiag,
        CT_Bits3FaceDiagFaceDiag,
        CT_Bits5FaceDiagFaceDiag,
        CT_Bits4SameFace,
        CT_Bits4FaceEdge,
        CT_Bits4FaceFaceDiagL,
        CT_Bits4FaceFaceDiagR,
        CT_Bits4FaceBoxDiag,
        CT_Bits4EdgeEdgePara,
        CT_Bits4EdgeEdgePerp,
        CT_NUM_TYPES
    };

    typedef void (MarchingCubes::*Function)(int[8]);

    struct GTE_IMPEXP Configuration
    {
        ConfigurationType type;
        Function F;
        int index[8];
    };

    static Configuration msConfiguration[256];
    static std::string msConfigurationString[CT_NUM_TYPES];
};

inline MarchingCubes::Topology const& MarchingCubes::GetTable(int entry) const
{
    return mTable[entry];
}

int const* MarchingCubes::GetTable() const
{
    return reinterpret_cast<int const*>(mTable);
}

}
