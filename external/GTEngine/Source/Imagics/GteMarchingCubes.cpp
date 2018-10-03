// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#include <GTEnginePCH.h>
#include <LowLevel/GteWrapper.h>
#include <Imagics/GteMarchingCubes.h>
#include <algorithm>
using namespace gte;

MarchingCubes::~MarchingCubes()
{
}

MarchingCubes::MarchingCubes()
{
    // Create the lookup table.
    for (mEntry = 0; mEntry < 256; ++mEntry)
    {
        (this->*msConfiguration[mEntry].F)(msConfiguration[mEntry].index);
    }
}

std::string MarchingCubes::GetConfigurationType(int entry)
{
    if (0 <= entry && entry < 256)
    {
        return msConfigurationString[msConfiguration[entry].type];
    }
    return "";
}

void MarchingCubes::SetTable(int numV, int const* vpair, int numT, int const* itriple)
{
    // The item is already zeroed in the constructor.
    Topology& topology = mTable[mEntry];
    topology.numVertices = numV;
    topology.numTriangles = numT;

    // Store vertex pairs with minimum index occurring first.
    for (int i = 0; i < numV; ++i, vpair += 2)
    {
        topology.vpair[i][0] = std::min(vpair[0], vpair[1]);
        topology.vpair[i][1] = std::max(vpair[0], vpair[1]);
    }

    // Store triangle triples as is.
    for (int i = 0; i < numT; ++i, itriple += 3)
    {
        topology.itriple[i] = { itriple[0], itriple[1], itriple[2] };
    }
}

void MarchingCubes::Bits0(int[8])
{
    SetTable(0, nullptr, 0, nullptr);
}

void MarchingCubes::Bits1(int index[8])
{
    int const numV = 3;
    int vpair[2 * numV] =
    {
        index[1], index[0],
        index[4], index[0],
        index[2], index[0]
    };

    int const numT = 1;
    int itriple[3 * numT] =
    {
        0, 1, 2
    };

    SetTable(numV, vpair, numT, itriple);
}

void MarchingCubes::Bits7(int index[8])
{
    int const numV = 3;
    int vpair[2 * numV] =
    {
        index[1], index[0],
        index[4], index[0],
        index[2], index[0]
    };

    int const numT = 1;
    int itriple[3 * numT] =
    {
        0, 2, 1
    };

    SetTable(numV, vpair, numT, itriple);
}

void MarchingCubes::Bits2Edge(int index[8])
{
    int const numV = 4;
    int vpair[2 * numV] =
    {
        index[4], index[0],
        index[2], index[0],
        index[3], index[1],
        index[5], index[1]
    };

    int const numT = 2;
    int itriple[3 * numT] =
    {
        0, 1, 2,
        0, 2, 3
    };

    SetTable(numV, vpair, numT, itriple);
}

void MarchingCubes::Bits6Edge(int index[8])
{
    int const numV = 4;
    int vpair[2 * numV] =
    {
        index[4], index[0],
        index[2], index[0],
        index[3], index[1],
        index[5], index[1]
    };

    int const numT = 2;
    int itriple[3 * numT] =
    {
        0, 2, 1,
        0, 3, 2
    };

    SetTable(numV, vpair, numT, itriple);
}

void MarchingCubes::Bits2FaceDiag(int index[8])
{
    int const numV = 6;
    int vpair[2 * numV] =
    {
        index[1], index[0],
        index[4], index[0],
        index[2], index[0],
        index[2], index[3],
        index[7], index[3],
        index[1], index[3]
    };

    int const numT = 2;
    int itriple[3 * numT] =
    {
        0, 1, 2,
        3, 4, 5
    };

    SetTable(numV, vpair, numT, itriple);
}

void MarchingCubes::Bits6FaceDiag(int index[8])
{
    int const numV = 6;
    int vpair[2 * numV] =
    {
        index[1], index[0],
        index[4], index[0],
        index[2], index[0],
        index[2], index[3],
        index[7], index[3],
        index[1], index[3]
    };

    // Not the reverse ordering from Bit2FaceDiag due to ambiguous face
    // handling.
    int const numT = 4;
    int itriple[3 * numT] =
    {
        1, 0, 5,
        1, 5, 4,
        1, 4, 3,
        1, 3, 2
    };

    SetTable(numV, vpair, numT, itriple);
}

void MarchingCubes::Bits2BoxDiag(int index[8])
{
    int const numV = 6;
    int vpair[2 * numV] =
    {
        index[1], index[0],
        index[4], index[0],
        index[2], index[0],
        index[3], index[7],
        index[6], index[7],
        index[5], index[7]
    };

    int const numT = 2;
    int itriple[3 * numT] =
    {
        0, 1, 2,
        3, 4, 5
    };

    SetTable(numV, vpair, numT, itriple);
}

void MarchingCubes::Bits6BoxDiag(int index[8])
{
    int const numV = 6;
    int vpair[2 * numV] =
    {
        index[1], index[0],
        index[4], index[0],
        index[2], index[0],
        index[3], index[7],
        index[6], index[7],
        index[5], index[7]
    };

    int const numT = 2;
    int itriple[3 * numT] =
    {
        0, 2, 1,
        3, 5, 4
    };

    SetTable(numV, vpair, numT, itriple);
}

void MarchingCubes::Bits3SameFace(int index[8])
{
    int const numV = 5;
    int vpair[2 * numV] =
    {
        index[4], index[0],
        index[2], index[6],
        index[2], index[3],
        index[1], index[3],
        index[1], index[5]
    };

    int const numT = 3;
    int itriple[3 * numT] =
    {
        0, 1, 2,
        0, 2, 3,
        0, 3, 4
    };

    SetTable(numV, vpair, numT, itriple);
}

void MarchingCubes::Bits5SameFace(int index[8])
{
    int const numV = 5;
    int vpair[2 * numV] =
    {
        index[4], index[0],
        index[2], index[6],
        index[2], index[3],
        index[1], index[3],
        index[1], index[5]
    };

    int const numT = 3;
    int itriple[3 * numT] =
    {
        0, 2, 1,
        0, 3, 2,
        0, 4, 3
    };

    SetTable(numV, vpair, numT, itriple);
}

void MarchingCubes::Bits3EdgeFaceDiag(int index[8])
{
    int const numV = 7;
    int vpair[2 * numV] =
    {
        index[0], index[1],
        index[4], index[5],
        index[4], index[6],
        index[0], index[2],
        index[2], index[3],
        index[3], index[7],
        index[1], index[3]
    };

    int const numT = 3;
    int itriple[3 * numT] =
    {
        0, 1, 2,
        0, 2, 3,
        4, 5, 6
    };

    SetTable(numV, vpair, numT, itriple);
}

void MarchingCubes::Bits5EdgeFaceDiag(int index[8])
{
    int const numV = 7;
    int vpair[2 * numV] =
    {
        index[0], index[1],
        index[4], index[5],
        index[4], index[6],
        index[0], index[2],
        index[2], index[3],
        index[3], index[7],
        index[1], index[3]
    };

    // Not the reverse ordering from Bit3EdgeFaceDiag due to ambiguous face
    // handling.
    int const numT = 5;
    int itriple[3 * numT] =
    {
        5, 0, 6,
        5, 1, 0,
        5, 2, 1,
        5, 3, 2,
        5, 4, 3
    };

    SetTable(numV, vpair, numT, itriple);
}

void MarchingCubes::Bits3FaceDiagFaceDiag(int index[8])
{
    int const numV = 9;
    int vpair[2 * numV] =
    {
        index[0], index[1],
        index[0], index[4],
        index[0], index[2],
        index[2], index[3],
        index[3], index[7],
        index[1], index[3],
        index[1], index[5],
        index[5], index[7],
        index[4], index[5]
    };

    int const numT = 3;
    int itriple[3 * numT] =
    {
        0, 1, 2,
        3, 4, 5,
        6, 7, 8
    };

    SetTable(numV, vpair, numT, itriple);
}

void MarchingCubes::Bits5FaceDiagFaceDiag(int index[8])
{
    int const numV = 9;
    int vpair[2 * numV] =
    {
        index[0], index[1],
        index[0], index[4],
        index[0], index[2],
        index[2], index[3],
        index[3], index[7],
        index[1], index[3],
        index[1], index[5],
        index[5], index[7],
        index[4], index[5]
    };

    // Not the reverse ordering from Bit3FaceDiagFaceDiag due to ambiguous face
    // handling.
    int const numT = 5;
    int itriple[3 * numT] =
    {
        1, 3, 2,
        1, 4, 3,
        1, 7, 4,
        1, 8, 7,
        0, 5, 6
    };

    SetTable(numV, vpair, numT, itriple);
}

void MarchingCubes::Bits4SameFace(int index[8])
{
    int const numV = 4;
    int vpair[2 * numV] =
    {
        index[0], index[4],
        index[2], index[6],
        index[3], index[7],
        index[1], index[5]
    };

    int const numT = 2;
    int itriple[3 * numT] =
    {
        0, 1, 2,
        0, 2, 3
    };

    SetTable(numV, vpair, numT, itriple);
}

void MarchingCubes::Bits4FaceEdge(int index[8])
{
    int const numV = 6;
    int vpair[2 * numV] =
    {
        index[4], index[5],
        index[4], index[6],
        index[2], index[6],
        index[2], index[3],
        index[1], index[3],
        index[1], index[5]
    };

    int const numT = 4;
    int itriple[3 * numT] =
    {
        0, 1, 2,
        0, 2, 3,
        0, 3, 4,
        0, 4, 5
    };

    SetTable(numV, vpair, numT, itriple);
}

void MarchingCubes::Bits4FaceFaceDiagL(int index[8])
{
    int const numV = 6;
    int vpair[2 * numV] =
    {
        index[4], index[5],
        index[0], index[4],
        index[2], index[6],
        index[2], index[3],
        index[1], index[3],
        index[5], index[7]
    };

    int const numT = 4;
    int itriple[3 * numT] =
    {
        0, 1, 2,
        0, 2, 3,
        0, 3, 4,
        0, 4, 5
    };

    SetTable(numV, vpair, numT, itriple);
}

void MarchingCubes::Bits4FaceFaceDiagR(int index[8])
{
    int const numV = 6;
    int vpair[2 * numV] =
    {
        index[4], index[6],
        index[6], index[7],
        index[2], index[3],
        index[1], index[3],
        index[1], index[5],
        index[0], index[4]
    };

    int const numT = 4;
    int itriple[3 * numT] =
    {
        0, 1, 2,
        0, 2, 3,
        0, 3, 4,
        0, 4, 5
    };

    SetTable(numV, vpair, numT, itriple);
}

void MarchingCubes::Bits4FaceBoxDiag(int index[8])
{
    int const numV = 8;
    int vpair[2 * numV] =
    {
        index[0], index[4],
        index[2], index[6],
        index[2], index[3],
        index[1], index[3],
        index[1], index[5],
        index[6], index[7],
        index[5], index[7],
        index[3], index[7]
    };

    int const numT = 4;
    int itriple[3 * numT] =
    {
        0, 1, 2,
        0, 2, 3,
        0, 3, 4,
        5, 6, 7
    };

    SetTable(numV, vpair, numT, itriple);
}

void MarchingCubes::Bits4EdgeEdgePara(int index[8])
{
    int const numV = 8;
    int vpair[2 * numV] =
    {
        index[0], index[4],
        index[0], index[2],
        index[1], index[3],
        index[1], index[5],
        index[2], index[6],
        index[4], index[6],
        index[5], index[7],
        index[3], index[7]
    };

    int const numT = 4;
    int itriple[3 * numT] =
    {
        0, 1, 2,
        0, 2, 3,
        4, 5, 6,
        4, 6, 7
    };

    SetTable(numV, vpair, numT, itriple);
}

void MarchingCubes::Bits4EdgeEdgePerp(int index[8])
{
    int const numV = 12;
    int vpair[2 * numV] =
    {
        index[0], index[1],
        index[0], index[4],
        index[0], index[2],
        index[2], index[6],
        index[4], index[6],
        index[6], index[7],
        index[2], index[3],
        index[3], index[7],
        index[1], index[3],
        index[1], index[5],
        index[5], index[7],
        index[4], index[5]
    };

    int const numT = 4;
    int itriple[3 * numT] =
    {
        0, 1, 2,
        3, 4, 5,
        6, 7, 8,
        9, 10, 11
    };

    SetTable(numV, vpair, numT, itriple);

}

MarchingCubes::Topology::Topology()
    :
    numVertices(0),
    numTriangles(0)
{
    std::fill(vpair.begin(), vpair.end(), std::array<int, 2>{ 0, 0 });
    std::fill(itriple.begin(), itriple.end(), std::array<int, 3>{ 0, 0, 0 });
}

#define MC_ENTRY(name) CT_##name, &MarchingCubes::name

MarchingCubes::Configuration MarchingCubes::msConfiguration[256] =
{
    /*00000000*/{ MC_ENTRY(Bits0), 0, 1, 2, 3, 4, 5, 6, 7 },
    /*00000001*/{ MC_ENTRY(Bits1), 0, 1, 2, 3, 4, 5, 6, 7 },
    /*00000010*/{ MC_ENTRY(Bits1), 1, 3, 0, 2, 5, 7, 4, 6 },
    /*00000011*/{ MC_ENTRY(Bits2Edge), 0, 1, 2, 3, 4, 5, 6, 7 },
    /*00000100*/{ MC_ENTRY(Bits1), 2, 0, 3, 1, 6, 4, 7, 5 },
    /*00000101*/{ MC_ENTRY(Bits2Edge), 2, 0, 3, 1, 6, 4, 7, 5 },
    /*00000110*/{ MC_ENTRY(Bits2FaceDiag), 1, 3, 0, 2, 5, 7, 4, 6 },
    /*00000111*/{ MC_ENTRY(Bits3SameFace), 0, 1, 2, 3, 4, 5, 6, 7 },
    /*00001000*/{ MC_ENTRY(Bits1), 3, 2, 1, 0, 7, 6, 5, 4 },
    /*00001001*/{ MC_ENTRY(Bits2FaceDiag), 0, 1, 2, 3, 4, 5, 6, 7 },
    /*00001010*/{ MC_ENTRY(Bits2Edge), 1, 3, 0, 2, 5, 7, 4, 6 },
    /*00001011*/{ MC_ENTRY(Bits3SameFace), 1, 3, 0, 2, 5, 7, 4, 6 },
    /*00001100*/{ MC_ENTRY(Bits2Edge), 3, 2, 1, 0, 7, 6, 5, 4 },
    /*00001101*/{ MC_ENTRY(Bits3SameFace), 2, 0, 3, 1, 6, 4, 7, 5 },
    /*00001110*/{ MC_ENTRY(Bits3SameFace), 3, 2, 1, 0, 7, 6, 5, 4 },
    /*00001111*/{ MC_ENTRY(Bits4SameFace), 0, 1, 2, 3, 4, 5, 6, 7 },
    /*00010000*/{ MC_ENTRY(Bits1), 4, 5, 0, 1, 6, 7, 2, 3 },
    /*00010001*/{ MC_ENTRY(Bits2Edge), 4, 0, 6, 2, 5, 1, 7, 3 },
    /*00010010*/{ MC_ENTRY(Bits2FaceDiag), 1, 0, 5, 4, 3, 2, 7, 6 },
    /*00010011*/{ MC_ENTRY(Bits3SameFace), 0, 4, 1, 5, 2, 6, 3, 7 },
    /*00010100*/{ MC_ENTRY(Bits2FaceDiag), 4, 0, 6, 2, 5, 1, 7, 3 },
    /*00010101*/{ MC_ENTRY(Bits3SameFace), 0, 2, 4, 6, 1, 3, 5, 7 },
    /*00010110*/{ MC_ENTRY(Bits3FaceDiagFaceDiag),  2, 0, 3, 1, 6, 4, 7, 5 },
    /*00010111*/{ MC_ENTRY(Bits4FaceEdge), 0, 1, 2, 3, 4, 5, 6, 7 },
    /*00011000*/{ MC_ENTRY(Bits2BoxDiag), 3, 2, 1, 0, 7, 6, 5, 4 },
    /*00011001*/{ MC_ENTRY(Bits3EdgeFaceDiag), 0, 1, 2, 3, 4, 5, 6, 7 },
    /*00011010*/{ MC_ENTRY(Bits3EdgeFaceDiag), 1, 0, 5, 4, 3, 2, 7, 6 },
    /*00011011*/{ MC_ENTRY(Bits4FaceFaceDiagR), 1, 3, 0, 2, 5, 7, 4, 6 },
    /*00011100*/{ MC_ENTRY(Bits3EdgeFaceDiag), 2, 6, 0, 4, 3, 7, 1, 5 },
    /*00011101*/{ MC_ENTRY(Bits4FaceFaceDiagL), 2, 0, 3, 1, 6, 4, 7, 5 },
    /*00011110*/{ MC_ENTRY(Bits4FaceBoxDiag), 3, 2, 1, 0, 7, 6, 5, 4 },
    /*00011111*/{ MC_ENTRY(Bits5SameFace), 7, 5, 6, 4, 3, 1, 2, 0 },
    /*00100000*/{ MC_ENTRY(Bits1), 5, 7, 1, 3, 4, 6, 0, 2 },
    /*00100001*/{ MC_ENTRY(Bits2FaceDiag), 0, 4, 1, 5, 2, 6, 3, 7 },
    /*00100010*/{ MC_ENTRY(Bits2Edge), 5, 1, 4, 0, 7, 3, 6, 2 },
    /*00100011*/{ MC_ENTRY(Bits3SameFace), 1, 0, 5, 4, 3, 2, 7, 6 },
    /*00100100*/{ MC_ENTRY(Bits2BoxDiag), 2, 0, 3, 1, 6, 4, 7, 5 },
    /*00100101*/{ MC_ENTRY(Bits3EdgeFaceDiag), 0, 4, 1, 5, 2, 6, 3, 7 },
    /*00100110*/{ MC_ENTRY(Bits3EdgeFaceDiag), 1, 3, 0, 2, 5, 7, 4, 6 },
    /*00100111*/{ MC_ENTRY(Bits4FaceFaceDiagL), 0, 1, 2, 3, 4, 5, 6, 7 },
    /*00101000*/{ MC_ENTRY(Bits2FaceDiag), 5, 7, 1, 3, 4, 6, 0, 2 },
    /*00101001*/{ MC_ENTRY(Bits3FaceDiagFaceDiag), 0, 1, 2, 3, 4, 5, 6, 7 },
    /*00101010*/{ MC_ENTRY(Bits3SameFace), 1, 5, 3, 7, 0, 4, 2, 6 },
    /*00101011*/{ MC_ENTRY(Bits4FaceEdge), 1, 3, 0, 2, 5, 7, 4, 6 },
    /*00101100*/{ MC_ENTRY(Bits3EdgeFaceDiag), 3, 1, 7, 5, 2, 0, 6, 4 },
    /*00101101*/{ MC_ENTRY(Bits4FaceBoxDiag), 2, 0, 3, 1, 6, 4, 7, 5 },
    /*00101110*/{ MC_ENTRY(Bits4FaceFaceDiagR), 3, 2, 1, 0, 7, 6, 5, 4 },
    /*00101111*/{ MC_ENTRY(Bits5SameFace), 6, 7, 4, 5, 2, 3, 0, 1 },
    /*00110000*/{ MC_ENTRY(Bits2Edge), 4, 5, 0, 1, 6, 7, 2, 3 },
    /*00110001*/{ MC_ENTRY(Bits3SameFace), 4, 5, 0, 1, 6, 7, 2, 3 },
    /*00110010*/{ MC_ENTRY(Bits3SameFace), 5, 1, 4, 0, 7, 3, 6, 2 },
    /*00110011*/{ MC_ENTRY(Bits4SameFace), 0, 4, 1, 5, 2, 6, 3, 7 },
    /*00110100*/{ MC_ENTRY(Bits3EdgeFaceDiag), 4, 0, 6, 2, 5, 1, 7, 3 },
    /*00110101*/{ MC_ENTRY(Bits4FaceFaceDiagR), 0, 2, 4, 6, 1, 3, 5, 7 },
    /*00110110*/{ MC_ENTRY(Bits4FaceBoxDiag), 5, 1, 4, 0, 7, 3, 6, 2 },
    /*00110111*/{ MC_ENTRY(Bits5SameFace), 7, 6, 3, 2, 5, 4, 1, 0 },
    /*00111000*/{ MC_ENTRY(Bits3EdgeFaceDiag), 5, 7, 1, 3, 4, 6, 0, 2 },
    /*00111001*/{ MC_ENTRY(Bits4FaceBoxDiag), 4, 5, 0, 1, 6, 7, 2, 3 },
    /*00111010*/{ MC_ENTRY(Bits4FaceFaceDiagL), 5, 1, 4, 0, 7, 3, 6, 2 },
    /*00111011*/{ MC_ENTRY(Bits5SameFace), 6, 2, 7, 3, 4, 0, 5, 1 },
    /*00111100*/{ MC_ENTRY(Bits4EdgeEdgePara), 3, 2, 1, 0, 7, 6, 5, 4 },
    /*00111101*/{ MC_ENTRY(Bits5EdgeFaceDiag), 7, 3, 5, 1, 6, 2, 4, 0 },
    /*00111110*/{ MC_ENTRY(Bits5EdgeFaceDiag), 6, 4, 2, 0, 7, 5, 3, 1 },
    /*00111111*/{ MC_ENTRY(Bits6Edge), 6, 7, 4, 5, 2, 3, 0, 1 },
    /*01000000*/{ MC_ENTRY(Bits1), 6, 7, 4, 5, 2, 3, 0, 1 },
    /*01000001*/{ MC_ENTRY(Bits2FaceDiag), 0, 2, 4, 6, 1, 3, 5, 7 },
    /*01000010*/{ MC_ENTRY(Bits2BoxDiag), 1, 3, 0, 2, 5, 7, 4, 6 },
    /*01000011*/{ MC_ENTRY(Bits3EdgeFaceDiag), 0, 2, 4, 6, 1, 3, 5, 7 },
    /*01000100*/{ MC_ENTRY(Bits2Edge), 6, 2, 7, 3, 4, 0, 5, 1 },
    /*01000101*/{ MC_ENTRY(Bits3SameFace), 2, 6, 0, 4, 3, 7, 1, 5 },
    /*01000110*/{ MC_ENTRY(Bits3EdgeFaceDiag), 2, 0, 3, 1, 6, 4, 7, 5 },
    /*01000111*/{ MC_ENTRY(Bits4FaceFaceDiagR), 0, 1, 2, 3, 4, 5, 6, 7 },
    /*01001000*/{ MC_ENTRY(Bits2FaceDiag), 3, 7, 2, 6, 1, 5, 0, 4 },
    /*01001001*/{ MC_ENTRY(Bits3FaceDiagFaceDiag), 3, 2, 1, 0, 7, 6, 5, 4 },
    /*01001010*/{ MC_ENTRY(Bits3EdgeFaceDiag), 3, 7, 2, 6, 1, 5, 0, 4 },
    /*01001011*/{ MC_ENTRY(Bits4FaceBoxDiag), 1, 3, 0, 2, 5, 7, 4, 6 },
    /*01001100*/{ MC_ENTRY(Bits3SameFace), 2, 3, 6, 7, 0, 1, 4, 5 },
    /*01001101*/{ MC_ENTRY(Bits4FaceEdge), 2, 0, 3, 1, 6, 4, 7, 5 },
    /*01001110*/{ MC_ENTRY(Bits4FaceFaceDiagL), 3, 2, 1, 0, 7, 6, 5, 4 },
    /*01001111*/{ MC_ENTRY(Bits5SameFace), 5, 4, 7, 6, 1, 0, 3, 2 },
    /*01010000*/{ MC_ENTRY(Bits2Edge), 6, 4, 2, 0, 7, 5, 3, 1 },
    /*01010001*/{ MC_ENTRY(Bits3SameFace), 4, 0, 6, 2, 5, 1, 7, 3 },
    /*01010010*/{ MC_ENTRY(Bits3EdgeFaceDiag), 4, 5, 0, 1, 6, 7, 2, 3 },
    /*01010011*/{ MC_ENTRY(Bits4FaceFaceDiagL), 0, 4, 1, 5, 2, 6, 3, 7 },
    /*01010100*/{ MC_ENTRY(Bits3SameFace), 6, 4, 2, 0, 7, 5, 3, 1 },
    /*01010101*/{ MC_ENTRY(Bits4SameFace), 0, 2, 4, 6, 1, 3, 5, 7 },
    /*01010110*/{ MC_ENTRY(Bits4FaceBoxDiag), 6, 4, 2, 0, 7, 5, 3, 1 },
    /*01010111*/{ MC_ENTRY(Bits5SameFace), 7, 3, 5, 1, 6, 2, 4, 0 },
    /*01011000*/{ MC_ENTRY(Bits3EdgeFaceDiag), 6, 2, 7, 3, 4, 0, 5, 1 },
    /*01011001*/{ MC_ENTRY(Bits4FaceBoxDiag), 4, 0, 6, 2, 5, 1, 7, 3 },
    /*01011010*/{ MC_ENTRY(Bits4EdgeEdgePara), 1, 3, 0, 2, 5, 7, 4, 6 },
    /*01011011*/{ MC_ENTRY(Bits5EdgeFaceDiag), 7, 6, 3, 2, 5, 4, 1, 0 },
    /*01011100*/{ MC_ENTRY(Bits4FaceFaceDiagR), 6, 4, 2, 0, 7, 5, 3, 1 },
    /*01011101*/{ MC_ENTRY(Bits5SameFace), 5, 7, 1, 3, 4, 6, 0, 2 },
    /*01011110*/{ MC_ENTRY(Bits5EdgeFaceDiag), 5, 1, 4, 0, 7, 3, 6, 2 },
    /*01011111*/{ MC_ENTRY(Bits6Edge), 5, 7, 1, 3, 4, 6, 0, 2 },
    /*01100000*/{ MC_ENTRY(Bits2FaceDiag), 5, 4, 7, 6, 1, 0, 3, 2 },
    /*01100001*/{ MC_ENTRY(Bits3FaceDiagFaceDiag), 5, 4, 7, 6, 1, 0, 3, 2 },
    /*01100010*/{ MC_ENTRY(Bits3EdgeFaceDiag), 5, 4, 7, 6, 1, 0, 3, 2 },
    /*01100011*/{ MC_ENTRY(Bits4FaceBoxDiag), 1, 0, 5, 4, 3, 2, 7, 6 },
    /*01100100*/{ MC_ENTRY(Bits3EdgeFaceDiag), 6, 7, 4, 5, 2, 3, 0, 1 },
    /*01100101*/{ MC_ENTRY(Bits4FaceBoxDiag), 2, 6, 0, 4, 3, 7, 1, 5 },
    /*01100110*/{ MC_ENTRY(Bits4EdgeEdgePara), 6, 2, 7, 3, 4, 0, 5, 1 },
    /*01100111*/{ MC_ENTRY(Bits5EdgeFaceDiag), 7, 5, 6, 4, 3, 1, 2, 0 },
    /*01101000*/{ MC_ENTRY(Bits3FaceDiagFaceDiag), 6, 7, 4, 5, 2, 3, 0, 1 },
    /*01101001*/{ MC_ENTRY(Bits4EdgeEdgePerp), 0, 1, 2, 3, 4, 5, 6, 7 },
    /*01101010*/{ MC_ENTRY(Bits4FaceBoxDiag), 1, 5, 3, 7, 0, 4, 2, 6 },
    /*01101011*/{ MC_ENTRY(Bits5FaceDiagFaceDiag), 4, 6, 5, 7, 0, 2, 1, 3 },
    /*01101100*/{ MC_ENTRY(Bits4FaceBoxDiag), 2, 3, 6, 7, 0, 1, 4, 5 },
    /*01101101*/{ MC_ENTRY(Bits5FaceDiagFaceDiag), 7, 5, 6, 4, 3, 1, 2, 0 },
    /*01101110*/{ MC_ENTRY(Bits5EdgeFaceDiag), 4, 6, 5, 7, 0, 2, 1, 3 },
    /*01101111*/{ MC_ENTRY(Bits6FaceDiag), 7, 5, 6, 4, 3, 1, 2, 0 },
    /*01110000*/{ MC_ENTRY(Bits3SameFace), 4, 6, 5, 7, 0, 2, 1, 3 },
    /*01110001*/{ MC_ENTRY(Bits4FaceEdge), 4, 6, 5, 7, 0, 2, 1, 3 },
    /*01110010*/{ MC_ENTRY(Bits4FaceFaceDiagR), 5, 1, 4, 0, 7, 3, 6, 2 },
    /*01110011*/{ MC_ENTRY(Bits5SameFace), 3, 7, 2, 6, 1, 5, 0, 4 },
    /*01110100*/{ MC_ENTRY(Bits4FaceFaceDiagL), 4, 6, 5, 7, 0, 2, 1, 3 },
    /*01110101*/{ MC_ENTRY(Bits5SameFace), 3, 1, 7, 5, 2, 0, 6, 4 },
    /*01110110*/{ MC_ENTRY(Bits5EdgeFaceDiag), 3, 2, 1, 0, 7, 6, 5, 4 },
    /*01110111*/{ MC_ENTRY(Bits6Edge), 3, 7, 2, 6, 1, 5, 0, 4 },
    /*01111000*/{ MC_ENTRY(Bits4FaceBoxDiag), 4, 6, 5, 7, 0, 2, 1, 3 },
    /*01111001*/{ MC_ENTRY(Bits5FaceDiagFaceDiag), 1, 3, 0, 2, 5, 7, 4, 6 },
    /*01111010*/{ MC_ENTRY(Bits5EdgeFaceDiag), 2, 3, 6, 7, 0, 1, 4, 5 },
    /*01111011*/{ MC_ENTRY(Bits6FaceDiag), 2, 3, 6, 7, 0, 1, 4, 5 },
    /*01111100*/{ MC_ENTRY(Bits5EdgeFaceDiag), 1, 5, 3, 7, 0, 4, 2, 6 },
    /*01111101*/{ MC_ENTRY(Bits6FaceDiag), 1, 5, 3, 7, 0, 4, 2, 6 },
    /*01111110*/{ MC_ENTRY(Bits6BoxDiag), 0, 1, 2, 3, 4, 5, 6, 7 },
    /*01111111*/{ MC_ENTRY(Bits7), 7, 3, 5, 1, 6, 2, 4, 0 },
    /*10000000*/{ MC_ENTRY(Bits1), 7, 3, 5, 1, 6, 2, 4, 0 },
    /*10000001*/{ MC_ENTRY(Bits2BoxDiag), 0, 1, 2, 3, 4, 5, 6, 7 },
    /*10000010*/{ MC_ENTRY(Bits2FaceDiag), 1, 5, 3, 7, 0, 4, 2, 6 },
    /*10000011*/{ MC_ENTRY(Bits3EdgeFaceDiag), 1, 5, 3, 7, 0, 4, 2, 6 },
    /*10000100*/{ MC_ENTRY(Bits2FaceDiag), 2, 3, 6, 7, 0, 1, 4, 5 },
    /*10000101*/{ MC_ENTRY(Bits3EdgeFaceDiag), 2, 3, 6, 7, 0, 1, 4, 5 },
    /*10000110*/{ MC_ENTRY(Bits3FaceDiagFaceDiag), 1, 3, 0, 2, 5, 7, 4, 6 },
    /*10000111*/{ MC_ENTRY(Bits4FaceBoxDiag), 0, 1, 2, 3, 4, 5, 6, 7 },
    /*10001000*/{ MC_ENTRY(Bits2Edge), 3, 7, 2, 6, 1, 5, 0, 4 },
    /*10001001*/{ MC_ENTRY(Bits3EdgeFaceDiag), 3, 2, 1, 0, 7, 6, 5, 4 },
    /*10001010*/{ MC_ENTRY(Bits3SameFace), 3, 1, 7, 5, 2, 0, 6, 4 },
    /*10001011*/{ MC_ENTRY(Bits4FaceFaceDiagL), 1, 3, 0, 2, 5, 7, 4, 6 },
    /*10001100*/{ MC_ENTRY(Bits3SameFace), 3, 7, 2, 6, 1, 5, 0, 4 },
    /*10001101*/{ MC_ENTRY(Bits4FaceFaceDiagR), 2, 0, 3, 1, 6, 4, 7, 5 },
    /*10001110*/{ MC_ENTRY(Bits4FaceEdge), 3, 2, 1, 0, 7, 6, 5, 4 },
    /*10001111*/{ MC_ENTRY(Bits5SameFace), 4, 6, 5, 7, 0, 2, 1, 3 },
    /*10010000*/{ MC_ENTRY(Bits2FaceDiag), 7, 5, 6, 4, 3, 1, 2, 0 },
    /*10010001*/{ MC_ENTRY(Bits3EdgeFaceDiag), 4, 6, 5, 7, 0, 2, 1, 3 },
    /*10010010*/{ MC_ENTRY(Bits3FaceDiagFaceDiag), 7, 5, 6, 4, 3, 1, 2, 0 },
    /*10010011*/{ MC_ENTRY(Bits4FaceBoxDiag), 0, 4, 1, 5, 2, 6, 3, 7 },
    /*10010100*/{ MC_ENTRY(Bits3FaceDiagFaceDiag), 4, 6, 5, 7, 0, 2, 1, 3 },
    /*10010101*/{ MC_ENTRY(Bits4FaceBoxDiag), 0, 2, 4, 6, 1, 3, 5, 7 },
    /*10010110*/{ MC_ENTRY(Bits4EdgeEdgePerp), 1, 3, 0, 2, 5, 7, 4, 6 },
    /*10010111*/{ MC_ENTRY(Bits5FaceDiagFaceDiag), 6, 7, 4, 5, 2, 3, 0, 1 },
    /*10011000*/{ MC_ENTRY(Bits3EdgeFaceDiag), 7, 5, 6, 4, 3, 1, 2, 0 },
    /*10011001*/{ MC_ENTRY(Bits4EdgeEdgePara), 4, 0, 6, 2, 5, 1, 7, 3 },
    /*10011010*/{ MC_ENTRY(Bits4FaceBoxDiag), 3, 1, 7, 5, 2, 0, 6, 4 },
    /*10011011*/{ MC_ENTRY(Bits5EdgeFaceDiag), 6, 7, 4, 5, 2, 3, 0, 1 },
    /*10011100*/{ MC_ENTRY(Bits4FaceBoxDiag), 3, 7, 2, 6, 1, 5, 0, 4 },
    /*10011101*/{ MC_ENTRY(Bits5EdgeFaceDiag), 5, 4, 7, 6, 1, 0, 3, 2 },
    /*10011110*/{ MC_ENTRY(Bits5FaceDiagFaceDiag), 5, 4, 7, 6, 1, 0, 3, 2 },
    /*10011111*/{ MC_ENTRY(Bits6FaceDiag), 5, 4, 7, 6, 1, 0, 3, 2 },
    /*10100000*/{ MC_ENTRY(Bits2Edge), 5, 7, 1, 3, 4, 6, 0, 2 },
    /*10100001*/{ MC_ENTRY(Bits3EdgeFaceDiag), 5, 1, 4, 0, 7, 3, 6, 2 },
    /*10100010*/{ MC_ENTRY(Bits3SameFace), 5, 7, 1, 3, 4, 6, 0, 2 },
    /*10100011*/{ MC_ENTRY(Bits4FaceFaceDiagR), 1, 0, 5, 4, 3, 2, 7, 6 },
    /*10100100*/{ MC_ENTRY(Bits3EdgeFaceDiag), 7, 6, 3, 2, 5, 4, 1, 0 },
    /*10100101*/{ MC_ENTRY(Bits4EdgeEdgePara), 2, 0, 3, 1, 6, 4, 7, 5 },
    /*10100110*/{ MC_ENTRY(Bits4FaceBoxDiag), 5, 7, 1, 3, 4, 6, 0, 2 },
    /*10100111*/{ MC_ENTRY(Bits5EdgeFaceDiag), 6, 2, 7, 3, 4, 0, 5, 1 },
    /*10101000*/{ MC_ENTRY(Bits3SameFace), 7, 3, 5, 1, 6, 2, 4, 0 },
    /*10101001*/{ MC_ENTRY(Bits4FaceBoxDiag), 7, 3, 5, 1, 6, 2, 4, 0 },
    /*10101010*/{ MC_ENTRY(Bits4SameFace), 1, 5, 3, 7, 0, 4, 2, 6 },
    /*10101011*/{ MC_ENTRY(Bits5SameFace), 6, 4, 2, 0, 7, 5, 3, 1 },
    /*10101100*/{ MC_ENTRY(Bits4FaceFaceDiagL), 3, 7, 2, 6, 1, 5, 0, 4 },
    /*10101101*/{ MC_ENTRY(Bits5EdgeFaceDiag), 4, 5, 0, 1, 6, 7, 2, 3 },
    /*10101110*/{ MC_ENTRY(Bits5SameFace), 4, 0, 6, 2, 5, 1, 7, 3 },
    /*10101111*/{ MC_ENTRY(Bits6Edge), 6, 4, 2, 0, 7, 5, 3, 1 },
    /*10110000*/{ MC_ENTRY(Bits3SameFace), 5, 4, 7, 6, 1, 0, 3, 2 },
    /*10110001*/{ MC_ENTRY(Bits4FaceFaceDiagL), 4, 5, 0, 1, 6, 7, 2, 3 },
    /*10110010*/{ MC_ENTRY(Bits4FaceEdge), 5, 1, 4, 0, 7, 3, 6, 2 },
    /*10110011*/{ MC_ENTRY(Bits5SameFace), 2, 3, 6, 7, 0, 1, 4, 5 },
    /*10110100*/{ MC_ENTRY(Bits4FaceBoxDiag), 5, 4, 7, 6, 1, 0, 3, 2 },
    /*10110101*/{ MC_ENTRY(Bits5EdgeFaceDiag), 3, 7, 2, 6, 1, 5, 0, 4 },
    /*10110110*/{ MC_ENTRY(Bits5FaceDiagFaceDiag), 3, 2, 1, 0, 7, 6, 5, 4 },
    /*10110111*/{ MC_ENTRY(Bits6FaceDiag), 3, 7, 2, 6, 1, 5, 0, 4 },
    /*10111000*/{ MC_ENTRY(Bits4FaceFaceDiagR), 7, 3, 5, 1, 6, 2, 4, 0 },
    /*10111001*/{ MC_ENTRY(Bits5EdgeFaceDiag), 2, 0, 3, 1, 6, 4, 7, 5 },
    /*10111010*/{ MC_ENTRY(Bits5SameFace), 2, 6, 0, 4, 3, 7, 1, 5 },
    /*10111011*/{ MC_ENTRY(Bits6Edge), 6, 2, 7, 3, 4, 0, 5, 1 },
    /*10111100*/{ MC_ENTRY(Bits5EdgeFaceDiag), 0, 2, 4, 6, 1, 3, 5, 7 },
    /*10111101*/{ MC_ENTRY(Bits6BoxDiag), 1, 3, 0, 2, 5, 7, 4, 6 },
    /*10111110*/{ MC_ENTRY(Bits6FaceDiag), 0, 2, 4, 6, 1, 3, 5, 7 },
    /*10111111*/{ MC_ENTRY(Bits7), 6, 7, 4, 5, 2, 3, 0, 1 },
    /*11000000*/{ MC_ENTRY(Bits2Edge), 6, 7, 4, 5, 2, 3, 0, 1 },
    /*11000001*/{ MC_ENTRY(Bits3EdgeFaceDiag), 6, 4, 2, 0, 7, 5, 3, 1 },
    /*11000010*/{ MC_ENTRY(Bits3EdgeFaceDiag), 7, 3, 5, 1, 6, 2, 4, 0 },
    /*11000011*/{ MC_ENTRY(Bits4EdgeEdgePara), 0, 1, 2, 3, 4, 5, 6, 7 },
    /*11000100*/{ MC_ENTRY(Bits3SameFace), 6, 2, 7, 3, 4, 0, 5, 1 },
    /*11000101*/{ MC_ENTRY(Bits4FaceFaceDiagL), 2, 6, 0, 4, 3, 7, 1, 5 },
    /*11000110*/{ MC_ENTRY(Bits4FaceBoxDiag), 6, 2, 7, 3, 4, 0, 5, 1 },
    /*11000111*/{ MC_ENTRY(Bits5EdgeFaceDiag), 5, 7, 1, 3, 4, 6, 0, 2 },
    /*11001000*/{ MC_ENTRY(Bits3SameFace), 7, 6, 3, 2, 5, 4, 1, 0 },
    /*11001001*/{ MC_ENTRY(Bits4FaceBoxDiag), 7, 6, 3, 2, 5, 4, 1, 0 },
    /*11001010*/{ MC_ENTRY(Bits4FaceFaceDiagR), 7, 6, 3, 2, 5, 4, 1, 0 },
    /*11001011*/{ MC_ENTRY(Bits5EdgeFaceDiag), 4, 0, 6, 2, 5, 1, 7, 3 },
    /*11001100*/{ MC_ENTRY(Bits4SameFace), 2, 3, 6, 7, 0, 1, 4, 5 },
    /*11001101*/{ MC_ENTRY(Bits5SameFace), 5, 1, 4, 0, 7, 3, 6, 2 },
    /*11001110*/{ MC_ENTRY(Bits5SameFace), 4, 5, 0, 1, 6, 7, 2, 3 },
    /*11001111*/{ MC_ENTRY(Bits6Edge), 4, 5, 0, 1, 6, 7, 2, 3 },
    /*11010000*/{ MC_ENTRY(Bits3SameFace), 6, 7, 4, 5, 2, 3, 0, 1 },
    /*11010001*/{ MC_ENTRY(Bits4FaceFaceDiagR), 4, 0, 6, 2, 5, 1, 7, 3 },
    /*11010010*/{ MC_ENTRY(Bits4FaceBoxDiag), 6, 7, 4, 5, 2, 3, 0, 1 },
    /*11010011*/{ MC_ENTRY(Bits5EdgeFaceDiag), 3, 1, 7, 5, 2, 0, 6, 4 },
    /*11010100*/{ MC_ENTRY(Bits4FaceEdge), 6, 4, 2, 0, 7, 5, 3, 1 },
    /*11010101*/{ MC_ENTRY(Bits5SameFace), 1, 5, 3, 7, 0, 4, 2, 6 },
    /*11010110*/{ MC_ENTRY(Bits5FaceDiagFaceDiag), 0, 1, 2, 3, 4, 5, 6, 7 },
    /*11010111*/{ MC_ENTRY(Bits6FaceDiag), 5, 7, 1, 3, 4, 6, 0, 2 },
    /*11011000*/{ MC_ENTRY(Bits4FaceFaceDiagL), 6, 7, 4, 5, 2, 3, 0, 1 },
    /*11011001*/{ MC_ENTRY(Bits5EdgeFaceDiag), 1, 3, 0, 2, 5, 7, 4, 6 },
    /*11011010*/{ MC_ENTRY(Bits5EdgeFaceDiag), 0, 4, 1, 5, 2, 6, 3, 7 },
    /*11011011*/{ MC_ENTRY(Bits6BoxDiag), 2, 0, 3, 1, 6, 4, 7, 5 },
    /*11011100*/{ MC_ENTRY(Bits5SameFace), 1, 0, 5, 4, 3, 2, 7, 6 },
    /*11011101*/{ MC_ENTRY(Bits6Edge), 5, 1, 4, 0, 7, 3, 6, 2 },
    /*11011110*/{ MC_ENTRY(Bits6FaceDiag), 0, 4, 1, 5, 2, 6, 3, 7 },
    /*11011111*/{ MC_ENTRY(Bits7), 5, 7, 1, 3, 4, 6, 0, 2 },
    /*11100000*/{ MC_ENTRY(Bits3SameFace), 7, 5, 6, 4, 3, 1, 2, 0 },
    /*11100001*/{ MC_ENTRY(Bits4FaceBoxDiag), 7, 5, 6, 4, 3, 1, 2, 0 },
    /*11100010*/{ MC_ENTRY(Bits4FaceFaceDiagL), 7, 5, 6, 4, 3, 1, 2, 0 },
    /*11100011*/{ MC_ENTRY(Bits5EdgeFaceDiag), 2, 6, 0, 4, 3, 7, 1, 5 },
    /*11100100*/{ MC_ENTRY(Bits4FaceFaceDiagR), 7, 5, 6, 4, 3, 1, 2, 0 },
    /*11100101*/{ MC_ENTRY(Bits5EdgeFaceDiag), 1, 0, 5, 4, 3, 2, 7, 6 },
    /*11100110*/{ MC_ENTRY(Bits5EdgeFaceDiag), 0, 1, 2, 3, 4, 5, 6, 7 },
    /*11100111*/{ MC_ENTRY(Bits6BoxDiag), 3, 2, 1, 0, 7, 6, 5, 4 },
    /*11101000*/{ MC_ENTRY(Bits4FaceEdge), 7, 5, 6, 4, 3, 1, 2, 0 },
    /*11101001*/{ MC_ENTRY(Bits5FaceDiagFaceDiag), 2, 0, 3, 1, 6, 4, 7, 5 },
    /*11101010*/{ MC_ENTRY(Bits5SameFace), 0, 2, 4, 6, 1, 3, 5, 7 },
    /*11101011*/{ MC_ENTRY(Bits6FaceDiag), 4, 0, 6, 2, 5, 1, 7, 3 },
    /*11101100*/{ MC_ENTRY(Bits5SameFace), 0, 4, 1, 5, 2, 6, 3, 7 },
    /*11101101*/{ MC_ENTRY(Bits6FaceDiag), 1, 0, 5, 4, 3, 2, 7, 6 },
    /*11101110*/{ MC_ENTRY(Bits6Edge), 4, 0, 6, 2, 5, 1, 7, 3 },
    /*11101111*/{ MC_ENTRY(Bits7), 4, 5, 0, 1, 6, 7, 2, 3 },
    /*11110000*/{ MC_ENTRY(Bits4SameFace), 4, 6, 5, 7, 0, 2, 1, 3 },
    /*11110001*/{ MC_ENTRY(Bits5SameFace), 3, 2, 1, 0, 7, 6, 5, 4 },
    /*11110010*/{ MC_ENTRY(Bits5SameFace), 2, 0, 3, 1, 6, 4, 7, 5 },
    /*11110011*/{ MC_ENTRY(Bits6Edge), 3, 2, 1, 0, 7, 6, 5, 4 },
    /*11110100*/{ MC_ENTRY(Bits5SameFace), 1, 3, 0, 2, 5, 7, 4, 6 },
    /*11110101*/{ MC_ENTRY(Bits6Edge), 1, 3, 0, 2, 5, 7, 4, 6 },
    /*11110110*/{ MC_ENTRY(Bits6FaceDiag), 0, 1, 2, 3, 4, 5, 6, 7 },
    /*11110111*/{ MC_ENTRY(Bits7), 3, 2, 1, 0, 7, 6, 5, 4 },
    /*11111000*/{ MC_ENTRY(Bits5SameFace), 0, 1, 2, 3, 4, 5, 6, 7 },
    /*11111001*/{ MC_ENTRY(Bits6FaceDiag), 1, 3, 0, 2, 5, 7, 4, 6 },
    /*11111010*/{ MC_ENTRY(Bits6Edge), 2, 0, 3, 1, 6, 4, 7, 5 },
    /*11111011*/{ MC_ENTRY(Bits7), 2, 0, 3, 1, 6, 4, 7, 5 },
    /*11111100*/{ MC_ENTRY(Bits6Edge), 0, 1, 2, 3, 4, 5, 6, 7 },
    /*11111101*/{ MC_ENTRY(Bits7), 1, 3, 0, 2, 5, 7, 4, 6 },
    /*11111110*/{ MC_ENTRY(Bits7), 0, 1, 2, 3, 4, 5, 6, 7 },
    /*11111111*/{ MC_ENTRY(Bits0), 0, 1, 2, 3, 4, 5, 6, 7 }
};

std::string MarchingCubes::msConfigurationString[CT_NUM_TYPES] =
{
    "Bits0",
    "Bits1",
    "Bits7",
    "Bits2Edge",
    "Bits6Edge",
    "Bits2FaceDiag",
    "Bits6FaceDiag",
    "Bits2BoxDiag",
    "Bits6BoxDiag",
    "Bits3SameFace",
    "Bits5SameFace",
    "Bits3EdgeFaceDiag",
    "Bits5EdgeFaceDiag",
    "Bits3FaceDiagFaceDiag",
    "Bits5FaceDiagFaceDiag",
    "Bits4SameFace",
    "Bits4FaceEdge",
    "Bits4FaceFaceDiagL",
    "Bits4FaceFaceDiagR",
    "Bits4FaceBoxDiag",
    "Bits4EdgeEdgePara",
    "Bits4EdgeEdgePerp"
};
