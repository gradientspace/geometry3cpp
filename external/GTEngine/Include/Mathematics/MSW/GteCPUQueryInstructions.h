// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.13.0 (2018/04/02)

#pragma once

#include <array>
#include <bitset>
#include <string>
#include <vector>

namespace gte
{

// Determine the capabilities of Intel or AMD CPU processors.  This class is
// based on sample code found at
// https://msdn.microsoft.com/en-us/library/hskdteyh.aspx

class CPUQueryInstructions
{
public:
    CPUQueryInstructions();

    inline std::string GetVendor() const { return mVendor; }
    inline std::string GetBrand() const { return mBrand; }

    inline bool IsSSE3() const { return mF1_ECX[0]; }
    inline bool IsPCLMULQDQ() const { return mF1_ECX[1]; }
    inline bool IsMONITOR() const { return mF1_ECX[3]; }
    inline bool IsSSSE3() const { return mF1_ECX[9]; }
    inline bool IsFMA() const { return mF1_ECX[12]; }
    inline bool IsCMPXCHG16B() const { return mF1_ECX[13]; }
    inline bool IsSSE41() const { return mF1_ECX[19]; }
    inline bool IsSSE42() const { return mF1_ECX[20]; }
    inline bool IsMOVBE() const { return mF1_ECX[22]; }
    inline bool IsPOPCNT() const { return mF1_ECX[23]; }
    inline bool IsAES() const { return mF1_ECX[25]; }
    inline bool IsXSAVE() const { return mF1_ECX[26]; }
    inline bool IsOSXSAVE() const { return mF1_ECX[27]; }
    inline bool IsAVX() const { return mF1_ECX[28]; }
    inline bool IsF16C() const { return mF1_ECX[29]; }
    inline bool IsRDRAND() const { return mF1_ECX[30]; }

    inline bool IsMSR() const { return mF1_EDX[5]; }
    inline bool IsCX8() const { return mF1_EDX[8]; }
    inline bool IsSEP() const { return mF1_EDX[11]; }
    inline bool IsCMOV() const { return mF1_EDX[15]; }
    inline bool IsCLFSH() const { return mF1_EDX[19]; }
    inline bool IsMMX() const { return mF1_EDX[23]; }
    inline bool IsFXSR() const { return mF1_EDX[24]; }
    inline bool IsSSE() const { return mF1_EDX[25]; }
    inline bool IsSSE2() const { return mF1_EDX[26]; }

    inline bool IsFSGSBASE() const { return mF7_EBX[0]; }
    inline bool IsBMI1() const { return mF7_EBX[3]; }
    inline bool IsHLE() const { return mIsIntel && mF7_EBX[4]; }
    inline bool IsAVX2() const { return mF7_EBX[5]; }
    inline bool IsBMI2() const { return mF7_EBX[8]; }
    inline bool IsERMS() const { return mF7_EBX[9]; }
    inline bool IsINVPCID() const { return mF7_EBX[10]; }
    inline bool IsRTM() const { return mIsIntel && mF7_EBX[11]; }
    inline bool IsAVX512F() const { return mF7_EBX[16]; }
    inline bool IsRDSEED() const { return mF7_EBX[18]; }
    inline bool IsADX() const { return mF7_EBX[19]; }
    inline bool IsAVX512PF() const { return mF7_EBX[26]; }
    inline bool IsAVX512ER() const { return mF7_EBX[27]; }
    inline bool IsAVX512CD() const { return mF7_EBX[28]; }
    inline bool IsSHA() const { return mF7_EBX[29]; }

    inline bool IsPREFETCHWT1() const { return mF7_ECX[0]; }

    inline bool IsLAHF() const { return mF81_ECX[0]; }
    inline bool IsLZCNT() const { return mIsIntel && mF81_ECX[5]; }
    inline bool IsABM() const { return mIsAMD && mF81_ECX[5]; }
    inline bool IsSSE4a() const { return mIsAMD && mF81_ECX[6]; }
    inline bool IsXOP() const { return mIsAMD && mF81_ECX[11]; }
    inline bool IsTBM() const { return mIsAMD && mF81_ECX[21]; }

    inline bool IsSYSCALL() const { return mIsIntel && mF81_EDX[11]; }
    inline bool IsMMXEXT() const { return mIsAMD && mF81_EDX[22]; }
    inline bool IsRDTSCP() const { return mIsIntel && mF81_EDX[27]; }
    inline bool Is3DNOWEXT() const { return mIsAMD && mF81_EDX[30]; }
    inline bool Is3DNOW() const { return mIsAMD && mF81_EDX[31]; }

private:
    int mNumIds;
    int mNumExIds;
    std::string mVendor;
    std::string mBrand;
    bool mIsIntel;
    bool mIsAMD;
    std::bitset<32> mF1_ECX;
    std::bitset<32> mF1_EDX;
    std::bitset<32> mF7_EBX;
    std::bitset<32> mF7_ECX;
    std::bitset<32> mF81_ECX;
    std::bitset<32> mF81_EDX;
    std::vector<std::array<int, 4>> mData;
    std::vector<std::array<int, 4>> mExtData;
};

}
