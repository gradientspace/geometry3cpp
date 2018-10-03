// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.13.0 (2018/04/02)

#include <GTEnginePCH.h>
#include <Mathematics/MSW/GteCPUQueryInstructions.h>
#include <intrin.h>
using namespace gte;

CPUQueryInstructions::CPUQueryInstructions()
    :
    mNumIds(0),
    mNumExIds(0),
    mIsIntel(false),
    mIsAMD(false),
    mF1_ECX{ 0 },
    mF1_EDX{ 0 },
    mF7_EBX{ 0 },
    mF7_ECX{ 0 },
    mF81_ECX{ 0 },
    mF81_EDX{ 0 },
    mData{},
    mExtData{}
{
    std::array<int, 4> cpui;

    // Calling __cpuid with 0x0 as the function_id argument gets the number of
    // the highest valid function ID.
    __cpuid(cpui.data(), 0);
    mNumIds = cpui[0];

    for (int i = 0; i <= mNumIds; ++i)
    {
        __cpuidex(cpui.data(), i, 0);
        mData.push_back(cpui);
    }

    // Get the vendor string.
    char vendor[0x20];
    memset(vendor, 0, sizeof(vendor));
    *reinterpret_cast<int*>(vendor) = mData[0][1];
    *reinterpret_cast<int*>(vendor + 4) = mData[0][3];
    *reinterpret_cast<int*>(vendor + 8) = mData[0][2];
    mVendor = vendor;
    if (mVendor == "GenuineIntel")
    {
        mIsIntel = true;
    }
    else if (mVendor == "AuthenticAMD")
    {
        mIsAMD = true;
    }

    // Get the bitset with flags for function 0x00000001.
    if (mNumIds >= 1)
    {
        mF1_ECX = mData[1][2];
        mF1_EDX = mData[1][3];
    }

    // Get the bitset with flags for function 0x00000007.
    if (mNumIds >= 7)
    {
        mF7_EBX = mData[7][1];
        mF7_ECX = mData[7][2];
    }

    // Calling __cpuid with 0x80000000 as the function_id argument gets the
    // number of the highest valid extended ID.
    __cpuid(cpui.data(), 0x80000000);
    mNumExIds = cpui[0];

    char brand[0x40];
    memset(brand, 0, sizeof(brand));

    for (int i = 0x80000000; i <= mNumExIds; ++i)
    {
        __cpuidex(cpui.data(), i, 0);
        mExtData.push_back(cpui);
    }

    // Get the bitset with flags for function 0x80000001.
    if (mNumExIds >= 0x80000001)
    {
        mF81_ECX = mExtData[1][2];
        mF81_EDX = mExtData[1][3];
    }

    // Interpret the CPU brand string if reported.
    if (mNumExIds >= 0x80000004)
    {
        memcpy(brand, mExtData[2].data(), sizeof(cpui));
        memcpy(brand + 16, mExtData[3].data(), sizeof(cpui));
        memcpy(brand + 32, mExtData[4].data(), sizeof(cpui));
        mBrand = brand;
    }
}
