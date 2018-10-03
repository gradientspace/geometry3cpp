// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#include <GTEnginePCH.h>
#include <Mathematics/MSW/GteIntelSSE.h>
#include <Mathematics/GteConstants.h>
using namespace gte;

// Integer masks.
SIMD::Vector const SIMD::ZZZZ(0u);
SIMD::Vector const SIMD::ZZZF(0x00000000u, 0x00000000u, 0x00000000u, 0xFFFFFFFFu);
SIMD::Vector const SIMD::ZZFZ(0x00000000u, 0x00000000u, 0xFFFFFFFFu, 0x00000000u);
SIMD::Vector const SIMD::ZZFF(0x00000000u, 0x00000000u, 0xFFFFFFFFu, 0xFFFFFFFFu);
SIMD::Vector const SIMD::ZFZZ(0x00000000u, 0xFFFFFFFFu, 0x00000000u, 0x00000000u);
SIMD::Vector const SIMD::ZFZF(0x00000000u, 0xFFFFFFFFu, 0x00000000u, 0xFFFFFFFFu);
SIMD::Vector const SIMD::ZFFZ(0x00000000u, 0xFFFFFFFFu, 0xFFFFFFFFu, 0x00000000u);
SIMD::Vector const SIMD::ZFFF(0x00000000u, 0xFFFFFFFFu, 0xFFFFFFFFu, 0xFFFFFFFFu);
SIMD::Vector const SIMD::FZZZ(0xFFFFFFFFu, 0x00000000u, 0x00000000u, 0x00000000u);
SIMD::Vector const SIMD::FZZF(0xFFFFFFFFu, 0x00000000u, 0x00000000u, 0xFFFFFFFFu);
SIMD::Vector const SIMD::FZFZ(0xFFFFFFFFu, 0x00000000u, 0xFFFFFFFFu, 0x00000000u);
SIMD::Vector const SIMD::FZFF(0xFFFFFFFFu, 0x00000000u, 0xFFFFFFFFu, 0xFFFFFFFFu);
SIMD::Vector const SIMD::FFZZ(0xFFFFFFFFu, 0xFFFFFFFFu, 0x00000000u, 0x00000000u);
SIMD::Vector const SIMD::FFZF(0xFFFFFFFFu, 0xFFFFFFFFu, 0x00000000u, 0xFFFFFFFFu);
SIMD::Vector const SIMD::FFFZ(0xFFFFFFFFu, 0xFFFFFFFFu, 0xFFFFFFFFu, 0x00000000u);
SIMD::Vector const SIMD::FFFF(0xFFFFFFFFu);
SIMD::Vector const SIMD::SIGN(0x80000000u);
SIMD::Vector const SIMD::NSIGN(0x7FFFFFFFu);
SIMD::Vector const SIMD::NOFRC(0x00800000u);

// Numeric constants.
SIMD::Vector const SIMD::PZZZ(+1.0f,  0.0f,  0.0f,  0.0f);
SIMD::Vector const SIMD::ZPZZ( 0.0f, +1.0f,  0.0f,  0.0f);
SIMD::Vector const SIMD::ZZPZ( 0.0f,  0.0f, +1.0f,  0.0f);
SIMD::Vector const SIMD::ZZZP( 0.0f,  0.0f,  0.0f, +1.0f);
SIMD::Vector const SIMD::MZZZ(-1.0f,  0.0f,  0.0f,  0.0f);
SIMD::Vector const SIMD::ZMZZ( 0.0f, -1.0f,  0.0f,  0.0f);
SIMD::Vector const SIMD::ZZMZ( 0.0f,  0.0f, -1.0f,  0.0f);
SIMD::Vector const SIMD::ZZZM( 0.0f,  0.0f,  0.0f, -1.0f);
SIMD::Vector const SIMD::MMMM(-1.0f, -1.0f, -1.0f, -1.0f);
SIMD::Vector const SIMD::MMMP(-1.0f, -1.0f, -1.0f, +1.0f);
SIMD::Vector const SIMD::MMPM(-1.0f, -1.0f, +1.0f, -1.0f);
SIMD::Vector const SIMD::MMPP(-1.0f, -1.0f, +1.0f, +1.0f);
SIMD::Vector const SIMD::MPMM(-1.0f, +1.0f, -1.0f, -1.0f);
SIMD::Vector const SIMD::MPMP(-1.0f, +1.0f, -1.0f, +1.0f);
SIMD::Vector const SIMD::MPPM(-1.0f, +1.0f, +1.0f, -1.0f);
SIMD::Vector const SIMD::MPPP(-1.0f, +1.0f, +1.0f, +1.0f);
SIMD::Vector const SIMD::PMMM(+1.0f, -1.0f, -1.0f, -1.0f);
SIMD::Vector const SIMD::PMMP(+1.0f, -1.0f, -1.0f, +1.0f);
SIMD::Vector const SIMD::PMPM(+1.0f, -1.0f, +1.0f, -1.0f);
SIMD::Vector const SIMD::PMPP(+1.0f, -1.0f, +1.0f, +1.0f);
SIMD::Vector const SIMD::PPMM(+1.0f, +1.0f, -1.0f, -1.0f);
SIMD::Vector const SIMD::PPMP(+1.0f, +1.0f, -1.0f, +1.0f);
SIMD::Vector const SIMD::PPPM(+1.0f, +1.0f, +1.0f, -1.0f);
SIMD::Vector const SIMD::PPPP(+1.0f, +1.0f, +1.0f, +1.0f);
SIMD::Vector const SIMD::UNIT[4] = { PZZZ, ZPZZ, ZZPZ, ZZZP };

// Constants involving pi.
SIMD::Vector const SIMD::PI((float)GTE_C_PI);
SIMD::Vector const SIMD::HALF_PI((float)GTE_C_HALF_PI);
SIMD::Vector const SIMD::TWO_PI((float)GTE_C_TWO_PI);
SIMD::Vector const SIMD::INV_PI((float)GTE_C_INV_PI);
SIMD::Vector const SIMD::INV_TWO_PI((float)GTE_C_INV_TWO_PI);

// Constants to support approximations of sin(x).
SIMD::Vector const SIMD::C_SIN_APPR_DEG11_0((float)GTE_C_SIN_DEG11_C0);
SIMD::Vector const SIMD::C_SIN_APPR_DEG11_1((float)GTE_C_SIN_DEG11_C1);
SIMD::Vector const SIMD::C_SIN_APPR_DEG11_2((float)GTE_C_SIN_DEG11_C2);
SIMD::Vector const SIMD::C_SIN_APPR_DEG11_3((float)GTE_C_SIN_DEG11_C3);
SIMD::Vector const SIMD::C_SIN_APPR_DEG11_4((float)GTE_C_SIN_DEG11_C4);
SIMD::Vector const SIMD::C_SIN_APPR_DEG11_5((float)GTE_C_SIN_DEG11_C5);
SIMD::Vector const SIMD::C_SIN_APPR_DEG7_0((float)GTE_C_SIN_DEG7_C0);
SIMD::Vector const SIMD::C_SIN_APPR_DEG7_1((float)GTE_C_SIN_DEG7_C1);
SIMD::Vector const SIMD::C_SIN_APPR_DEG7_2((float)GTE_C_SIN_DEG7_C2);
SIMD::Vector const SIMD::C_SIN_APPR_DEG7_3((float)GTE_C_SIN_DEG7_C3);

// Constants to support approximations of cos(x).
SIMD::Vector const SIMD::C_COS_APPR_DEG10_0((float)GTE_C_COS_DEG10_C0);
SIMD::Vector const SIMD::C_COS_APPR_DEG10_1((float)GTE_C_COS_DEG10_C1);
SIMD::Vector const SIMD::C_COS_APPR_DEG10_2((float)GTE_C_COS_DEG10_C2);
SIMD::Vector const SIMD::C_COS_APPR_DEG10_3((float)GTE_C_COS_DEG10_C3);
SIMD::Vector const SIMD::C_COS_APPR_DEG10_4((float)GTE_C_COS_DEG10_C4);
SIMD::Vector const SIMD::C_COS_APPR_DEG10_5((float)GTE_C_COS_DEG10_C5);
SIMD::Vector const SIMD::C_COS_APPR_DEG6_0((float)GTE_C_COS_DEG6_C0);
SIMD::Vector const SIMD::C_COS_APPR_DEG6_1((float)GTE_C_COS_DEG6_C1);
SIMD::Vector const SIMD::C_COS_APPR_DEG6_2((float)GTE_C_COS_DEG6_C2);
SIMD::Vector const SIMD::C_COS_APPR_DEG6_3((float)GTE_C_COS_DEG6_C3);
