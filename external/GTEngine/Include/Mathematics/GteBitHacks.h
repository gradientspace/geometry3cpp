// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <GTEngineDEF.h>
#include <cstdint>

// Convenient macros for manipulating 64-bit integers.
#define GTE_I64(v) v##LL
#define GTE_U64(v) v##ULL
#define GTE_GET_LO_I64(v) (int32_t)((v) & 0x00000000ffffffffLL)
#define GTE_GET_HI_I64(v) (int32_t)(((v) >> 32) & 0x00000000ffffffffLL)
#define GTE_GET_LO_U64(v) (uint32_t)((v) & 0x00000000ffffffffULL)
#define GTE_GET_HI_U64(v) (uint32_t)(((v) >> 32) & 0x00000000ffffffffULL)
#define GTE_SET_LO_I64(v,lo) (((v) & 0xffffffff00000000LL) | (int64_t)(lo))
#define GTE_SET_HI_I64(v,hi) (((v) & 0x00000000ffffffffLL) | ((int64_t)(hi) << 32))
#define GTE_MAKE_I64(hi,lo)  ((int64_t)(lo) | ((int64_t)(hi) << 32))
#define GTE_SET_LO_U64(v,lo) (((v) & 0xffffffff00000000ULL) | (uint64_t)(lo))
#define GTE_SET_HI_U64(v,hi) (((v) & 0x00000000ffffffffULL) | ((uint64_t)(hi) << 32))
#define GTE_MAKE_U64(hi,lo)  ((uint64_t)(lo) | ((uint64_t)(hi) << 32))

namespace gte
{

GTE_IMPEXP bool IsPowerOfTwo(uint32_t value);
GTE_IMPEXP bool IsPowerOfTwo(int32_t value);

GTE_IMPEXP uint32_t Log2OfPowerOfTwo(uint32_t powerOfTwo);
GTE_IMPEXP int32_t Log2OfPowerOfTwo(int32_t powerOfTwo);

// Call these only for nonzero values.  If value is zero, then GetLeadingBit
// and GetTrailingBit return zero.
GTE_IMPEXP int32_t GetLeadingBit(uint32_t value);
GTE_IMPEXP int32_t GetLeadingBit(int32_t value);
GTE_IMPEXP int32_t GetLeadingBit(uint64_t value);
GTE_IMPEXP int32_t GetLeadingBit(int64_t value);
GTE_IMPEXP int32_t GetTrailingBit(uint32_t value);
GTE_IMPEXP int32_t GetTrailingBit(int32_t value);
GTE_IMPEXP int32_t GetTrailingBit(uint64_t value);
GTE_IMPEXP int32_t GetTrailingBit(int64_t value);

// Round up to a power of two.  If input is zero, the return is 1.  If input
// is larger than 2^{31}, the return is 2^{32}.
GTE_IMPEXP uint64_t RoundUpToPowerOfTwo(uint32_t value);

// Round down to a power of two.  If input is zero, the return is 0.
GTE_IMPEXP uint32_t RoundDownToPowerOfTwo(uint32_t value);

}
