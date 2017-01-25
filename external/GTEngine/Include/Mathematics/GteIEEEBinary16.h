// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <Mathematics/GteIEEEBinary.h>

namespace gte
{

struct _Float16 { int16_t _dummy; };

class IEEEBinary16 : public IEEEBinary<_Float16, uint16_t, 16, 11>
{
public:
    // Construction and destruction.  The base class destructor is hidden, but
    // this is safe because there are no side effects of the destruction.
    ~IEEEBinary16();
    IEEEBinary16();  // uninitialized
    IEEEBinary16(IEEEBinary16 const& object);
    IEEEBinary16(float number);
    IEEEBinary16(double number);
    IEEEBinary16(uint16_t encoding);

    // Implicit conversions.
    operator float() const;
    operator double() const;

    // Assignment.
    IEEEBinary16& operator=(IEEEBinary16 const& object);

    // Comparison.
    bool operator==(IEEEBinary16 const& object) const;
    bool operator!=(IEEEBinary16 const& object) const;
    bool operator< (IEEEBinary16 const& object) const;
    bool operator<=(IEEEBinary16 const& object) const;
    bool operator> (IEEEBinary16 const& object) const;
    bool operator>=(IEEEBinary16 const& object) const;

private:
    // Support for conversions between encodings.
    enum
    {
        F32_NUM_ENCODING_BITS = 32,
        F32_NUM_TRAILING_BITS = 23,
        F32_EXPONENT_BIAS = 127,
        F32_MAX_BIASED_EXPONENT = 255,
        F32_SIGN_MASK = 0x80000000,
        F32_NOT_SIGN_MASK = 0x7FFFFFFF,
        F32_BIASED_EXPONENT_MASK = 0x7F800000,
        F32_TRAILING_MASK = 0x007FFFFF,
        F16_AVR_MIN_SUBNORMAL_ZERO = 0x33000000,
        F16_MIN_SUBNORMAL = 0x33800000,
        F16_MIN_NORMAL = 0x38800000,
        F16_MAX_NORMAL = 0x477FE000,
        F16_AVR_MAX_NORMAL_INFINITY = 0x477FF000,
        DIFF_NUM_ENCODING_BITS = 16,
        DIFF_NUM_TRAILING_BITS = 13,
        DIFF_PAYLOAD_SHIFT = 13,
        INT_PART_MASK = 0x007FE000,
        FRC_PART_MASK = 0x00001FFF,
        FRC_HALF = 0x00001000
    };

    static uint16_t Convert32To16(uint32_t encoding);
    static uint32_t Convert16To32(uint16_t encoding);
};

// Arithmetic operations (high-precision).
IEEEBinary16 operator-(IEEEBinary16 x);
float operator+(IEEEBinary16 x, IEEEBinary16 y);
float operator-(IEEEBinary16 x, IEEEBinary16 y);
float operator*(IEEEBinary16 x, IEEEBinary16 y);
float operator/(IEEEBinary16 x, IEEEBinary16 y);
float operator+(IEEEBinary16 x, float y);
float operator-(IEEEBinary16 x, float y);
float operator*(IEEEBinary16 x, float y);
float operator/(IEEEBinary16 x, float y);
float operator+(float x, IEEEBinary16 y);
float operator-(float x, IEEEBinary16 y);
float operator*(float x, IEEEBinary16 y);
float operator/(float x, IEEEBinary16 y);

// Arithmetic updates.
IEEEBinary16& operator+=(IEEEBinary16& x, IEEEBinary16 y);
IEEEBinary16& operator-=(IEEEBinary16& x, IEEEBinary16 y);
IEEEBinary16& operator*=(IEEEBinary16& x, IEEEBinary16 y);
IEEEBinary16& operator/=(IEEEBinary16& x, IEEEBinary16 y);
IEEEBinary16& operator+=(IEEEBinary16& x, float y);
IEEEBinary16& operator-=(IEEEBinary16& x, float y);
IEEEBinary16& operator*=(IEEEBinary16& x, float y);
IEEEBinary16& operator/=(IEEEBinary16& x, float y);

}
