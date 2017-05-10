// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <Mathematics/GteBSNumber.h>
#include <Mathematics/GteBSRational.h>
#include <Mathematics/GteIEEEBinary16.h>
#include <Mathematics/GteUIntegerAP32.h>
#include <Mathematics/GteUIntegerFP32.h>
#include <cmath>
#include <type_traits>

namespace gte
{

class Arithmetic
{
public:
    // The tag-dispatch pattern is used for template-parameter-controlled
    // instantiation of mathematical functions.  See GteFunctions.h for
    // examples of how to use this type-trait system.

    enum Type
    {
        IS_INVALID,             // not an arithmetic type
        IS_FLOATING_POINT,      // 'float' or 'double' or 'long double'
        IS_FP16,                // IEEEBinary16
        IS_BINARY_SCIENTIFIC    // BSNumber or BSRational
    };

    typedef std::integral_constant<Type, IS_INVALID> IsInvalidType;
    typedef std::integral_constant<Type, IS_FLOATING_POINT> IsFPType;
    typedef std::integral_constant<Type, IS_FP16> IsFP16Type;
    typedef std::integral_constant<Type, IS_BINARY_SCIENTIFIC> IsBSType;

    template <typename T> struct WhichType : IsInvalidType{};
    template <typename U> struct WhichType<BSNumber<U>> : IsBSType{};
    template <typename U> struct WhichType<BSRational<U>> : IsBSType{};
};

template<> struct Arithmetic::WhichType<float> : IsFPType{};
template<> struct Arithmetic::WhichType<double> : IsFPType{};
template<> struct Arithmetic::WhichType<long double> : IsFPType{};
template<> struct Arithmetic::WhichType<IEEEBinary16> : IsFP16Type{};

}
