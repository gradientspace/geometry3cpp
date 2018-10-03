// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.9.0 (2017/06/18)

#pragma once

#include <Mathematics/GteFunctions.h>
#include <Mathematics/GteBSNumber.h>

namespace gte
{
    template <typename UInteger>
    class Function<BSNumber<UInteger>>
    {
    public:
        inline static BSNumber<UInteger> ACos(BSNumber<UInteger> const& x)
        {
            return (BSNumber<UInteger>)acos((double)x);
        }

        inline static BSNumber<UInteger> ACosh(BSNumber<UInteger> x)
        {
            double y = (double)x;
            return (BSNumber<UInteger>)log(y + sqrt(y * y - 1.0));
        }

        inline static BSNumber<UInteger> ASin(BSNumber<UInteger> x)
        {
            return (BSNumber<UInteger>)asin((double)x);
        }

        inline static BSNumber<UInteger> ASinh(BSNumber<UInteger> x)
        {
            double y = (double)x;
            return (BSNumber<UInteger>)log(y + sqrt(y * y + 1.0));
        }

        inline static BSNumber<UInteger> ATan(BSNumber<UInteger> x)
        {
            return (BSNumber<UInteger>)atan((double)x);
        }

        inline static BSNumber<UInteger> ATanh(BSNumber<UInteger> x)
        {
            double y = (double)x;
            return (BSNumber<UInteger>)(log((1.0 + y) / (1.0 - y)) * 0.5);
        }

        inline static BSNumber<UInteger> ATan2(BSNumber<UInteger> y, BSNumber<UInteger> x)
        {
            return (BSNumber<UInteger>)atan2((double)y, (double)x);

        }

        inline static BSNumber<UInteger> ATanpi(BSNumber<UInteger> x)
        {
            return (BSNumber<UInteger>)(atan((double)x) * GTE_C_INV_PI);
        }

        inline static BSNumber<UInteger> ATan2pi(BSNumber<UInteger> y, BSNumber<UInteger> x)
        {
            return (BSNumber<UInteger>)(atan2((double)y, (double)x) * GTE_C_INV_PI);
        }

        inline static BSNumber<UInteger> Ceil(BSNumber<UInteger> x)
        {
            return (BSNumber<UInteger>)ceil((double)x);
        }

        inline static BSNumber<UInteger> Cos(BSNumber<UInteger> x)
        {
            return (BSNumber<UInteger>)cos((double)x);
        }

        inline static BSNumber<UInteger> Cosh(BSNumber<UInteger> x)
        {
            return (BSNumber<UInteger>)cosh((double)x);
        }

        inline static BSNumber<UInteger> Cospi(BSNumber<UInteger> x)
        {
            return (BSNumber<UInteger>)cos((double)(x * (BSNumber<UInteger>)GTE_C_PI));
        }

        inline static BSNumber<UInteger> Exp(BSNumber<UInteger> x)
        {
            return (BSNumber<UInteger>)exp((double)x);
        }

        inline static BSNumber<UInteger> Exp2(BSNumber<UInteger> x)
        {
            return (BSNumber<UInteger>)exp((double)(x * (BSNumber<UInteger>)GTE_C_LN_2));
        }

        inline static BSNumber<UInteger> Exp10(BSNumber<UInteger> x)
        {
            return (BSNumber<UInteger>)exp((double)(x * (BSNumber<UInteger>)GTE_C_LN_10));
        }

        inline static BSNumber<UInteger> FAbs(BSNumber<UInteger> x)
        {
            return (x.GetSign() >= 0 ? x : -x);
        }

        inline static BSNumber<UInteger> Floor(BSNumber<UInteger> x)
        {
            return (BSNumber<UInteger>)floor((double)x);
        }

        inline static BSNumber<UInteger> FMod(BSNumber<UInteger> x, BSNumber<UInteger> y)
        {
            return (BSNumber<UInteger>)fmod((double)x, (double)y);
        }

        inline static BSNumber<UInteger> FRExp(BSNumber<UInteger> x, int* exponent)
        {
            return (BSNumber<UInteger>)frexp((double)x, exponent);
        }

        inline static BSNumber<UInteger> InvSqrt(BSNumber<UInteger> x)
        {
            return (BSNumber<UInteger>)(1.0 / sqrt((double)x));
        }

        inline static BSNumber<UInteger> LDExp(BSNumber<UInteger> x, int exponent)
        {
            return (BSNumber<UInteger>)ldexp((double)x, exponent);
        }

        inline static BSNumber<UInteger> Log(BSNumber<UInteger> x)
        {
            return (BSNumber<UInteger>)log((double)x);
        }

        inline static BSNumber<UInteger> Log2(BSNumber<UInteger> x)
        {
            return (BSNumber<UInteger>)(log((double)x) * GTE_C_INV_LN_2);
        }

        inline static BSNumber<UInteger> Log10(BSNumber<UInteger> x)
        {
            return (BSNumber<UInteger>)log10((double)x);
        }

        inline static BSNumber<UInteger> Pow(BSNumber<UInteger> x, BSNumber<UInteger> y)
        {
            return (BSNumber<UInteger>)pow((double)x, (double)y);
        }

        inline static BSNumber<UInteger> Sin(BSNumber<UInteger> x)
        {
            return (BSNumber<UInteger>)sin((double)x);
        }

        inline static BSNumber<UInteger> Sinh(BSNumber<UInteger> x)
        {
            return (BSNumber<UInteger>)sinh((double)x);
        }

        inline static BSNumber<UInteger> Sinpi(BSNumber<UInteger> x)
        {
            return (BSNumber<UInteger>)sin((double)(x * (BSNumber<UInteger>)GTE_C_PI));
        }

        inline static BSNumber<UInteger> Sqr(BSNumber<UInteger> x)
        {
            return x * x;
        }

        inline static BSNumber<UInteger> Sqrt(BSNumber<UInteger> x)
        {
            return (BSNumber<UInteger>)sqrt((double)x);
        }

        inline static BSNumber<UInteger> Tan(BSNumber<UInteger> x)
        {
            return (BSNumber<UInteger>)tan((double)x);
        }

        inline static BSNumber<UInteger> Tanh(BSNumber<UInteger> x)
        {
            return (BSNumber<UInteger>)tanh((double)x);
        }

        inline static BSNumber<UInteger> Sign(BSNumber<UInteger> x)
        {
            return (BSNumber<UInteger>)x.GetSign();
        }

        inline static int ISign(BSNumber<UInteger> x)
        {
            return x.GetSign();
        }

        inline static BSNumber<UInteger> Clamp(BSNumber<UInteger> x, BSNumber<UInteger> min, BSNumber<UInteger> max)
        {
            return (x <= min ? min : (x >= max ? max : x));
        }

        inline static BSNumber<UInteger> Saturate(BSNumber<UInteger> x)
        {
            BSNumber<UInteger> const zero(0), one(1);
            return (x <= zero ? zero : (x >= one ? one : x));
        }

        inline static unsigned int GetMaxBisections()
        {
            return std::numeric_limits<unsigned int>::max();
        }
    };
}
