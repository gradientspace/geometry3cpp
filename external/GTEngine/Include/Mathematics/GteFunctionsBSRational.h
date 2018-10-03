// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.9.0 (2017/06/18)

#pragma once

#include <Mathematics/GteFunctions.h>
#include <Mathematics/GteBSRational.h>

namespace gte
{
    template <typename UInteger>
    class Function<BSRational<UInteger>>
    {
    public:
        inline static BSRational<UInteger> ACos(BSRational<UInteger> const& x)
        {
            return (BSRational<UInteger>)acos((double)x);
        }

        inline static BSRational<UInteger> ACosh(BSRational<UInteger> x)
        {
            double y = (double)x;
            return (BSRational<UInteger>)log(y + sqrt(y * y - 1.0));
        }

        inline static BSRational<UInteger> ASin(BSRational<UInteger> x)
        {
            return (BSRational<UInteger>)asin((double)x);
        }

        inline static BSRational<UInteger> ASinh(BSRational<UInteger> x)
        {
            double y = (double)x;
            return (BSRational<UInteger>)log(y + sqrt(y * y + 1.0));
        }

        inline static BSRational<UInteger> ATan(BSRational<UInteger> x)
        {
            return (BSRational<UInteger>)atan((double)x);
        }

        inline static BSRational<UInteger> ATanh(BSRational<UInteger> x)
        {
            double y = (double)x;
            return (BSRational<UInteger>)(log((1.0 + y) / (1.0 - y)) * 0.5);
        }

        inline static BSRational<UInteger> ATan2(BSRational<UInteger> y, BSRational<UInteger> x)
        {
            return (BSRational<UInteger>)atan2((double)y, (double)x);

        }

        inline static BSRational<UInteger> ATanpi(BSRational<UInteger> x)
        {
            return (BSRational<UInteger>)(atan((double)x) * GTE_C_INV_PI);
        }

        inline static BSRational<UInteger> ATan2pi(BSRational<UInteger> y, BSRational<UInteger> x)
        {
            return (BSRational<UInteger>)(atan2((double)y, (double)x) * GTE_C_INV_PI);
        }

        inline static BSRational<UInteger> Ceil(BSRational<UInteger> x)
        {
            return (BSRational<UInteger>)ceil((double)x);
        }

        inline static BSRational<UInteger> Cos(BSRational<UInteger> x)
        {
            return (BSRational<UInteger>)cos((double)x);
        }

        inline static BSRational<UInteger> Cosh(BSRational<UInteger> x)
        {
            return (BSRational<UInteger>)cosh((double)x);
        }

        inline static BSRational<UInteger> Cospi(BSRational<UInteger> x)
        {
            return (BSRational<UInteger>)cos((double)(x * (BSRational<UInteger>)GTE_C_PI));
        }

        inline static BSRational<UInteger> Exp(BSRational<UInteger> x)
        {
            return (BSRational<UInteger>)exp((double)x);
        }

        inline static BSRational<UInteger> Exp2(BSRational<UInteger> x)
        {
            return (BSRational<UInteger>)exp((double)(x * (BSRational<UInteger>)GTE_C_LN_2));
        }

        inline static BSRational<UInteger> Exp10(BSRational<UInteger> x)
        {
            return (BSRational<UInteger>)exp((double)(x * (BSRational<UInteger>)GTE_C_LN_10));
        }

        inline static BSRational<UInteger> FAbs(BSRational<UInteger> x)
        {
            return (x.GetSign() >= 0 ? x : -x);
        }

        inline static BSRational<UInteger> Floor(BSRational<UInteger> x)
        {
            return (BSRational<UInteger>)floor((double)x);
        }

        inline static BSRational<UInteger> FMod(BSRational<UInteger> x, BSRational<UInteger> y)
        {
            return (BSRational<UInteger>)fmod((double)x, (double)y);
        }

        inline static BSRational<UInteger> FRExp(BSRational<UInteger> x, int* exponent)
        {
            return (BSRational<UInteger>)frexp((double)x, exponent);
        }

        inline static BSRational<UInteger> InvSqrt(BSRational<UInteger> x)
        {
            return (BSRational<UInteger>)(1.0 / sqrt((double)x));
        }

        inline static BSRational<UInteger> LDExp(BSRational<UInteger> x, int exponent)
        {
            return (BSRational<UInteger>)ldexp((double)x, exponent);
        }

        inline static BSRational<UInteger> Log(BSRational<UInteger> x)
        {
            return (BSRational<UInteger>)log((double)x);
        }

        inline static BSRational<UInteger> Log2(BSRational<UInteger> x)
        {
            return (BSRational<UInteger>)(log((double)x) * GTE_C_INV_LN_2);
        }

        inline static BSRational<UInteger> Log10(BSRational<UInteger> x)
        {
            return (BSRational<UInteger>)log10((double)x);
        }

        inline static BSRational<UInteger> Pow(BSRational<UInteger> x, BSRational<UInteger> y)
        {
            return (BSRational<UInteger>)pow((double)x, (double)y);
        }

        inline static BSRational<UInteger> Sin(BSRational<UInteger> x)
        {
            return (BSRational<UInteger>)sin((double)x);
        }

        inline static BSRational<UInteger> Sinh(BSRational<UInteger> x)
        {
            return (BSRational<UInteger>)sinh((double)x);
        }

        inline static BSRational<UInteger> Sinpi(BSRational<UInteger> x)
        {
            return (BSRational<UInteger>)sin((double)(x * (BSRational<UInteger>)GTE_C_PI));
        }

        inline static BSRational<UInteger> Sqr(BSRational<UInteger> x)
        {
            return x * x;
        }

        inline static BSRational<UInteger> Sqrt(BSRational<UInteger> x)
        {
            return (BSRational<UInteger>)sqrt((double)x);
        }

        inline static BSRational<UInteger> Tan(BSRational<UInteger> x)
        {
            return (BSRational<UInteger>)tan((double)x);
        }

        inline static BSRational<UInteger> Tanh(BSRational<UInteger> x)
        {
            return (BSRational<UInteger>)tanh((double)x);
        }

        inline static BSRational<UInteger> Sign(BSRational<UInteger> x)
        {
            return (BSRational<UInteger>)x.GetSign();
        }

        inline static int ISign(BSRational<UInteger> x)
        {
            return x.GetSign();
        }

        inline static BSRational<UInteger> Clamp(BSRational<UInteger> x, BSRational<UInteger> min, BSRational<UInteger> max)
        {
            return (x <= min ? min : (x >= max ? max : x));
        }

        inline static BSRational<UInteger> Saturate(BSRational<UInteger> x)
        {
            BSRational<UInteger> const zero(0), one(1);
            return (x <= zero ? zero : (x >= one ? one : x));
        }

        inline static unsigned int GetMaxBisections()
        {
            return std::numeric_limits<unsigned int>::max();
        }
    };
}
