// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.9.0 (2017/06/18)

#pragma once

#include <Mathematics/GteFunctions.h>
#include <Mathematics/GteIEEEBinary16.h>

namespace gte
{
    template <>
    class Function<IEEEBinary16>
    {
    public:
        inline static IEEEBinary16 ACos(IEEEBinary16 x)
        {
            return (IEEEBinary16)acos((float)x);
        }

        inline static IEEEBinary16 ACosh(IEEEBinary16 x)
        {
            float y = (float)x;
            return (IEEEBinary16)log(y + sqrt(y * y - 1.0f));
        }

        inline static IEEEBinary16 ASin(IEEEBinary16 x)
        {
            return (IEEEBinary16)asin((float)x);
        }

        inline static IEEEBinary16 ASinh(IEEEBinary16 x)
        {
            float y = (float)x;
            return (IEEEBinary16)log(y + sqrt(y * y + 1.0f));
        }

        inline static IEEEBinary16 ATan(IEEEBinary16 x)
        {
            return (IEEEBinary16)atan((float)x);
        }

        inline static IEEEBinary16 ATanh(IEEEBinary16 x)
        {
            float y = (float)x;
            return (IEEEBinary16)(log((1.0f + y) / (1.0f - y)) * 0.5f);
        }

        inline static IEEEBinary16 ATan2(IEEEBinary16 y, IEEEBinary16 x)
        {
            return (IEEEBinary16)atan2((float)y, (float)x);
        }

        inline static IEEEBinary16 ATanpi(IEEEBinary16 x)
        {
            return (IEEEBinary16)(atan((float)x) * (float)GTE_C_INV_PI);
        }

        inline static IEEEBinary16 ATan2pi(IEEEBinary16 y, IEEEBinary16 x)
        {
            return (IEEEBinary16)(atan2((float)y, (float)x) * (float)GTE_C_INV_PI);
        }

        inline static IEEEBinary16 Ceil(IEEEBinary16 x)
        {
            return (IEEEBinary16)ceil((float)x);
        }

        inline static IEEEBinary16 Cos(IEEEBinary16 x)
        {
            return (IEEEBinary16)cos((float)x);
        }

        inline static IEEEBinary16 Cosh(IEEEBinary16 x)
        {
            return (IEEEBinary16)cosh((float)x);
        }

        inline static IEEEBinary16 Cospi(IEEEBinary16 x)
        {
            return (IEEEBinary16)cos((float)x * (float)GTE_C_PI);
        }

        static IEEEBinary16 Exp(IEEEBinary16 x)
        {
            return (IEEEBinary16)exp((float)x);
        }

        static IEEEBinary16 Exp2(IEEEBinary16 x)
        {
            return (IEEEBinary16)exp((float)x * (float)GTE_C_LN_2);
        }

        static IEEEBinary16 Exp10(IEEEBinary16 x)
        {
            return (IEEEBinary16)exp((float)x * (float)GTE_C_LN_10);
        }

        inline static IEEEBinary16 FAbs(IEEEBinary16 x)
        {
            return (IEEEBinary16)fabs((float)x);
        }

        inline static IEEEBinary16 Floor(IEEEBinary16 x)
        {
            return (IEEEBinary16)floor((double)x);
        }

        inline static IEEEBinary16 FMod(IEEEBinary16 x, IEEEBinary16 y)
        {
            return (IEEEBinary16)fmod((float)x, (float)y);
        }

        inline static IEEEBinary16 FRExp(IEEEBinary16 x, int* exponent)
        {
            return (IEEEBinary16)frexp((float)x, exponent);
        }

        inline static IEEEBinary16 InvSqrt(IEEEBinary16 x)
        {
            return (IEEEBinary16)(1.0f / sqrt((float)x));
        }

        inline static IEEEBinary16 LDExp(IEEEBinary16 x, int exponent)
        {
            return (IEEEBinary16)ldexp((float)x, exponent);
        }

        inline static IEEEBinary16 Log(IEEEBinary16 x)
        {
            return (IEEEBinary16)log((float)x);
        }

        inline static IEEEBinary16 Log2(IEEEBinary16 x)
        {
            return (IEEEBinary16)(log((float)x) * (float)GTE_C_INV_LN_2);
        }

        inline static IEEEBinary16 Log10(IEEEBinary16 x)
        {
            return (IEEEBinary16)log10((float)x);
        }

        inline static IEEEBinary16 Pow(IEEEBinary16 x, IEEEBinary16 y)
        {
            return (IEEEBinary16)pow((float)x, (float)y);
        }

        inline static IEEEBinary16 Sin(IEEEBinary16 x)
        {
            return (IEEEBinary16)sin((float)x);
        }

        inline static IEEEBinary16 Sinh(IEEEBinary16 x)
        {
            return (IEEEBinary16)sinh((float)x);
        }

        inline static IEEEBinary16 Sinpi(IEEEBinary16 x)
        {
            return (IEEEBinary16)sin((float)x * (float)GTE_C_PI);
        }

        inline static IEEEBinary16 Sqr(IEEEBinary16 x)
        {
            return x * x;
        }

        inline static IEEEBinary16 Sqrt(IEEEBinary16 x)
        {
            return (IEEEBinary16)sqrt((float)x);
        }

        inline static IEEEBinary16 Tan(IEEEBinary16 x)
        {
            return (IEEEBinary16)tan((float)x);
        }

        inline static IEEEBinary16 Tanh(IEEEBinary16 x)
        {
            return (IEEEBinary16)tanh((float)x);
        }

        inline static IEEEBinary16 Sign(IEEEBinary16 x)
        {
            float y = (IEEEBinary16)x;
            return (y > 0.0f ? 1.0f : (y < 0.0f ? -1.0f : 0.0f));
        }

        inline static int ISign(IEEEBinary16 x)
        {
            float y = (IEEEBinary16)x;
            return (y > 0.0f ? 1 : (y < 0 ? -1 : 0));
        }

        inline static IEEEBinary16 Clamp(IEEEBinary16 x, IEEEBinary16 min, IEEEBinary16 max)
        {
            return (x <= min ? min : (x >= max ? max : x));
        }

        inline static IEEEBinary16 Saturate(IEEEBinary16 x)
        {
            float y = (float)x;
            return (y <= 0.0f ? 0.0f : (y >= 1.0f ? 1.0f : y));
        }

        inline static unsigned int GetMaxBisections()
        {
            // std::numeric_limits<IEEEBinary16>::digits = 11
            // std::numeric_limits<IEEEBinary16>::min_exponent = -13
            // 27 = 3 + digits - min_exponent
            return 27;
        }
    };
}
