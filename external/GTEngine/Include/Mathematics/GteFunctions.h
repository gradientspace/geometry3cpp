// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <Mathematics/GteArithmetic.h>
#include <Mathematics/GteConstants.h>
#include <limits>

namespace gte
{

template <typename Real>
class Function
{
public:
    static Real ACos(Real const& x);                    // acos(x)
    static Real ACosh(Real const& x);                   // acosh(x)
    static Real ASin(Real const& x);                    // asin(x)
    static Real ASinh(Real const& x);                   // asinh(x)
    static Real ATan(Real const& x);                    // atan(x)
    static Real ATanh(Real const& x);                   // atanh(x)
    static Real ATan2(Real const& y, Real const& x);    // atan2(y, x)
    static Real ATanpi(Real const& x);                  // atan(x)/pi
    static Real ATan2pi(Real const& y, Real const& x);  // atan2(y,x)/pi
    static Real Ceil(Real const& x);                    // ceil(x)
    static Real Cos(Real const& x);                     // cos(x)
    static Real Cosh(Real const& x);                    // cosh(x)
    static Real Cospi(Real const& x);                   // cos(pi*x)
    static Real Exp(Real const& x);                     // e^x
    static Real Exp2(Real const& x);                    // 2^x
    static Real Exp10(Real const& x);                   // 10^x
    static Real FAbs(Real const& x);                    // |x|
    static Real Floor(Real const& x);                   // floor(x)
    static Real FMod(Real const& x, Real const& y);     // fmod(x,y)
    static Real InvSqrt(Real const& x);                 // x^{-1/2}
    static Real Log(Real const& x);                     // log_e(x)
    static Real Log2(Real const& x);                    // log_2(x)
    static Real Log10(Real const& x);                   // log_10(x)
    static Real Pow(Real const& x, Real const& y);      // x^y
    static Real Sin(Real const& x);                     // sin(x)
    static Real Sinh(Real const& x);                    // sinh(x)
    static Real Sinpi(Real const& x);                   // sin(pi*x)
    static Real Sqr(Real const& x);                     // x^2
    static Real Sqrt(Real const& x);                    // x^{1/2}
    static Real Tan(Real const& x);                     // tan(x)
    static Real Tanh(Real const& x);                    // tanh(x)

    static Real Sign(Real const& x);    // sign of x as a Real number
    static int ISign(Real const& x);    // sign of x as an integer

    // Clamp x to the interval [min,max].
    static Real Clamp(Real const& x, Real const& min, Real const & max);

    // Clamp x to the interval [0,1].
    static Real Saturate(Real const& x);

    // The maximum number of bisections required for root finding on an
    // interval of finite floating-point numbers.  After this number of
    // iterations, the root-bounding interval consists of two consecutive
    // floating-point numbers.  Thus, you have reached the limit of the
    // precision of the numbers.  WARNING:  The return value for BSNumber
    // or BSRational types is std::numeric_limits<unsigned int>::max(),
    // because there is no limit on precision (other than the practical
    // limitations of memory usage and computational costs).  We recommend
    // you use a reasonable number of the BS* types.
    static unsigned int GetMaxBisections();

private:
    // The tag-dispatch pattern is used for template-parameter-controlled
    // instantiation of mathematics functions.

    template <typename T> static T ACosImpl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T ACosImpl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T ACosImpl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T ACoshImpl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T ACoshImpl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T ACoshImpl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T ASinImpl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T ASinImpl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T ASinImpl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T ASinhImpl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T ASinhImpl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T ASinhImpl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T ATanImpl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T ATanImpl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T ATanImpl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T ATanhImpl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T ATanhImpl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T ATanhImpl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T ATan2Impl(T y, T x, Arithmetic::IsFPType tag);
    template <typename T> static T ATan2Impl(T y, T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T ATan2Impl(T const& y, T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T ATanpiImpl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T ATanpiImpl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T ATanpiImpl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T ATan2piImpl(T y, T x, Arithmetic::IsFPType tag);
    template <typename T> static T ATan2piImpl(T y, T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T ATan2piImpl(T const& y, T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T CeilImpl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T CeilImpl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T CeilImpl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T CosImpl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T CosImpl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T CosImpl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T CoshImpl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T CoshImpl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T CoshImpl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T CospiImpl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T CospiImpl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T CospiImpl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T ExpImpl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T ExpImpl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T ExpImpl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T Exp2Impl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T Exp2Impl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T Exp2Impl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T Exp10Impl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T Exp10Impl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T Exp10Impl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T FAbsImpl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T FAbsImpl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T FAbsImpl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T FloorImpl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T FloorImpl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T FloorImpl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T FModImpl(T x, T y, Arithmetic::IsFPType tag);
    template <typename T> static T FModImpl(T x, T y, Arithmetic::IsFP16Type tag);
    template <typename T> static T FModImpl(T const& x, T const& y, Arithmetic::IsBSType tag);

    template <typename T> static T InvSqrtImpl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T InvSqrtImpl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T InvSqrtImpl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T LogImpl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T LogImpl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T LogImpl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T Log2Impl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T Log2Impl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T Log2Impl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T Log10Impl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T Log10Impl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T Log10Impl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T PowImpl(T x, T y, Arithmetic::IsFPType tag);
    template <typename T> static T PowImpl(T x, T y, Arithmetic::IsFP16Type tag);
    template <typename T> static T PowImpl(T const& x, T const& y, Arithmetic::IsBSType tag);

    template <typename T> static T SinImpl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T SinImpl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T SinImpl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T SinhImpl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T SinhImpl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T SinhImpl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T SinpiImpl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T SinpiImpl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T SinpiImpl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T SqrtImpl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T SqrtImpl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T SqrtImpl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T TanImpl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T TanImpl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T TanImpl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T TanhImpl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T TanhImpl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T TanhImpl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T SignImpl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T SignImpl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T SignImpl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static int ISignImpl(T x, Arithmetic::IsFPType tag);
    template <typename T> static int ISignImpl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static int ISignImpl(T const& x, Arithmetic::IsBSType tag);

    template <typename T> static T SaturateImpl(T x, Arithmetic::IsFPType tag);
    template <typename T> static T SaturateImpl(T x, Arithmetic::IsFP16Type tag);
    template <typename T> static T SaturateImpl(T const& x, Arithmetic::IsBSType tag);

    template <typename T>
    static unsigned int GetMaxBisectionsImpl(Arithmetic::IsFPType tag);

    template <typename T>
    static unsigned int GetMaxBisectionsImpl(Arithmetic::IsFP16Type tag);

    template <typename T>
    static unsigned int GetMaxBisectionsImpl(Arithmetic::IsBSType tag);
};


template <typename Real>
Real Function<Real>::ACos(Real const& x)
{
    return ACosImpl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::ACosImpl(T x, Arithmetic::IsFPType)
{
    return acos(x);
}

template <typename Real> template <typename T>
T Function<Real>::ACosImpl(T x, Arithmetic::IsFP16Type)
{
    return (T)acos((float)x);
}

template <typename Real> template <typename T>
T Function<Real>::ACosImpl(T const& x, Arithmetic::IsBSType)
{
    return (T)acos((double)x);
}



template <typename Real>
Real Function<Real>::ACosh(Real const& x)
{
    return ACoshImpl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::ACoshImpl(T x, Arithmetic::IsFPType)
{
    return log(x + sqrt(x * x - (Real)1));
}

template <typename Real> template <typename T>
T Function<Real>::ACoshImpl(T x, Arithmetic::IsFP16Type)
{
    float y = (float)x;
    return (T)log(y + sqrt(y * y - 1.0f));
}

template <typename Real> template <typename T>
T Function<Real>::ACoshImpl(T const& x, Arithmetic::IsBSType)
{
    double y = (double)x;
    return (T)log(y + sqrt(y * y - 1.0));
}



template <typename Real>
Real Function<Real>::ASin(Real const& x)
{
    return ASinImpl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::ASinImpl(T x, Arithmetic::IsFPType)
{
    return asin(x);
}

template <typename Real> template <typename T>
T Function<Real>::ASinImpl(T x, Arithmetic::IsFP16Type)
{
    return (T)asin((float)x);
}

template <typename Real> template <typename T>
T Function<Real>::ASinImpl(T const& x, Arithmetic::IsBSType)
{
    return (T)asin((double)x);
}



template <typename Real>
Real Function<Real>::ASinh(Real const& x)
{
    return ASinhImpl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::ASinhImpl(T x, Arithmetic::IsFPType)
{
    return log(x + sqrt(x * x + (Real)1));
}

template <typename Real> template <typename T>
T Function<Real>::ASinhImpl(T x, Arithmetic::IsFP16Type)
{
    float y = (float)x;
    return (T)log(y + sqrt(y * y + 1.0f));
}

template <typename Real> template <typename T>
T Function<Real>::ASinhImpl(T const& x, Arithmetic::IsBSType)
{
    double y = (double)x;
    return (T)log(y + sqrt(y * y + 1.0));
}



template <typename Real>
Real Function<Real>::ATan(Real const& x)
{
    return ATanImpl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::ATanImpl(T x, Arithmetic::IsFPType)
{
    return atan(x);
}

template <typename Real> template <typename T>
T Function<Real>::ATanImpl(T x, Arithmetic::IsFP16Type)
{
    return (T)atan((float)x);
}

template <typename Real> template <typename T>
T Function<Real>::ATanImpl(T const& x, Arithmetic::IsBSType)
{
    return (T)atan((double)x);
}



template <typename Real>
Real Function<Real>::ATanh(Real const& x)
{
    return ATanhImpl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::ATanhImpl(T x, Arithmetic::IsFPType)
{
    return log( ((Real)1 + x) / ((Real)1 - x) ) * (Real)0.5;
}

template <typename Real> template <typename T>
T Function<Real>::ATanhImpl(T x, Arithmetic::IsFP16Type)
{
    float y = (float)x;
    return (T)(log( (1.0f + y) / (1.0f - y) ) * 0.5f);
}

template <typename Real> template <typename T>
T Function<Real>::ATanhImpl(T const& x, Arithmetic::IsBSType)
{
    double y = (double)x;
    return (T)(log( (1.0 + y) / (1.0 - y) ) * 0.5);
}



template <typename Real>
Real Function<Real>::ATan2(Real const& y, Real const& x)
{
    return ATan2Impl<Real>(y, x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::ATan2Impl(T y, T x, Arithmetic::IsFPType)
{
    return atan2(y, x);
}

template <typename Real> template <typename T>
T Function<Real>::ATan2Impl(T y, T x, Arithmetic::IsFP16Type)
{
    return (T)atan2((float)y, (float)x);
}

template <typename Real> template <typename T>
T Function<Real>::ATan2Impl(T const& y, T const& x, Arithmetic::IsBSType)
{
    return (T)atan2((double)y, (double)x);
}



template <typename Real>
Real Function<Real>::ATanpi(Real const& x)
{
    return ATanpiImpl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::ATanpiImpl(T x, Arithmetic::IsFPType)
{
    return atan(x) * (Real)GTE_C_INV_PI;
}

template <typename Real> template <typename T>
T Function<Real>::ATanpiImpl(T x, Arithmetic::IsFP16Type)
{
    return (T)(atan((float)x) * (float)GTE_C_INV_PI);
}

template <typename Real> template <typename T>
T Function<Real>::ATanpiImpl(T const& x, Arithmetic::IsBSType)
{
    return (T)(atan((double)x) * GTE_C_INV_PI);
}



template <typename Real>
Real Function<Real>::ATan2pi(Real const& y, Real const& x)
{
    return ATan2piImpl<Real>(y, x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::ATan2piImpl(T y, T x, Arithmetic::IsFPType)
{
    return atan2(y, x) * (Real)GTE_C_INV_PI;
}

template <typename Real> template <typename T>
T Function<Real>::ATan2piImpl(T y, T x, Arithmetic::IsFP16Type)
{
    return (T)(atan2((float)y, (float)x) * (float)GTE_C_INV_PI);
}

template <typename Real> template <typename T>
T Function<Real>::ATan2piImpl(T const& y, T const& x, Arithmetic::IsBSType)
{
    return (T)(atan2((double)y, (double)x) * GTE_C_INV_PI);
}



template <typename Real>
Real Function<Real>::Ceil(Real const& x)
{
    return CeilImpl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::CeilImpl(T x, Arithmetic::IsFPType)
{
    return ceil(x);
}

template <typename Real> template <typename T>
T Function<Real>::CeilImpl(T x, Arithmetic::IsFP16Type)
{
    return (T)ceil((float)x);
}

template <typename Real> template <typename T>
T Function<Real>::CeilImpl(T const& x, Arithmetic::IsBSType)
{
    return (T)ceil((double)x);
}



template <typename Real>
Real Function<Real>::Cos(Real const& x)
{
    return CosImpl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::CosImpl(T x, Arithmetic::IsFPType)
{
    return cos(x);
}

template <typename Real> template <typename T>
T Function<Real>::CosImpl(T x, Arithmetic::IsFP16Type)
{
    return (T)cos((float)x);
}

template <typename Real> template <typename T>
T Function<Real>::CosImpl(T const& x, Arithmetic::IsBSType)
{
    return (T)cos((double)x);
}



template <typename Real>
Real Function<Real>::Cosh(Real const& x)
{
    return CoshImpl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::CoshImpl(T x, Arithmetic::IsFPType)
{
    return cosh(x);
}

template <typename Real> template <typename T>
T Function<Real>::CoshImpl(T x, Arithmetic::IsFP16Type)
{
    return (T)cosh((float)x);
}

template <typename Real> template <typename T>
T Function<Real>::CoshImpl(T const& x, Arithmetic::IsBSType)
{
    return (T)cosh((double)x);
}



template <typename Real>
Real Function<Real>::Cospi(Real const& x)
{
    return CospiImpl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::CospiImpl(T x, Arithmetic::IsFPType)
{
    return cos(x * (Real)GTE_C_PI);
}

template <typename Real> template <typename T>
T Function<Real>::CospiImpl(T x, Arithmetic::IsFP16Type)
{
    return (T)cos((float)x * (float)GTE_C_PI);
}

template <typename Real> template <typename T>
T Function<Real>::CospiImpl(T const& x, Arithmetic::IsBSType)
{
    return (T)cos((double)(x * (T)GTE_C_PI));
}



template <typename Real>
Real Function<Real>::Exp(Real const& x)
{
    return ExpImpl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::ExpImpl(T x, Arithmetic::IsFPType)
{
    return exp(x);
}

template <typename Real> template <typename T>
T Function<Real>::ExpImpl(T x, Arithmetic::IsFP16Type)
{
    return (T)exp((float)x);
}

template <typename Real> template <typename T>
T Function<Real>::ExpImpl(T const& x, Arithmetic::IsBSType)
{
    return (T)exp((double)x);
}



template <typename Real>
Real Function<Real>::Exp2(Real const& x)
{
    return Exp2Impl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::Exp2Impl(T x, Arithmetic::IsFPType)
{
    return exp(x * (Real)GTE_C_LN_2);
}

template <typename Real> template <typename T>
T Function<Real>::Exp2Impl(T x, Arithmetic::IsFP16Type)
{
    return (T)exp((float)x * (float)GTE_C_LN_2);
}

template <typename Real> template <typename T>
T Function<Real>::Exp2Impl(T const& x, Arithmetic::IsBSType)
{
    return (T)exp((double)(x * (T)GTE_C_LN_2));
}



template <typename Real>
Real Function<Real>::Exp10(Real const& x)
{
    return Exp10Impl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::Exp10Impl(T x, Arithmetic::IsFPType)
{
    return exp(x * (Real)GTE_C_LN_10);
}

template <typename Real> template <typename T>
T Function<Real>::Exp10Impl(T x, Arithmetic::IsFP16Type)
{
    return (T)exp((float)x * (float)GTE_C_LN_10);
}

template <typename Real> template <typename T>
T Function<Real>::Exp10Impl(T const& x, Arithmetic::IsBSType)
{
    return (T)exp((double)(x * (T)GTE_C_LN_10));
}



template <typename Real>
Real Function<Real>::FAbs(Real const& x)
{
    return FAbsImpl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::FAbsImpl(T x, Arithmetic::IsFPType)
{
    return fabs(x);
}

template <typename Real> template <typename T>
T Function<Real>::FAbsImpl(T x, Arithmetic::IsFP16Type)
{
    return (Real)fabs((float)x);
}

template <typename Real> template <typename T>
T Function<Real>::FAbsImpl(T const& x, Arithmetic::IsBSType)
{
    return (x.GetSign() >= 0 ? x : -x);
}



template <typename Real>
Real Function<Real>::Floor(Real const& x)
{
    return FloorImpl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::FloorImpl(T x, Arithmetic::IsFPType)
{
    return floor(x);
}

template <typename Real> template <typename T>
T Function<Real>::FloorImpl(T x, Arithmetic::IsFP16Type)
{
    return (T)floor((float)x);
}

template <typename Real> template <typename T>
T Function<Real>::FloorImpl(T const& x, Arithmetic::IsBSType)
{
    return (T)floor((double)x);
}



template <typename Real>
Real Function<Real>::FMod(Real const& x, Real const& y)
{
    return FModImpl<Real>(x, y, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::FModImpl(T x, T y, Arithmetic::IsFPType)
{
    return fmod(x, y);
}

template <typename Real> template <typename T>
T Function<Real>::FModImpl(T x, T y, Arithmetic::IsFP16Type)
{
    return (T)fmod((float)x, (float)y);
}

template <typename Real> template <typename T>
T Function<Real>::FModImpl(T const& x, T const& y, Arithmetic::IsBSType)
{
    return (T)fmod((double)x, (double)y);
}



template <typename Real>
Real Function<Real>::InvSqrt(Real const& x)
{
    return InvSqrtImpl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::InvSqrtImpl(T x, Arithmetic::IsFPType)
{
    return ((T)1) / sqrt(x);
}

template <typename Real> template <typename T>
T Function<Real>::InvSqrtImpl(T x, Arithmetic::IsFP16Type)
{
    return (T)(1.0f / sqrt((float)x));
}

template <typename Real> template <typename T>
T Function<Real>::InvSqrtImpl(T const& x, Arithmetic::IsBSType)
{
    return (T)(1.0 / sqrt((double)x));
}



template <typename Real>
Real Function<Real>::Log(Real const& x)
{
    return LogImpl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::LogImpl(T x, Arithmetic::IsFPType)
{
    return log(x);
}

template <typename Real> template <typename T>
T Function<Real>::LogImpl(T x, Arithmetic::IsFP16Type)
{
    return (T)log((float)x);
}

template <typename Real> template <typename T>
T Function<Real>::LogImpl(T const& x, Arithmetic::IsBSType)
{
    return (T)log((double)x);
}



template <typename Real>
Real Function<Real>::Log2(Real const& x)
{
    return Log2Impl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::Log2Impl(T x, Arithmetic::IsFPType)
{
    return log(x) * (Real)GTE_C_INV_LN_2;
}

template <typename Real> template <typename T>
T Function<Real>::Log2Impl(T x, Arithmetic::IsFP16Type)
{
    return (T)(log((float)x) * (float)GTE_C_INV_LN_2);
}

template <typename Real> template <typename T>
T Function<Real>::Log2Impl(T const& x, Arithmetic::IsBSType)
{
    return (T)(log((double)x) * GTE_C_INV_LN_2);
}



template <typename Real>
Real Function<Real>::Log10(Real const& x)
{
    return Log10Impl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::Log10Impl(T x, Arithmetic::IsFPType)
{
    return log10(x);
}

template <typename Real> template <typename T>
T Function<Real>::Log10Impl(T x, Arithmetic::IsFP16Type)
{
    return (T)log10((float)x);
}

template <typename Real> template <typename T>
T Function<Real>::Log10Impl(T const& x, Arithmetic::IsBSType)
{
    return (T)log10((double)x);
}



template <typename Real>
Real Function<Real>::Pow(Real const& x, Real const& y)
{
    return PowImpl<Real>(x, y, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::PowImpl(T x, T y, Arithmetic::IsFPType)
{
    return pow(x, y);
}

template <typename Real> template <typename T>
T Function<Real>::PowImpl(T x, T y, Arithmetic::IsFP16Type)
{
    return (T)pow((float)x, (float)y);
}

template <typename Real> template <typename T>
T Function<Real>::PowImpl(T const& x, T const& y, Arithmetic::IsBSType)
{
    return (T)pow((double)x, (double)y);
}



template <typename Real>
Real Function<Real>::Sin(Real const& x)
{
    return SinImpl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::SinImpl(T x, Arithmetic::IsFPType)
{
    return sin(x);
}

template <typename Real> template <typename T>
T Function<Real>::SinImpl(T x, Arithmetic::IsFP16Type)
{
    return (T)sin((float)x);
}

template <typename Real> template <typename T>
T Function<Real>::SinImpl(T const& x, Arithmetic::IsBSType)
{
    return (T)sin((double)x);
}



template <typename Real>
Real Function<Real>::Sinh(Real const& x)
{
    return SinhImpl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::SinhImpl(T x, Arithmetic::IsFPType)
{
    return sinh(x);
}

template <typename Real> template <typename T>
T Function<Real>::SinhImpl(T x, Arithmetic::IsFP16Type)
{
    return (T)sinh((float)x);
}

template <typename Real> template <typename T>
T Function<Real>::SinhImpl(T const& x, Arithmetic::IsBSType)
{
    return (T)sinh((double)x);
}



template <typename Real>
Real Function<Real>::Sinpi(Real const& x)
{
    return SinpiImpl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::SinpiImpl(T x, Arithmetic::IsFPType)
{
    return sin(x * (Real)GTE_C_PI);
}

template <typename Real> template <typename T>
T Function<Real>::SinpiImpl(T x, Arithmetic::IsFP16Type)
{
    return (T)sin((float)x * (float)GTE_C_PI);
}

template <typename Real> template <typename T>
T Function<Real>::SinpiImpl(T const& x, Arithmetic::IsBSType)
{
    return (T)sin((double)(x * (T)GTE_C_PI));
}



template <typename Real>
Real Function<Real>::Sqr(Real const& x)
{
    return x * x;
}



template <typename Real>
Real Function<Real>::Sqrt(Real const& x)
{
    return SqrtImpl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::SqrtImpl(T x, Arithmetic::IsFPType)
{
    return sqrt(x);
}

template <typename Real> template <typename T>
T Function<Real>::SqrtImpl(T x, Arithmetic::IsFP16Type)
{
    return (T)sqrt((float)x);
}

template <typename Real> template <typename T>
T Function<Real>::SqrtImpl(T const& x, Arithmetic::IsBSType)
{
    return (T)sqrt((double)x);
}



template <typename Real>
Real Function<Real>::Tan(Real const& x)
{
    return TanImpl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::TanImpl(T x, Arithmetic::IsFPType)
{
    return tan(x);
}

template <typename Real> template <typename T>
T Function<Real>::TanImpl(T x, Arithmetic::IsFP16Type)
{
    return (T)tan((float)x);
}

template <typename Real> template <typename T>
T Function<Real>::TanImpl(T const& x, Arithmetic::IsBSType)
{
    return (T)tan((double)x);
}



template <typename Real>
Real Function<Real>::Tanh(Real const& x)
{
    return TanhImpl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::TanhImpl(T x, Arithmetic::IsFPType)
{
    return tanh(x);
}

template <typename Real> template <typename T>
T Function<Real>::TanhImpl(T x, Arithmetic::IsFP16Type)
{
    return (T)tanh((float)x);
}

template <typename Real> template <typename T>
T Function<Real>::TanhImpl(T const& x, Arithmetic::IsBSType)
{
    return (T)tanh((double)x);
}



template <typename Real>
Real Function<Real>::Sign(Real const& x)
{
    return SignImpl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::SignImpl(T x, Arithmetic::IsFPType)
{
    return (x > (Real)0 ? (Real)1 : (x < (Real)0 ? (Real)-1 : (Real)0));
}

template <typename Real> template <typename T>
T Function<Real>::SignImpl(T x, Arithmetic::IsFP16Type)
{
    float y = (Real)x;
    return (y > 0.0f ? 1.0f : (y < 0.0f ? -1.0f : 0.0f));
}

template <typename Real> template <typename T>
T Function<Real>::SignImpl(T const& x, Arithmetic::IsBSType)
{
    return (Real)x.GetSign();
}



template <typename Real>
int Function<Real>::ISign(Real const& x)
{
    return ISignImpl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
int Function<Real>::ISignImpl(T x, Arithmetic::IsFPType)
{
    return (x > (Real)0 ? 1 : (x < (Real)0 ? -1 : 0));
}

template <typename Real> template <typename T>
int Function<Real>::ISignImpl(T x, Arithmetic::IsFP16Type)
{
    float y = (Real)x;
    return (y > 0.0f ? 1 : (y < 0 ? -1 : 0));
}

template <typename Real> template <typename T>
int Function<Real>::ISignImpl(T const& x, Arithmetic::IsBSType)
{
    return x.GetSign();
}



template <typename Real>
Real Function<Real>::Clamp(Real const& x, Real const& min, Real const& max)
{
    return (x <= min ? min : (x >= max ? max : x));
}

template <typename Real>
Real Function<Real>::Saturate(Real const& x)
{
    return SaturateImpl<Real>(x, Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
T Function<Real>::SaturateImpl(T x, Arithmetic::IsFPType)
{
    return (x <= (Real)0 ? (Real)0 : (x >= (Real)1 ? (Real)1 : x));
}

template <typename Real> template <typename T>
T Function<Real>::SaturateImpl(T x, Arithmetic::IsFP16Type)
{
    float y = (float)x;
    return (y <= 0.0f ? 0.0f : (y >= 1.0f ? 1.0f : y));
}

template <typename Real> template <typename T>
T Function<Real>::SaturateImpl(T const& x, Arithmetic::IsBSType)
{
    Real const zero(0), one(1);
    return (x <= zero ? zero : (x >= one ? one : x));
}



template <typename Real>
unsigned int Function<Real>::GetMaxBisections()
{
    return GetMaxBisectionsImpl<Real>(Arithmetic::WhichType<Real>());
}

template <typename Real> template <typename T>
unsigned int Function<Real>::GetMaxBisectionsImpl(Arithmetic::IsFPType)
{
    return (unsigned int)(3 + std::numeric_limits<Real>::digits -
        std::numeric_limits<Real>::min_exponent);
}

template <typename Real> template <typename T>
unsigned int Function<Real>::GetMaxBisectionsImpl(Arithmetic::IsFP16Type)
{
    // std::numeric_limits<IEEEBinary16>::digits = 11
    // std::numeric_limits<IEEEBinary16>::min_exponent = -13
    // 27 = 3 + digits - min_exponent
    return 27;
}

template <typename Real> template <typename T>
unsigned int Function<Real>::GetMaxBisectionsImpl(Arithmetic::IsBSType)
{
    return std::numeric_limits<unsigned int>::max();
}


}
