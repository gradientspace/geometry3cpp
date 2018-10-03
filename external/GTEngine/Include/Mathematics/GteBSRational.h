// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.1 (2017/11/26)

#pragma once

#include <Mathematics/GteBSNumber.h>
#include <limits>

// See the comments in GteBSNumber.h about the UIntegerType requirements.  The
// denominator of a BSRational is chosen to be positive, which allows some
// simplification of comparisons.  Also see the comments about exposing the
// GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE conditional define.

namespace gte
{

template <typename UIntegerType>
class BSRational
{
public:
    // Construction.  The default constructor generates the zero BSRational.
    // The constructors that take only numerators set the denominators to one.
    BSRational();
    BSRational(BSRational const& r);
    BSRational(float numerator);
    BSRational(double numerator);
    BSRational(int32_t numerator);
    BSRational(uint32_t numerator);
    BSRational(int64_t numerator);
    BSRational(uint64_t numerator);
    BSRational(BSNumber<UIntegerType> const& numerator);
    BSRational(float numerator, float denominator);
    BSRational(double numerator, double denominator);
    BSRational(BSNumber<UIntegerType> const& numerator,
        BSNumber<UIntegerType> const& denominator);

    // Implicit conversions.
    operator float() const;
    operator double() const;

    // Assignment.
    BSRational& operator=(BSRational const& r);

    // Support for std::move.
    BSRational(BSRational&& r);
    BSRational& operator=(BSRational&& r);

    // Member access.
    inline int GetSign() const;
    inline BSNumber<UIntegerType> const& GetNumerator() const;
    inline BSNumber<UIntegerType> const& GetDenominator() const;

    // Comparisons.
    bool operator==(BSRational const& r) const;
    bool operator!=(BSRational const& r) const;
    bool operator< (BSRational const& r) const;
    bool operator<=(BSRational const& r) const;
    bool operator> (BSRational const& r) const;
    bool operator>=(BSRational const& r) const;

    // Unary operations.
    BSRational operator+() const;
    BSRational operator-() const;

    // Arithmetic.
    BSRational operator+(BSRational const& r) const;
    BSRational operator-(BSRational const& r) const;
    BSRational operator*(BSRational const& r) const;
    BSRational operator/(BSRational const& r) const;
    BSRational& operator+=(BSRational const& r);
    BSRational& operator-=(BSRational const& r);
    BSRational& operator*=(BSRational const& r);
    BSRational& operator/=(BSRational const& r);

    // Disk input/output.  The fstream objects should be created using
    // std::ios::binary.  The return value is 'true' iff the operation
    // was successful.
    bool Write(std::ofstream& output) const;
    bool Read(std::ifstream& input);

private:
    // Generic conversion code that converts to the correctly
    // rounded result using round-to-nearest-ties-to-even.
    template <typename UIntType, typename RealType>
    RealType Convert() const;

#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
public:
    // List this first so that it shows up first in the debugger watch window.
    double mValue;
private:
#endif

    BSNumber<UIntegerType> mNumerator, mDenominator;

    friend class UnitTestBSRational;
};


template <typename UIntegerType>
BSRational<UIntegerType>::BSRational()
    :
    mNumerator(0),
    mDenominator(1)
{
#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
    mValue = (double)*this;
#endif
}

template <typename UIntegerType>
BSRational<UIntegerType>::BSRational(BSRational const& r)
{
    *this = r;
}

template <typename UIntegerType>
BSRational<UIntegerType>::BSRational(float numerator)
    :
    mNumerator(numerator),
    mDenominator(1.0f)
{
#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
    mValue = (double)*this;
#endif
}

template <typename UIntegerType>
BSRational<UIntegerType>::BSRational(double numerator)
    :
    mNumerator(numerator),
    mDenominator(1.0)
{
#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
    mValue = (double)*this;
#endif
}

template <typename UIntegerType>
BSRational<UIntegerType>::BSRational(int32_t numerator)
    :
    mNumerator(numerator),
    mDenominator(1)
{
#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
    mValue = (double)*this;
#endif
}

template <typename UIntegerType>
BSRational<UIntegerType>::BSRational(uint32_t numerator)
    :
    mNumerator(numerator),
    mDenominator(1)
{
#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
    mValue = (double)*this;
#endif
}

template <typename UIntegerType>
BSRational<UIntegerType>::BSRational(int64_t numerator)
    :
    mNumerator(numerator),
    mDenominator(1)
{
#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
    mValue = (double)*this;
#endif
}

template <typename UIntegerType>
BSRational<UIntegerType>::BSRational(uint64_t numerator)
    :
    mNumerator(numerator),
    mDenominator(1)
{
#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
    mValue = (double)*this;
#endif
}

template <typename UIntegerType>
BSRational<UIntegerType>::BSRational(BSNumber<UIntegerType> const& numerator)
    :
    mNumerator(numerator),
    mDenominator(1)
{
#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
    mValue = (double)*this;
#endif
}

template <typename UIntegerType>
BSRational<UIntegerType>::BSRational(float numerator, float denominator)
    :
    mNumerator(numerator),
    mDenominator(denominator)
{
    LogAssert(mDenominator.mSign != 0, "Division by zero not allowed.");
    if (mDenominator.mSign < 0)
    {
        mNumerator.mSign = -mNumerator.mSign;
        mDenominator.mSign = 1;
    }
#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
    mValue = (double)*this;
#endif
}

template <typename UIntegerType>
BSRational<UIntegerType>::BSRational(double numerator, double denominator)
    :
    mNumerator(numerator),
    mDenominator(denominator)
{
    LogAssert(mDenominator.mSign != 0, "Division by zero not allowed.");
    if (mDenominator.mSign < 0)
    {
        mNumerator.mSign = -mNumerator.mSign;
        mDenominator.mSign = 1;
    }
#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
    mValue = (double)*this;
#endif
}

template <typename UIntegerType>
BSRational<UIntegerType>::BSRational(BSNumber<UIntegerType> const& numerator,
    BSNumber<UIntegerType> const& denominator)
    :
    mNumerator(numerator),
    mDenominator(denominator)
{
    LogAssert(mDenominator.mSign != 0, "Division by zero not allowed.");
    if (mDenominator.mSign < 0)
    {
        mNumerator.mSign = -mNumerator.mSign;
        mDenominator.mSign = 1;
    }

    // Set the exponent of the denominator to zero, but you can do so only
    // by modifying the biased exponent.  Adjust the numerator accordingly.
    // This prevents large growth of the exponents in both numerator and
    // denominator simultaneously.
    mNumerator.mBiasedExponent -= mDenominator.GetExponent();
    mDenominator.mBiasedExponent =
        -(mDenominator.GetUInteger().GetNumBits() - 1);

#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
    mValue = (double)*this;
#endif
}

template <typename UIntegerType>
BSRational<UIntegerType>::operator float() const
{
    return Convert<uint32_t, float>();
}

template <typename UIntegerType>
BSRational<UIntegerType>::operator double() const
{
    return Convert<uint64_t, double>();
}

template <typename UIntegerType>
BSRational<UIntegerType>& BSRational<UIntegerType>::operator=(
    BSRational const& r)
{
    mNumerator = r.mNumerator;
    mDenominator = r.mDenominator;
#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
    mValue = (double)*this;
#endif
    return *this;
}

template <typename UIntegerType>
BSRational<UIntegerType>::BSRational(BSRational&& r)
{
    *this = std::move(r);
}

template <typename UIntegerType>
BSRational<UIntegerType>& BSRational<UIntegerType>::operator=(BSRational&& r)
{
    mNumerator = std::move(r.mNumerator);
    mDenominator = std::move(r.mDenominator);
#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
    mValue = (double)*this;
#endif
    return *this;
}

template <typename UIntegerType> inline
int BSRational<UIntegerType>::GetSign() const
{
    return mNumerator.GetSign() * mDenominator.GetSign();
}

template <typename UIntegerType> inline
BSNumber<UIntegerType> const& BSRational<UIntegerType>::GetNumerator() const
{
    return mNumerator;
}

template <typename UIntegerType> inline
BSNumber<UIntegerType> const& BSRational<UIntegerType>::GetDenominator() const
{
    return mDenominator;
}

template <typename UIntegerType>
bool BSRational<UIntegerType>::operator==(BSRational const& r) const
{
    // Do inexpensive sign tests first for optimum performance.
    if (mNumerator.mSign != r.mNumerator.mSign)
    {
        return false;
    }
    if (mNumerator.mSign == 0)
    {
        // The numbers are both zero.
        return true;
    }

    return mNumerator * r.mDenominator == mDenominator * r.mNumerator;
}

template <typename UIntegerType>
bool BSRational<UIntegerType>::operator!=(BSRational const& r) const
{
    return !operator==(r);
}

template <typename UIntegerType>
bool BSRational<UIntegerType>::operator< (BSRational const& r) const
{
    // Do inexpensive sign tests first for optimum performance.
    if (mNumerator.mSign > 0)
    {
        if (r.mNumerator.mSign <= 0)
        {
            return false;
        }
    }
    else if (mNumerator.mSign == 0)
    {
        return r.mNumerator.mSign > 0;
    }
    else if (mNumerator.mSign < 0)
    {
        if (r.mNumerator.mSign >= 0)
        {
            return true;
        }
    }

    return mNumerator * r.mDenominator < mDenominator * r.mNumerator;
}

template <typename UIntegerType>
bool BSRational<UIntegerType>::operator<=(BSRational const& r) const
{
    return !operator>(r);
}

template <typename UIntegerType>
bool BSRational<UIntegerType>::operator> (BSRational const& r) const
{
    return r.operator<(*this);
}

template <typename UIntegerType>
bool BSRational<UIntegerType>::operator>=(BSRational const& r) const
{
    return !operator<(r);
}

template <typename UIntegerType>
BSRational<UIntegerType> BSRational<UIntegerType>::operator+() const
{
    return *this;
}

template <typename UIntegerType>
BSRational<UIntegerType> BSRational<UIntegerType>::operator-() const
{
    return BSRational(-mNumerator, mDenominator);
}

template <typename UIntegerType>
BSRational<UIntegerType> BSRational<UIntegerType>::operator+(
    BSRational const& r) const
{
    BSNumber<UIntegerType> product0 = mNumerator * r.mDenominator;
    BSNumber<UIntegerType> product1 = mDenominator * r.mNumerator;
    BSNumber<UIntegerType> numerator = product0 + product1;

    // Complex expressions can lead to 0/denom, where denom is not 1.
    if (numerator.mSign != 0)
    {
        BSNumber<UIntegerType> denominator = mDenominator * r.mDenominator;
        return BSRational(numerator, denominator);
    }
    else
    {
        return BSRational(0);
    }
}

template <typename UIntegerType>
BSRational<UIntegerType> BSRational<UIntegerType>::operator-(
    BSRational const& r) const
{
    BSNumber<UIntegerType> product0 = mNumerator * r.mDenominator;
    BSNumber<UIntegerType> product1 = mDenominator * r.mNumerator;
    BSNumber<UIntegerType> numerator = product0 - product1;

    // Complex expressions can lead to 0/denom, where denom is not 1.
    if (numerator.mSign != 0)
    {
        BSNumber<UIntegerType> denominator = mDenominator * r.mDenominator;
        return BSRational(numerator, denominator);
    }
    else
    {
        return BSRational(0);
    }
}

template <typename UIntegerType>
BSRational<UIntegerType> BSRational<UIntegerType>::operator*(
    BSRational const& r) const
{
    BSNumber<UIntegerType> numerator = mNumerator * r.mNumerator;

    // Complex expressions can lead to 0/denom, where denom is not 1.
    if (numerator.mSign != 0)
    {
        BSNumber<UIntegerType> denominator = mDenominator * r.mDenominator;
        return BSRational(numerator, denominator);
    }
    else
    {
        return BSRational(0);
    }
}

template <typename UIntegerType>
BSRational<UIntegerType> BSRational<UIntegerType>::operator/(
    BSRational const& r) const
{
    LogAssert(r.mNumerator.mSign != 0, "Division by zero not allowed.");

    BSNumber<UIntegerType> numerator = mNumerator * r.mDenominator;

    // Complex expressions can lead to 0/denom, where denom is not 1.
    if (numerator.mSign != 0)
    {
        BSNumber<UIntegerType> denominator = mDenominator * r.mNumerator;
        if (denominator.mSign < 0)
        {
            numerator.mSign = -numerator.mSign;
            denominator.mSign = 1;
        }
        return BSRational(numerator, denominator);
    }
    else
    {
        return BSRational(0);
    }
}

template <typename UIntegerType>
BSRational<UIntegerType>& BSRational<UIntegerType>::operator+=(
    BSRational const& r)
{
    *this = operator+(r);
    return *this;
}

template <typename UIntegerType>
BSRational<UIntegerType>& BSRational<UIntegerType>::operator-=(
    BSRational const& r)
{
    *this = operator-(r);
    return *this;
}

template <typename UIntegerType>
BSRational<UIntegerType>& BSRational<UIntegerType>::operator*=(
    BSRational const& r)
{
    *this = operator*(r);
    return *this;
}

template <typename UIntegerType>
BSRational<UIntegerType>& BSRational<UIntegerType>::operator/=(
    BSRational const& r)
{
    *this = operator/(r);
    return *this;
}

template <typename UIntegerType>
bool BSRational<UIntegerType>::Write(std::ofstream& output) const
{
    return mNumerator.Write(output) && mDenominator.Write(output);
}

template <typename UIntegerType>
bool BSRational<UIntegerType>::Read(std::ifstream& input)
{
    return mNumerator.Read(input) && mDenominator.Read(input);
}

template <typename UIntegerType>
template <typename UIntType, typename RealType>
RealType BSRational<UIntegerType>::Convert() const
{
    if (mNumerator.mSign == 0)
    {
        return (RealType)0;
    }

    // The ratio is abstractly of the form (1.u*2^p)/(1.v*2^q).  Convert to
    // the form (1.u/1.v)*2^{p-q}, if 1.u >= 1.v, or to the form
    // (2*(1.u)/1.v)*2*{p-q-1}) if 1.u < 1.v.  The final form n/d must be in
    // the interval [1,2).
    BSNumber<UIntegerType> n = mNumerator, d = mDenominator;
    int32_t sign = n.mSign * d.mSign;
    n.mSign = 1;
    d.mSign = 1;
    int32_t pmq = n.GetExponent() - d.GetExponent();
    n.mBiasedExponent = 1 - n.GetUInteger().GetNumBits();
    d.mBiasedExponent = 1 - d.GetUInteger().GetNumBits();
    if (BSNumber<UIntegerType>::LessThanIgnoreSign(n, d))
    {
        ++n.mBiasedExponent;
        --pmq;
    }

    // At this time, n/d = 1.c in [1,2).  Define the sequence of bits
    // w = 1c = w_{imax} w_{imax-1} ... w_0 w_{-1} w_{-2} ... where
    // imax = precision(RealType)-1 and w_{imax} = 1.

    // Compute 'precision' bits for w, the leading bit guaranteed to be 1
    // and occurring at index (1 << (precision-1)).
    BSNumber<UIntegerType> one(1), two(2);
    int const imax = std::numeric_limits<RealType>::digits - 1;
    UIntType w = 0;
    UIntType mask = ((UIntType)1 << imax);
    for (int i = imax; i >= 0; --i, mask >>= 1)
    {
        if (BSNumber<UIntegerType>::LessThanIgnoreSign(n, d))
        {
            n = two * n;
        }
        else
        {
            n = two * (n - d);
            w |= mask;
        }
    }

    // Apply the mode round-to-nearest-ties-to-even to decide whether to
    // round down or up.  We computed w = w_{imax} ... w_0.  The remainder
    // is n/d = w_{imax+1}.w_{imax+2}... in [0,2).  Compute n'/d = (n-d)/d
    // in [-1,1).  Round-to-nearest-ties-to-even mode is the following,
    // where we need only test the sign of n'.  A remainder of "half" is
    // the case n' = 0.
    //   Round down when n' < 0 or (n' = 0 and w_0 = 0):  use w
    //   Round up when n' > 0 or (n' = 0 and w_0 == 1):  use w+1
    n = n - d;
    if (n.mSign > 0 || (n.mSign == 0 && (w & 1) == 1))
    {
        ++w;
    }

    if (w > 0)
    {
        // Ensure that the low-order bit of w is 1 (required for BSNumber
        // integer part).
        int32_t trailing = GetTrailingBit(w);
        w >>= trailing;
        pmq += trailing;

        // Compute a BSNumber with integer part w and the appropriate
        // number of bits and exponents.
        BSNumber<UIntegerType> result(w);
        result.mBiasedExponent = pmq - imax;
        RealType converted = (RealType)result;
        if (sign < 0)
        {
            converted = -converted;
        }
        return converted;
    }
    else
    {
        return (RealType)0;
    }
}


}

namespace std
{
    template <typename UIntegerType> inline
    gte::BSRational<UIntegerType> abs(
        gte::BSRational<UIntegerType> const& number)
    {
        return (number.GetSign() >= 0 ? number : -number);
    }
}
