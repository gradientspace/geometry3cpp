// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <LowLevel/GteLogger.h>
#include <Mathematics/GteBitHacks.h>
#include <Mathematics/GteIEEEBinary.h>
#include <algorithm>
#include <fstream>
#include <limits>

// The class BSNumber (binary scientific number) is designed to provide exact
// arithmetic for robust algorithms, typically those for which we need to know
// the exact sign of determinants.  The template parameter UIntegerType must
// have support for at least the following public interface.  The fstream
// objects for Write/Read must be created using std::ios::binary.  The return
// value of Write/Read is 'true' iff the operation was successful.
//
//      class UIntegerType
//      {
//      public:
//          UIntegerType();
//          UIntegerType(UIntegerType const& number);
//          UIntegerType(uint32_t number);
//          UIntegerType(uint64_t number);
//          UIntegerType(int numBits);
//          UIntegerType& operator=(UIntegerType const& number);
//          UIntegerType(UIntegerType&& number);
//          UIntegerType& operator=(UIntegerType&& number);
//          int32_t GetNumBits() const;
//          bool operator==(UIntegerType const& number) const;
//          bool operator< (UIntegerType const& number) const;
//          void Add(UIntegerType const& n0, UIntegerType const& n1);
//          void Sub(UIntegerType const& n0, UIntegerType const& n1);
//          void Mul(UIntegerType const& n0, UIntegerType const& n1);
//          void ShiftLeft(UIntegerType const& number, int shift);
//          int32_t ShiftRightToOdd(UIntegerType const& number);
//          uint64_t GetPrefix(int numRequested) const;
//          bool Write(std::ofstream& output) const;
//          bool Read(std::ifstream& input);
//      };
//
// GTEngine currently has 32-bits-per-word storage for UIntegerType.  See the
// classes UIntegerAP32 (arbitrary precision), UIntegerFP32<N> (fixed
// precision), and UIntegerALU32 (arithmetic logic unit shared by the previous
// two classes).  The document at the following link describes the design,
// implementation, and use of BSNumber and BSRational.
//   http://www.geometrictools.com/Documentation/ArbitraryPrecision.pdf
//
// Support for debugging algorithms that use exact rational arithmetic.  Each
// BSNumber and BSRational has a double-precision member that is exposed when
// the conditional define is enabled.  Be aware that this can be very slow
// because of the conversion to double-precision whenever new objects are
// created by arithmetic operations.  As a faster alternative, you can add
// temporary code in your algorithms that explicitly convert specific rational
// numbers to double precision.
//
//#define GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE

namespace gte
{

template <typename UIntegerType> class BSRational;

template <typename UIntegerType>
class BSNumber
{
public:
    // Construction.  The default constructor generates the zero BSNumber.
    BSNumber();
    BSNumber(BSNumber const& number);
    BSNumber(float number);
    BSNumber(double number);
    BSNumber(int32_t number);
    BSNumber(uint32_t number);
    BSNumber(int64_t number);
    BSNumber(uint64_t number);

    // Implicit conversions.
    inline operator float() const;
    inline operator double() const;

    // Assignment.
    BSNumber& operator=(BSNumber const& number);

    // Support for std::move.
    BSNumber(BSNumber&& number);
    BSNumber& operator=(BSNumber&& number);

    // Member access.
    inline int32_t GetSign() const;
    inline int32_t GetBiasedExponent() const;
    inline int32_t GetExponent() const;
    inline UIntegerType const& GetUInteger() const;

    // Comparisons.
    bool operator==(BSNumber const& number) const;
    bool operator!=(BSNumber const& number) const;
    bool operator< (BSNumber const& number) const;
    bool operator<=(BSNumber const& number) const;
    bool operator> (BSNumber const& number) const;
    bool operator>=(BSNumber const& number) const;

    // Unary operations.
    BSNumber operator+() const;
    BSNumber operator-() const;

    // Arithmetic.
    BSNumber operator+(BSNumber const& number) const;
    BSNumber operator-(BSNumber const& number) const;
    BSNumber operator*(BSNumber const& number) const;
    BSNumber& operator+=(BSNumber const& number);
    BSNumber& operator-=(BSNumber const& number);
    BSNumber& operator*=(BSNumber const& number);

    // Disk input/output.  The fstream objects should be created using
    // std::ios::binary.  The return value is 'true' iff the operation
    // was successful.
    bool Write(std::ofstream& output) const;
    bool Read(std::ifstream& input);

private:
    // Helpers for operator==, operator<, operator+, operator-.
    static bool EqualIgnoreSign(BSNumber const& n0, BSNumber const& n1);
    static bool LessThanIgnoreSign(BSNumber const& n0, BSNumber const& n1);

    // Add two positive numbers.
    static BSNumber AddIgnoreSign(BSNumber const& n0, BSNumber const& n1,
        int32_t resultSign);

    // Subtract two positive numbers where n0 > n1.
    static BSNumber SubIgnoreSign(BSNumber const& n0, BSNumber const& n1,
        int32_t resultSign);

    // Support for conversions from floating-point numbers to BSNumber.
    template <typename IEEE>
    void ConvertFrom(typename IEEE::FloatType number);

    // Support for conversions from BSNumber to floating-point numbers.
    template <typename IEEE>
    typename IEEE::FloatType ConvertTo() const;

    template <typename IEEE>
    typename IEEE::UIntType GetTrailing(int32_t normal, int32_t sigma) const;

#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
public:
    // List this first so that it shows up first in the debugger watch window.
    double mValue;
private:
#endif

    // The number 0 is represented by: mSign = 0, mBiasedExponent = 0, and
    // mUInteger = 0.  For nonzero numbers, mSign != 0 and mUInteger > 0.
    int32_t mSign;
    int32_t mBiasedExponent;
    UIntegerType mUInteger;

    // Access to members to avoid exposing them publically when they are
    // needed only internally.
    friend class BSRational<UIntegerType>;
    friend class UnitTestBSNumber;
};


template <typename UIntegerType>
BSNumber<UIntegerType>::BSNumber()
    :
    mSign(0),
    mBiasedExponent(0)
{
#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
    mValue = (double)*this;
#endif
}

template <typename UIntegerType>
BSNumber<UIntegerType>::BSNumber(BSNumber const& number)
{
    *this = number;
}

template <typename UIntegerType>
BSNumber<UIntegerType>::BSNumber(float number)
{
    ConvertFrom<IEEEBinary32>(number);
#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
    mValue = (double)*this;
#endif
}

template <typename UIntegerType>
BSNumber<UIntegerType>::BSNumber(double number)
{
    ConvertFrom<IEEEBinary64>(number);
#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
    mValue = (double)*this;
#endif
}

template <typename UIntegerType>
BSNumber<UIntegerType>::BSNumber(int32_t number)
{
    if (number == 0)
    {
        mSign = 0;
        mBiasedExponent = 0;
    }
    else
    {
        if (number < 0)
        {
            mSign = -1;
            number = -number;
        }
        else
        {
            mSign = 1;
        }

        mBiasedExponent = GetTrailingBit(number);
        mUInteger = (uint32_t)number;
    }
#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
    mValue = (double)*this;
#endif
}

template <typename UIntegerType>
BSNumber<UIntegerType>::BSNumber(uint32_t number)
{
    if (number == 0)
    {
        mSign = 0;
        mBiasedExponent = 0;
    }
    else
    {
        mSign = 1;
        mBiasedExponent = GetTrailingBit(number);
        mUInteger = number;
    }
#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
    mValue = (double)*this;
#endif
}

template <typename UIntegerType>
BSNumber<UIntegerType>::BSNumber(int64_t number)
{
    if (number == 0)
    {
        mSign = 0;
        mBiasedExponent = 0;
    }
    else
    {
        if (number < 0)
        {
            mSign = -1;
            number = -number;
        }
        else
        {
            mSign = 1;
        }

        mBiasedExponent = GetTrailingBit(number);
        mUInteger = (uint64_t)number;
    }
#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
    mValue = (double)*this;
#endif
}

template <typename UIntegerType>
BSNumber<UIntegerType>::BSNumber(uint64_t number)
{
    if (number == 0)
    {
        mSign = 0;
        mBiasedExponent = 0;
    }
    else
    {
        mSign = 1;
        mBiasedExponent = GetTrailingBit(number);
        mUInteger = number;
    }
#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
    mValue = (double)*this;
#endif
}

template <typename UIntegerType> inline
BSNumber<UIntegerType>::operator float() const
{
    return ConvertTo<IEEEBinary32>();
}

template <typename UIntegerType>
inline BSNumber<UIntegerType>::operator double() const
{
    return ConvertTo<IEEEBinary64>();
}

template <typename UIntegerType>
BSNumber<UIntegerType>& BSNumber<UIntegerType>::operator=(
    BSNumber const& number)
{
    mSign = number.mSign;
    mBiasedExponent = number.mBiasedExponent;
    mUInteger = number.mUInteger;
#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
    mValue = number.mValue;
#endif
    return *this;
}

template <typename UIntegerType>
BSNumber<UIntegerType>::BSNumber(BSNumber&& number)
{
    *this = std::move(number);
}

template <typename UIntegerType>
BSNumber<UIntegerType>& BSNumber<UIntegerType>::operator=(BSNumber&& number)
{
    mSign = number.mSign;
    mBiasedExponent = number.mBiasedExponent;
    mUInteger = std::move(number.mUInteger);
    number.mSign = 0;
    number.mBiasedExponent = 0;
#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
    mValue = number.mValue;
#endif
    return *this;
}

template <typename UIntegerType> inline
int32_t BSNumber<UIntegerType>::GetSign() const
{
    return mSign;
}

template <typename UIntegerType> inline
int32_t BSNumber<UIntegerType>::GetBiasedExponent() const
{
    return mBiasedExponent;
}

template <typename UIntegerType> inline
int32_t BSNumber<UIntegerType>::GetExponent() const
{
    return mBiasedExponent + mUInteger.GetNumBits() - 1;
}

template <typename UIntegerType> inline
UIntegerType const& BSNumber<UIntegerType>::GetUInteger() const
{
    return mUInteger;
}

template <typename UIntegerType>
bool BSNumber<UIntegerType>::operator==(BSNumber const& number) const
{
    return (mSign == number.mSign ? EqualIgnoreSign(*this, number) : false);
}

template <typename UIntegerType>
bool BSNumber<UIntegerType>::operator!=(BSNumber const& number) const
{
    return !operator==(number);
}

template <typename UIntegerType>
bool BSNumber<UIntegerType>::operator<(BSNumber const& number) const
{
    if (mSign > 0)
    {
        if (number.mSign <= 0)
        {
            return false;
        }

        // Both numbers are positive.
        return LessThanIgnoreSign(*this, number);
    }
    else if (mSign < 0)
    {
        if (number.mSign >= 0)
        {
            return true;
        }

        // Both numbers are negative.
        return LessThanIgnoreSign(number, *this);
    }
    else
    {
        return number.mSign > 0;
    }
}

template <typename UIntegerType>
bool BSNumber<UIntegerType>::operator<=(BSNumber const& number) const
{
    return operator<(number) || operator==(number);
}

template <typename UIntegerType>
bool BSNumber<UIntegerType>::operator>(BSNumber const& number) const
{
    return !operator<=(number);
}

template <typename UIntegerType>
bool BSNumber<UIntegerType>::operator>=(BSNumber const& number) const
{
    return !operator<(number);
}

template <typename UIntegerType>
BSNumber<UIntegerType> BSNumber<UIntegerType>::operator+() const
{
    return *this;
}

template <typename UIntegerType>
BSNumber<UIntegerType> BSNumber<UIntegerType>::operator-() const
{
    BSNumber result = *this;
    result.mSign = -result.mSign;
#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
    result.mValue = (double)result;
#endif
    return result;
}

template <typename UIntegerType>
BSNumber<UIntegerType> BSNumber<UIntegerType>::operator+(BSNumber const& n1)
const
{
    BSNumber const& n0 = *this;

    if (n0.mSign == 0)
    {
        return n1;
    }

    if (n1.mSign == 0)
    {
        return n0;
    }

    if (n0.mSign > 0)
    {
        if (n1.mSign > 0)
        {
            // n0 + n1 = |n0| + |n1|
            return AddIgnoreSign(n0, n1, +1);
        }
        else // n1.mSign < 0
        {
            if (!EqualIgnoreSign(n0, n1))
            {
                if (LessThanIgnoreSign(n1, n0))
                {
                    // n0 + n1 = |n0| - |n1| > 0
                    return SubIgnoreSign(n0, n1, +1);
                }
                else
                {
                    // n0 + n1 = -(|n1| - |n0|) < 0
                    return SubIgnoreSign(n1, n0, -1);
                }
            }
            // else n0 + n1 = 0
        }
    }
    else // n0.mSign < 0
    {
        if (n1.mSign < 0)
        {
            // n0 + n1 = -(|n0| + |n1|)
            return AddIgnoreSign(n0, n1, -1);
        }
        else // n1.mSign > 0
        {
            if (!EqualIgnoreSign(n0, n1))
            {
                if (LessThanIgnoreSign(n1, n0))
                {
                    // n0 + n1 = -(|n0| - |n1|) < 0
                    return SubIgnoreSign(n0, n1, -1);
                }
                else
                {
                    // n0 + n1 = |n1| - |n0| > 0
                    return SubIgnoreSign(n1, n0, +1);
                }
            }
            // else n0 + n1 = 0
        }
    }

    return BSNumber();  // = 0
}

template <typename UIntegerType>
BSNumber<UIntegerType> BSNumber<UIntegerType>::operator-(BSNumber const& n1)
const
{
    BSNumber const& n0 = *this;

    if (n0.mSign == 0)
    {
        return -n1;
    }

    if (n1.mSign == 0)
    {
        return n0;
    }

    if (n0.mSign > 0)
    {
        if (n1.mSign < 0)
        {
            // n0 - n1 = |n0| + |n1|
            return AddIgnoreSign(n0, n1, +1);
        }
        else // n1.mSign > 0
        {
            if (!EqualIgnoreSign(n0, n1))
            {
                if (LessThanIgnoreSign(n1, n0))
                {
                    // n0 - n1 = |n0| - |n1| > 0
                    return SubIgnoreSign(n0, n1, +1);
                }
                else
                {
                    // n0 - n1 = -(|n1| - |n0|) < 0
                    return SubIgnoreSign(n1, n0, -1);
                }
            }
            // else n0 - n1 = 0
        }
    }
    else // n0.mSign < 0
    {
        if (n1.mSign > 0)
        {
            // n0 - n1 = -(|n0| + |n1|)
            return AddIgnoreSign(n0, n1, -1);
        }
        else // n1.mSign < 0
        {
            if (!EqualIgnoreSign(n0, n1))
            {
                if (LessThanIgnoreSign(n1, n0))
                {
                    // n0 - n1 = -(|n0| - |n1|) < 0
                    return SubIgnoreSign(n0, n1, -1);
                }
                else
                {
                    // n0 - n1 = |n1| - |n0| > 0
                    return SubIgnoreSign(n1, n0, +1);
                }
            }
            // else n0 - n1 = 0
        }
    }

    return BSNumber();  // = 0
}

template <typename UIntegerType>
BSNumber<UIntegerType> BSNumber<UIntegerType>::operator*(
    BSNumber const& number) const
{
    BSNumber result;  // = 0
    int sign = mSign * number.mSign;
    if (sign != 0)
    {
        result.mSign = sign;
        result.mBiasedExponent = mBiasedExponent + number.mBiasedExponent;
        result.mUInteger.Mul(mUInteger, number.mUInteger);
    }
#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
    result.mValue = (double)result;
#endif
    return result;
}

template <typename UIntegerType>
BSNumber<UIntegerType>& BSNumber<UIntegerType>::operator+=(
    BSNumber const& number)
{
    *this = operator+(number);
    return *this;
}

template <typename UIntegerType>
BSNumber<UIntegerType>& BSNumber<UIntegerType>::operator-=(
    BSNumber const& number)
{
    *this = operator-(number);
    return *this;
}

template <typename UIntegerType>
BSNumber<UIntegerType>& BSNumber<UIntegerType>::operator*= (
    BSNumber const& number)
{
    *this = operator*(number);
    return *this;
}

template <typename UIntegerType>
bool BSNumber<UIntegerType>::Write(std::ofstream& output) const
{
    if (output.write((char const*)&mSign, sizeof(mSign)).bad())
    {
        return false;
    }

    if (output.write((char const*)&mBiasedExponent,
        sizeof(mBiasedExponent)).bad())
    {
        return false;
    }

    return mUInteger.Write(output);
}

template <typename UIntegerType>
bool BSNumber<UIntegerType>::Read(std::ifstream& input)
{
    if (input.read((char*)&mSign, sizeof(mSign)).bad())
    {
        return false;
    }

    if (input.read((char*)&mBiasedExponent, sizeof(mBiasedExponent)).bad())
    {
        return false;
    }

    return mUInteger.Read(input);
}

template <typename UIntegerType>
bool BSNumber<UIntegerType>::EqualIgnoreSign(BSNumber const& n0,
    BSNumber const& n1)
{
    return n0.mBiasedExponent == n1.mBiasedExponent
        && n0.mUInteger == n1.mUInteger;
}

template <typename UIntegerType>
bool BSNumber<UIntegerType>::LessThanIgnoreSign(BSNumber const& n0,
    BSNumber const& n1)
{
    int32_t e0 = n0.GetExponent(), e1 = n1.GetExponent();
    if (e0 < e1)
    {
        return true;
    }
    if (e0 > e1)
    {
        return false;
    }
    return n0.mUInteger < n1.mUInteger;
}

template <typename UIntegerType>
BSNumber<UIntegerType> BSNumber<UIntegerType>::AddIgnoreSign(
    BSNumber const& n0, BSNumber const& n1, int32_t resultSign)
{
    BSNumber result, temp;

    int32_t diff = n0.mBiasedExponent - n1.mBiasedExponent;
    if (diff > 0)
    {
        temp.mUInteger.ShiftLeft(n0.mUInteger, diff);
        result.mUInteger.Add(temp.mUInteger, n1.mUInteger);
        result.mBiasedExponent = n1.mBiasedExponent;
    }
    else if (diff < 0)
    {
        temp.mUInteger.ShiftLeft(n1.mUInteger, -diff);
        result.mUInteger.Add(n0.mUInteger, temp.mUInteger);
        result.mBiasedExponent = n0.mBiasedExponent;
    }
    else
    {
        temp.mUInteger.Add(n0.mUInteger, n1.mUInteger);
        int32_t shift = result.mUInteger.ShiftRightToOdd(temp.mUInteger);
        result.mBiasedExponent = n0.mBiasedExponent + shift;
    }

    result.mSign = resultSign;
#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
    result.mValue = (double)result;
#endif
    return result;
}

template <typename UIntegerType>
BSNumber<UIntegerType> BSNumber<UIntegerType>::SubIgnoreSign(
    BSNumber const& n0, BSNumber const& n1, int32_t resultSign)
{
    BSNumber result, temp;

    int32_t diff = n0.mBiasedExponent - n1.mBiasedExponent;
    if (diff > 0)
    {
        temp.mUInteger.ShiftLeft(n0.mUInteger, diff);
        result.mUInteger.Sub(temp.mUInteger, n1.mUInteger);
        result.mBiasedExponent = n1.mBiasedExponent;
    }
    else if (diff < 0)
    {
        temp.mUInteger.ShiftLeft(n1.mUInteger, -diff);
        result.mUInteger.Sub(n0.mUInteger, temp.mUInteger);
        result.mBiasedExponent = n0.mBiasedExponent;
    }
    else
    {
        temp.mUInteger.Sub(n0.mUInteger, n1.mUInteger);
        int32_t shift = result.mUInteger.ShiftRightToOdd(temp.mUInteger);
        result.mBiasedExponent = n0.mBiasedExponent + shift;
    }

    result.mSign = resultSign;
#if defined(GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE)
    result.mValue = (double)result;
#endif
    return result;
}

template <typename UIntegerType>
template <typename IEEE>
void BSNumber<UIntegerType>::ConvertFrom(typename IEEE::FloatType number)
{
    IEEE x(number);

    // Extract sign s, biased exponent e, and trailing significand t.
    typename IEEE::UIntType s = x.GetSign();
    typename IEEE::UIntType e = x.GetBiased();
    typename IEEE::UIntType t = x.GetTrailing();

    if (e == 0)
    {
        if (t == 0)  // zeros
        {
            // x = (-1)^s * 0
            mSign = 0;
            mBiasedExponent = 0;
        }
        else  // subnormal numbers
        {
            // x = (-1)^s * 0.t * 2^{1-EXPONENT_BIAS}
            int32_t last = GetTrailingBit(t);
            int32_t diff = IEEE::NUM_TRAILING_BITS - last;
            mSign = (s > 0 ? -1 : 1);
            mBiasedExponent = IEEE::MIN_SUB_EXPONENT - diff;
            mUInteger = (t >> last);
        }
    }
    else if (e < IEEE::MAX_BIASED_EXPONENT)  // normal numbers
    {
        // x = (-1)^s * 1.t * 2^{e-EXPONENT_BIAS}
        if (t > 0)
        {
            int32_t last = GetTrailingBit(t);
            int32_t diff = IEEE::NUM_TRAILING_BITS - last;
            mSign = (s > 0 ? -1 : 1);
            mBiasedExponent =
                static_cast<int32_t>(e)-IEEE::EXPONENT_BIAS - diff;
            mUInteger = ((t | IEEE::SUP_TRAILING) >> last);
        }
        else
        {
            mSign = (s > 0 ? -1 : 1);
            mBiasedExponent = static_cast<int32_t>(e)-IEEE::EXPONENT_BIAS;
            mUInteger = (typename IEEE::UIntType)1;
        }
    }
    else  // e == MAX_BIASED_EXPONENT, special numbers
    {
        if (t == 0)  // infinities
        {
            // x = (-1)^s * infinity
            LogWarning("Input is " + std::string(s > 0 ? "-" : "+") +
                "infinity.");

            // Return (-1)^s * 2^{1+EXPONENT_BIAS} for a graceful exit.
            mSign = (s > 0 ? -1 : 1);
            mBiasedExponent = 1 + IEEE::EXPONENT_BIAS;
            mUInteger = (typename IEEE::UIntType)1;
        }
        else  // not-a-number (NaN)
        {
            LogError("Input is a " +
                std::string(t & IEEE::NAN_QUIET_MASK ?
                "quiet" : "signaling") + " NaN with payload " +
                std::to_string(t & IEEE::NAN_PAYLOAD_MASK) + ".");

            // Return 0 for a graceful exit.
            mSign = 0;
            mBiasedExponent = 0;
        }
    }
}

template <typename UIntegerType>
template <typename IEEE>
typename IEEE::FloatType BSNumber<UIntegerType>::ConvertTo() const
{
    typename IEEE::UIntType s = (mSign < 0 ? 1 : 0);
    typename IEEE::UIntType e, t;

    if (mSign != 0)
    {
        // The conversions use round-to-nearest-ties-to-even semantics.
        int32_t exponent = GetExponent();
        if (exponent < IEEE::MIN_EXPONENT)
        {
            if (exponent < IEEE::MIN_EXPONENT - 1
                || mUInteger.GetNumBits() == 1)  // x = 1.0*2^{MIN_EXPONENT-1}
            {
                // Round to zero.
                e = 0;
                t = 0;
            }
            else
            {
                // Round to min subnormal.
                e = 0;
                t = 1;
            }
        }
        else if (exponent < IEEE::MIN_SUB_EXPONENT)
        {
            // The second input is in {0, ..., NUM_TRAILING_BITS-1}.
            t = GetTrailing<IEEE>(0, IEEE::MIN_SUB_EXPONENT - exponent - 1);
            if (t & IEEE::SUP_TRAILING)
            {
                // Leading NUM_SIGNIFICAND_BITS bits were all 1, so round to
                // min normal.
                e = 1;
                t = 0;
            }
            else
            {
                e = 0;
            }
        }
        else if (exponent <= IEEE::EXPONENT_BIAS)
        {
            e = static_cast<uint32_t>(exponent + IEEE::EXPONENT_BIAS);
            t = GetTrailing<IEEE>(1, 0);
            if (t & (IEEE::SUP_TRAILING << 1))
            {
                // Carry-out occurred, so increase exponent by 1 and
                // shift right to compensate.
                ++e;
                t >>= 1;
            }
            // Eliminate the leading 1 (implied for normals).
            t &= ~IEEE::SUP_TRAILING;
        }
        else
        {
            // Set to infinity.
            e = IEEE::MAX_BIASED_EXPONENT;
            t = 0;
        }
    }
    else
    {
        // The input is zero.
        e = 0;
        t = 0;
    }

    IEEE x(s, e, t);
    return x.number;
}

template <typename UIntegerType>
template <typename IEEE>
typename IEEE::UIntType BSNumber<UIntegerType>::GetTrailing(int32_t normal,
    int32_t sigma) const
{
    int32_t const numRequested = IEEE::NUM_SIGNIFICAND_BITS + normal;

    // We need numRequested bits to determine rounding direction.  These are
    // stored in the high-order bits of 'prefix'.
    uint64_t prefix = mUInteger.GetPrefix(numRequested);

    // The first bit index after the implied binary point for rounding.
    int32_t diff = numRequested - sigma;
    int32_t roundBitIndex = 64 - diff;

    // Determine rounding value based on round-to-nearest-ties-to-even.
    uint64_t mask = (GTE_U64(1) << roundBitIndex);
    uint64_t round;
    if (prefix & mask)
    {
        // The first bit of the remainder is 1.
        if (mUInteger.GetNumBits() == diff)
        {
            // The first bit of the remainder is the lowest-order bit of
            // mBits[0].  Apply the ties-to-even rule.
            if (prefix & (mask << 1))
            {
                // The last bit of the trailing significand is odd, so
                // round up.
                round = 1;
            }
            else
            {
                // The last bit of the trailing significand is even, so
                // round down.
                round = 0;
            }
        }
        else
        {
            // The first bit of the remainder is not the lowest-order bit of
            // mBits[0].  The remainder as a fraction is larger than 1/2, so
            // round up.
            round = 1;
        }
    }
    else
    {
        // The first bit of the remainder is 0, so round down.
        round = 0;
    }

    // Get the unrounded trailing significand.
    uint64_t trailing = prefix >> (roundBitIndex + 1);

    // Apply the rounding.
    trailing += round;
    return static_cast<typename IEEE::UIntType>(trailing);
}


}

namespace std
{
    template <typename UIntegerType> inline
    gte::BSNumber<UIntegerType> abs(gte::BSNumber<UIntegerType> const& number)
    {
        return (number.GetSign() >= 0 ? number : -number);
    }
}
