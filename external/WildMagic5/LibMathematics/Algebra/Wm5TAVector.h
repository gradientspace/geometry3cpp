// Geometric Tools, LLC
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.16.1 (2017/10/12)

#ifndef WM5TAVECTOR_H
#define WM5TAVECTOR_H

#include "Wm5MathematicsLIB.h"
#include "Wm5THPoint.h"
#include "Wm5Vector3.h"

namespace Wm5
{

template <typename Real>
class TAVector : public THPoint<Real>
{
public:
    // Construction and destruction.  TAVector represents an affine vector of
    // the form (x,y,z,0).  The destructor hides the HPoint destructor, which
    // is not a problem because there are no side effects that must occur in
    // the base class.
    TAVector()
    {
        mTuple[0] = (Real)0;
        mTuple[1] = (Real)0;
        mTuple[2] = (Real)0;
        mTuple[3] = (Real)0;
    }

    TAVector(const TAVector& vec)
    {
        mTuple[0] = vec.mTuple[0];
        mTuple[1] = vec.mTuple[1];
        mTuple[2] = vec.mTuple[2];
        mTuple[3] = (Real)0;
    }

    TAVector(Real x, Real y, Real z)
    {
        mTuple[0] = x;
        mTuple[1] = y;
        mTuple[2] = z;
        mTuple[3] = (Real)0;
    }

    TAVector(Vector3<Real> const& vec)
    {
        mTuple[0] = vec[0];
        mTuple[1] = vec[1];
        mTuple[2] = vec[2];
        mTuple[3] = (Real)0;
    }

    ~TAVector()
    {
    }

    // Implicit conversions.
    operator const Vector3<Real>& () const
    {
        return *(const Vector3<Real>*)mTuple;
    }

    operator Vector3<Real>& ()
    {
        return *(Vector3<Real>*)mTuple;
    }

    // Assignment.
    TAVector& operator= (const TAVector& vec)
    {
        mTuple[0] = vec.mTuple[0];
        mTuple[1] = vec.mTuple[1];
        mTuple[2] = vec.mTuple[2];
        mTuple[3] = (Real)0;
        return *this;
    }

    // Arithmetic operations.
    TAVector operator+ (const TAVector& vec) const
    {
        return TAVector
        (
            mTuple[0] + vec.mTuple[0],
            mTuple[1] + vec.mTuple[1],
            mTuple[2] + vec.mTuple[2]
        );
    }

    TAVector operator- (const TAVector& vec) const
    {
        return TAVector
        (
            mTuple[0] - vec.mTuple[0],
            mTuple[1] - vec.mTuple[1],
            mTuple[2] - vec.mTuple[2]
        );
    }

    TAVector operator* (Real scalar) const
    {
        return TAVector
        (
            scalar * mTuple[0],
            scalar * mTuple[1],
            scalar * mTuple[2]
        );
    }

    TAVector operator/ (Real scalar) const
    {
        if (scalar != (Real)0)
        {
            Real invScalar = (Real)1 / scalar;
            return TAVector
            (
                invScalar * mTuple[0],
                invScalar * mTuple[1],
                invScalar * mTuple[2]
            );
        }

        Real infinity = std::numeric_limits<Real>::infinity();
        return TAVector(infinity, infinity, infinity);
    }

    TAVector operator- () const
    {
        return TAVector(-mTuple[0], -mTuple[1], -mTuple[2]);
    }

    // Arithmetic updates.
    TAVector& operator+= (const TAVector& vec)
    {
        mTuple[0] += vec[0];
        mTuple[1] += vec[1];
        mTuple[2] += vec[2];
        return *this;
    }

    TAVector& operator-= (const TAVector& vec)
    {
        mTuple[0] -= vec[0];
        mTuple[1] -= vec[1];
        mTuple[2] -= vec[2];
        return *this;
    }

    TAVector& operator*= (Real scalar)
    {
        mTuple[0] *= scalar;
        mTuple[1] *= scalar;
        mTuple[2] *= scalar;
        return *this;
    }

    TAVector& operator/= (Real scalar)
    {
        mTuple[0] /= scalar;
        mTuple[1] /= scalar;
        mTuple[2] /= scalar;
        return *this;
    }

    // Vector operations.
    Real Length() const
    {
        Real sqrLength = mTuple[0] * mTuple[0] + mTuple[1] * mTuple[1] +
            mTuple[2] * mTuple[2];

        return sqrt(sqrLength);
    }

    Real SquaredLength() const
    {
        Real sqrLength = mTuple[0] * mTuple[0] + mTuple[1] * mTuple[1] +
            mTuple[2] * mTuple[2];

        return sqrLength;
    }

    Real Dot(const TAVector& vec) const
    {
        Real dotProduct = mTuple[0] * vec.mTuple[0] + mTuple[1] * vec.mTuple[1] +
            mTuple[2] * vec.mTuple[2];

        return dotProduct;
    }

    Real Normalize(const Real epsilon = (Real)0)
    {
        Real length = Length();

        if (length > epsilon)
        {
            mTuple[0] /= length;
            mTuple[1] /= length;
            mTuple[2] /= length;
        }
        else
        {
            length = (Real)0;
            mTuple[0] = (Real)0;
            mTuple[1] = (Real)0;
            mTuple[2] = (Real)0;
        }

        return length;
    }

    TAVector Cross(const TAVector& vec) const
    {
        return TAVector
        (
            mTuple[1] * vec.mTuple[2] - mTuple[2] * vec.mTuple[1],
            mTuple[2] * vec.mTuple[0] - mTuple[0] * vec.mTuple[2],
            mTuple[0] * vec.mTuple[1] - mTuple[1] * vec.mTuple[0]
        );
    }

    TAVector UnitCross(const TAVector& vec) const
    {
        TAVector cross
        (
            mTuple[1] * vec.mTuple[2] - mTuple[2] * vec.mTuple[1],
            mTuple[2] * vec.mTuple[0] - mTuple[0] * vec.mTuple[2],
            mTuple[0] * vec.mTuple[1] - mTuple[1] * vec.mTuple[0]
        );

        cross.Normalize();
        return cross;
    }

    // Inputs must be initialized nonzero vectors.
    static void Orthonormalize(TAVector& vec0, TAVector& vec1, TAVector& vec2)
    {
        // If the input vectors are v0, v1, and v2, then the Gram-Schmidt
        // orthonormalization produces vectors u0, u1, and u2 as follows,
        //
        //   u0 = v0/|v0|
        //   u1 = (v1-(u0*v1)u0)/|v1-(u0*v1)u0|
        //   u2 = (v2-(u0*v2)u0-(u1*v2)u1)/|v2-(u0*v2)u0-(u1*v2)u1|
        //
        // where |A| indicates length of vector A and A*B indicates dot
        // product of vectors A and B.

        // Compute u0.
        vec0.Normalize();

        // Compute u1.
        Real dot0 = vec0.Dot(vec1);
        vec1 -= dot0*vec0;
        vec1.Normalize();

        // Compute u2.
        Real dot1 = vec1.Dot(vec2);
        dot0 = vec0.Dot(vec2);
        vec2 -= dot0*vec0 + dot1*vec1;
        vec2.Normalize();
    }

    static void Orthonormalize(TAVector* vec)
    {
        Orthonormalize(vec[0], vec[1], vec[2]);
    }

    // Input vec2 must be a nonzero vector. The output is an orthonormal
    // basis {vec0,vec1,vec2}.  The input vec2 is normalized by this function.
    // If you know that vec2 is already unit length, use the function
    // GenerateComplementBasis to compute vec0 and vec1.
    static void GenerateOrthonormalBasis(TAVector& vec0, TAVector& vec1,
        TAVector& vec2)
    {
        vec2.Normalize();
        GenerateComplementBasis(vec0, vec1, vec2);
    }

    // Input vec0 must be a unit-length vector.  The output vectors
    // {vec0,vec1} are unit length and mutually perpendicular, and
    // {vec0,vec1,vec2} is an orthonormal basis.
    static void GenerateComplementBasis(TAVector& vec0, TAVector& vec1,
        const TAVector& vec2)
    {
        Real invLength;

        if (fabs(vec2.mTuple[0]) >= fabs(vec2.mTuple[1]))
        {
            // vec2.x or vec2.z is the largest magnitude component, swap them
            invLength = (Real)1 / sqrt(vec2.mTuple[0] * vec2.mTuple[0] +
                vec2.mTuple[2] * vec2.mTuple[2]);
            vec0.mTuple[0] = -vec2.mTuple[2] * invLength;
            vec0.mTuple[1] = (Real)0;
            vec0.mTuple[2] = +vec2.mTuple[0] * invLength;
            vec1.mTuple[0] = vec2.mTuple[1] * vec0.mTuple[2];
            vec1.mTuple[1] = vec2.mTuple[2] * vec0.mTuple[0] -
                vec2.mTuple[0] * vec0.mTuple[2];
            vec1.mTuple[2] = -vec2.mTuple[1] * vec0.mTuple[0];
        }
        else
        {
            // vec2.y or vec2.z is the largest magnitude component, swap them
            invLength = (Real)1 / sqrt(vec2.mTuple[1] * vec2.mTuple[1] +
                vec2.mTuple[2] * vec2.mTuple[2]);
            vec0.mTuple[0] = (Real)0;
            vec0.mTuple[1] = +vec2.mTuple[2] * invLength;
            vec0.mTuple[2] = -vec2.mTuple[1] * invLength;
            vec1.mTuple[0] = vec2.mTuple[1] * vec0.mTuple[2] -
                vec2.mTuple[2] * vec0.mTuple[1];
            vec1.mTuple[1] = -vec2.mTuple[0] * vec0.mTuple[2];
            vec1.mTuple[2] = vec2.mTuple[0] * vec0.mTuple[1];
        }
    }

    // Special vectors.
    static const TAVector ZERO()
    {
        return TAVector((Real)0, (Real)0, (Real)0);
    }

    static const TAVector UNIT_X()
    {
        return TAVector((Real)1, (Real)0, (Real)0);
    }

    static const TAVector UNIT_Y()
    {
        return TAVector((Real)0, (Real)1, (Real)0);
    }

    static const TAVector UNIT_Z()
    {
        return TAVector((Real)0, (Real)0, (Real)1);
    }

protected:
    using THPoint<Real>::mTuple;
};

template <typename Real>
TAVector<Real> operator* (Real scalar, const TAVector<Real>& vec)
{
    return vec * scalar;
}

}

#endif
