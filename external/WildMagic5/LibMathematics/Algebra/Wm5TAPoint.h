// Geometric Tools, LLC
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.16.1 (2017/10/12)

#ifndef WM5TAPOINT_H
#define WM5TAPOINT_H

#include "Wm5MathematicsLIB.h"
#include "Wm5TAVector.h"
#include "Wm5Vector3.h"

namespace Wm5
{

template <typename Real>
class TAPoint : public THPoint<Real>
{
public:
    // Construction and destruction.  TAPoint represents an affine point of the
    // form (x,y,z,1).  The destructor hides the HPoint destructor, which is
    // not a problem because there are no side effects that must occur in the
    // base class.
    TAPoint()
    {
        mTuple[0] = (Real)0;
        mTuple[1] = (Real)0;
        mTuple[2] = (Real)0;
        mTuple[3] = (Real)1;
    }

    TAPoint(const TAPoint& pnt)
    {
        mTuple[0] = pnt.mTuple[0];
        mTuple[1] = pnt.mTuple[1];
        mTuple[2] = pnt.mTuple[2];
        mTuple[3] = (Real)1;
    }

    TAPoint(Real x, Real y, Real z)
    {
        mTuple[0] = x;
        mTuple[1] = y;
        mTuple[2] = z;
        mTuple[3] = (Real)1;
    }

    TAPoint(Vector3<Real> const& pnt)
    {
        mTuple[0] = pnt[0];
        mTuple[1] = pnt[1];
        mTuple[2] = pnt[2];
        mTuple[3] = (Real)1;
    }

    ~TAPoint()
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
    TAPoint& operator= (const TAPoint& pnt)
    {
        mTuple[0] = pnt.mTuple[0];
        mTuple[1] = pnt.mTuple[1];
        mTuple[2] = pnt.mTuple[2];
        mTuple[3] = (Real)1;
        return *this;
    }

    // Arithmetic operations supported by affine algebra.

    // A point minus a point is a vector.
    TAVector<Real> operator- (const TAPoint& pnt) const
    {
        return TAVector<Real>
        (
            mTuple[0] - pnt.mTuple[0],
            mTuple[1] - pnt.mTuple[1],
            mTuple[2] - pnt.mTuple[2]
        );
    }

    // A point plus or minus a vector is a point.
    TAPoint operator+ (const TAVector<Real>& vec) const
    {
        return TAPoint
        (
            mTuple[0] + vec[0],
            mTuple[1] + vec[1],
            mTuple[2] + vec[2]
        );
    }

    TAPoint operator- (const TAVector<Real>& vec) const
    {
        return TAPoint
        (
            mTuple[0] - vec[0],
            mTuple[1] - vec[1],
            mTuple[2] - vec[2]
        );
    }

    TAPoint& operator+= (const TAVector<Real>& vec)
    {
        mTuple[0] += vec[0];
        mTuple[1] += vec[1];
        mTuple[2] += vec[2];
        return *this;
    }

    TAPoint& operator-= (const TAVector<Real>& vec)
    {
        mTuple[0] -= vec[0];
        mTuple[1] -= vec[1];
        mTuple[2] -= vec[2];
        return *this;
    }

    // In affine algebra, points cannot be added arbitrarily.  However,
    // affine sums and affine differences are allowed.  You are responsible
    // for ensuring that you are computing one or the other.
    //
    // An affine sum is of the form
    //   r = s1*p1 + s2*p2 + ... + sn*pn
    // where p1 through pn are homogeneous points (w-values are 1) and
    // s1 through sn are scalars for which s1 + s2 + ... + sn = 1.  The
    // result r is a homogenous point.
    //
    // An affine difference is of the form
    //   r = d1*p1 + d2*p2 + ... + dn*pn
    // where p1 through pn are homogeneous points (w-values are 1) and
    // d1 through dn are scalars for which d1 + d2 + ... + dn = 0.  The
    // result r is a homogeneous vector.  NOTE:  The arithemtic operations
    // of this class return TAPoint objects, but the affine difference needs
    // to be a THPoint object.  The following code shows how to accomplish
    // this:
    //   TAPoint p1, p2, p3;  // initialized to whatever
    //   THPoint r = 1.5*p1 + (-0.2)*p2 + (-0.3)*p3;
    // The right-hand side is computed using TAPoint operations, so it is an
    // TAPoint object.  THPoint has a constructor that takes a 'const Real*'.
    // TAPoint has an implicit conversion to 'const Real*'.  The code
    //   TAPoint somePoint;  // initialized to whatever
    //   THPoint r = somePoint;
    // involves an TAPoint implicit conversion to 'const Real*' followed by
    // an THPoint(const Real*) constructor call.  The latter copies only the
    // x, y, and z components and sets the w component to zero.
    TAPoint operator+ (const TAPoint& pnt) const
    {
        return TAPoint
        (
            mTuple[0] + pnt.mTuple[0],
            mTuple[1] + pnt.mTuple[1],
            mTuple[2] + pnt.mTuple[2]
        );
    }

    TAPoint operator* (Real scalar) const
    {
        return TAPoint
        (
            scalar * mTuple[0],
            scalar * mTuple[1],
            scalar * mTuple[2]
        );
    }

    TAPoint operator/ (Real scalar) const
    {
        return TAPoint
        (
            mTuple[0] / scalar,
            mTuple[1] / scalar,
            mTuple[2] / scalar
        );
    }

    TAPoint& operator+= (const TAPoint& pnt)
    {
        mTuple[0] += pnt[0];
        mTuple[1] += pnt[1];
        mTuple[2] += pnt[2];
        return *this;
    }

    TAPoint& operator-= (const TAPoint& pnt)
    {
        mTuple[0] -= pnt[0];
        mTuple[1] -= pnt[1];
        mTuple[2] -= pnt[2];
        return *this;
    }

    TAPoint& operator*= (Real scalar)
    {
        mTuple[0] *= scalar;
        mTuple[1] *= scalar;
        mTuple[2] *= scalar;
        return *this;
    }

    TAPoint& operator/= (Real scalar)
    {
        mTuple[0] /= scalar;
        mTuple[1] /= scalar;
        mTuple[2] /= scalar;
        return *this;
    }

    TAPoint operator- () const
    {
        return TAPoint(-mTuple[0], -mTuple[1], -mTuple[2]);
    }

    // The dot product between a point and a vector is not allowed in affine
    // algebra.  However, it is convenient to have one defined when dealing
    // with planes.  Specifically, a plane is Dot(N,X-P) = 0, where N is a
    // vector, P is a specific point on the plane, and X is any point on the
    // plane.  The difference X-P is a vector, so Dot(N,X-P) is well defined.
    // If the plane is rewritten as Dot(N,X) = Dot(N,P), this is not supported
    // by affine algebra, but we sometimes need to compute Dot(N,P) anyway.
    // In the following, the w-component of the TAPoint is ignored, which means
    // the TAPoint is treated as a vector.
    Real Dot(const TAVector<Real>& vec) const
    {
        return mTuple[0] * vec[0] + mTuple[1] * vec[1] + mTuple[2] * vec[2];
    }

    // Special vector.
    static const TAPoint ORIGIN() { return TAPoint((Real)0, (Real)0, (Real)0); }

protected:
    using THPoint<Real>::mTuple;
};

template <typename Real>
TAPoint<Real> operator* (Real scalar, const TAPoint<Real>& pnt)
{
    return pnt*scalar;
}

}

#endif
