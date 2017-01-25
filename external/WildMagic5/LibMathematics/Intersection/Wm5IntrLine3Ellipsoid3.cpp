// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.2 (2012/11/03)

#include "Wm5MathematicsPCH.h"
#include "Wm5IntrLine3Ellipsoid3.h"

namespace Wm5
{
//----------------------------------------------------------------------------
template <typename Real>
IntrLine3Ellipsoid3<Real>::IntrLine3Ellipsoid3 (const Line3<Real>& line,
    const Ellipsoid3<Real>& ellipsoid)
    :
    mLine(&line),
    mEllipsoid(&ellipsoid),
    mNegativeThreshold((Real)0),
    mPositiveThreshold((Real)0)
{
}
//----------------------------------------------------------------------------
template <typename Real>
const Line3<Real>& IntrLine3Ellipsoid3<Real>::GetLine () const
{
    return *mLine;
}
//----------------------------------------------------------------------------
template <typename Real>
const Ellipsoid3<Real>& IntrLine3Ellipsoid3<Real>::GetEllipsoid () const
{
    return *mEllipsoid;
}
//----------------------------------------------------------------------------
template <typename Real>
bool IntrLine3Ellipsoid3<Real>::Test ()
{
    // The ellipsoid is (X-K)^T*M*(X-K)-1 = 0 and the line is X = P+t*D.
    // Substitute the line equation into the ellipsoid equation to obtain
    // a quadratic equation
    //   Q(t) = a2*t^2 + 2*a1*t + a0 = 0
    // where a2 = D^T*M*D, a1 = D^T*M*(P-K), and a0 = (P-K)^T*M*(P-K)-1.

    Matrix3<Real> M;
    mEllipsoid->GetM(M);

    Vector3<Real> diff = mLine->Origin - mEllipsoid->Center;
    Vector3<Real> matDir = M*mLine->Direction;
    Vector3<Real> matDiff = M*diff;
    Real a2 = mLine->Direction.Dot(matDir);
    Real a1 = mLine->Direction.Dot(matDiff);
    Real a0 = diff.Dot(matDiff) - (Real)1;

    // Intersection occurs if Q(t) has real roots.
    Real discr = a1*a1 - a0*a2;
    return discr >= mNegativeThreshold;
}
//----------------------------------------------------------------------------
template <typename Real>
bool IntrLine3Ellipsoid3<Real>::Find ()
{
    // The ellipsoid is (X-K)^T*M*(X-K)-1 = 0 and the line is X = P+t*D.
    // Substitute the line equation into the ellipsoid equation to obtain
    // a quadratic equation
    //   Q(t) = a2*t^2 + 2*a1*t + a0 = 0
    // where a2 = D^T*M*D, a1 = D^T*M*(P-K), and a0 = (P-K)^T*M*(P-K)-1.

    Matrix3<Real> M;
    mEllipsoid->GetM(M);

    Vector3<Real> diff = mLine->Origin - mEllipsoid->Center;
    Vector3<Real> matDir = M*mLine->Direction;
    Vector3<Real> matDiff = M*diff;
    Real a2 = mLine->Direction.Dot(matDir);
    Real a1 = mLine->Direction.Dot(matDiff);
    Real a0 = diff.Dot(matDiff) - (Real)1;

    // Intersection occurs if Q(t) has real roots.
    Real discr = a1*a1 - a0*a2;
    Real t[2];
    if (discr < mNegativeThreshold)
    {
        mIntersectionType = IT_EMPTY;
        mQuantity = 0;
    }
    else if (discr > mPositiveThreshold)
    {
        mIntersectionType = IT_SEGMENT;
        mQuantity = 2;

        Real root = Math<Real>::Sqrt(discr);
        Real inv = ((Real)1)/a2;
        t[0] = (-a1 - root)*inv;
        t[1] = (-a1 + root)*inv;
        mPoint[0] = mLine->Origin + t[0]*mLine->Direction;
        mPoint[1] = mLine->Origin + t[1]*mLine->Direction;
    }
    else
    {
        mIntersectionType = IT_POINT;
        mQuantity = 1;

        t[0] = -a1/a2;
        mPoint[0] = mLine->Origin + t[0]*mLine->Direction;
    }

    return mIntersectionType != IT_EMPTY;
}
//----------------------------------------------------------------------------
template <typename Real>
int IntrLine3Ellipsoid3<Real>::GetQuantity () const
{
    return mQuantity;
}
//----------------------------------------------------------------------------
template <typename Real>
const Vector3<Real>& IntrLine3Ellipsoid3<Real>::GetPoint (int i) const
{
    return mPoint[i];
}
//----------------------------------------------------------------------------
template <typename Real>
void IntrLine3Ellipsoid3<Real>::SetNegativeThreshold (Real negThreshold)
{
    if (negThreshold <= (Real)0)
    {
        mNegativeThreshold = negThreshold;
        return;
    }

    assertion(false, "Negative threshold must be nonpositive.");
}
//----------------------------------------------------------------------------
template <typename Real>
Real IntrLine3Ellipsoid3<Real>::GetNegativeThreshold () const
{
    return mNegativeThreshold;
}
//----------------------------------------------------------------------------
template <typename Real>
void IntrLine3Ellipsoid3<Real>::SetPositiveThreshold (Real posThreshold)
{
    if (posThreshold >= (Real)0)
    {
        mPositiveThreshold = posThreshold;
        return;
    }

    assertion(false, "Positive threshold must be nonnegative.");
}
//----------------------------------------------------------------------------
template <typename Real>
Real IntrLine3Ellipsoid3<Real>::GetPositiveThreshold () const
{
    return mPositiveThreshold;
}
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Explicit instantiation.
//----------------------------------------------------------------------------
template WM5_MATHEMATICS_ITEM
class IntrLine3Ellipsoid3<float>;

template WM5_MATHEMATICS_ITEM
class IntrLine3Ellipsoid3<double>;
//----------------------------------------------------------------------------
}
