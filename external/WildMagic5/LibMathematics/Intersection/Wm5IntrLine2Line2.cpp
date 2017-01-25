// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.3 (2013/11/12)

#include "Wm5MathematicsPCH.h"
#include "Wm5IntrLine2Line2.h"

namespace Wm5
{
//----------------------------------------------------------------------------
template <typename Real>
IntrLine2Line2<Real>::IntrLine2Line2 (const Line2<Real>& line0,
    const Line2<Real>& line1)
    :
    mLine0(&line0),
    mLine1(&line1),
    mDotThreshold(Math<Real>::ZERO_TOLERANCE)
{
}
//----------------------------------------------------------------------------
template <typename Real>
const Line2<Real>& IntrLine2Line2<Real>::GetLine0 () const
{
    return *mLine0;
}
//----------------------------------------------------------------------------
template <typename Real>
const Line2<Real>& IntrLine2Line2<Real>::GetLine1 () const
{
    return *mLine1;
}
//----------------------------------------------------------------------------
template <typename Real>
bool IntrLine2Line2<Real>::Test ()
{
    mIntersectionType = Classify(mLine0->Origin, mLine0->Direction,
        mLine1->Origin, mLine1->Direction, mDotThreshold);

    if (mIntersectionType == IT_POINT)
    {
        mQuantity = 1;
    }
    else if (mIntersectionType == IT_LINE)
    {
        mQuantity = INT_MAX;
    }
    else
    {
        mQuantity = 0;
    }

    return mIntersectionType != IT_EMPTY;
}
//----------------------------------------------------------------------------
template <typename Real>
bool IntrLine2Line2<Real>::Find ()
{
    Real s[2];
    mIntersectionType = Classify(mLine0->Origin, mLine0->Direction,
        mLine1->Origin, mLine1->Direction, mDotThreshold, s);

    if (mIntersectionType == IT_POINT)
    {
        mQuantity = 1;
        mPoint = mLine0->Origin + s[0]*mLine0->Direction;
    }
    else if (mIntersectionType == IT_LINE)
    {
        mQuantity = INT_MAX;
    }
    else
    {
        mQuantity = 0;
    }

    return mIntersectionType != IT_EMPTY;
}
//----------------------------------------------------------------------------
template <typename Real>
int IntrLine2Line2<Real>::Classify (const Vector2<Real>& P0,
    const Vector2<Real>& D0, const Vector2<Real>& P1, const Vector2<Real>& D1,
    Real dotThreshold, Real* s)
{
    // Ensure dotThreshold is nonnegative.
    dotThreshold = std::max(dotThreshold, (Real)0);

    // The intersection of two lines is a solution to P0+s0*D0 = P1+s1*D1.
    // Rewrite this as s0*D0 - s1*D1 = P1 - P0 = Q.  If D0.Dot(Perp(D1)) = 0,
    // the lines are parallel.  Additionally, if Q.Dot(Perp(D1)) = 0, the
    // lines are the same.  If D0.Dot(Perp(D1)) is not zero, then
    //   s0 = Q.Dot(Perp(D1))/D0.Dot(Perp(D1))
    // produces the point of intersection.  Also,
    //   s1 = Q.Dot(Perp(D0))/D0.Dot(Perp(D1))

    Vector2<Real> diff = P1 - P0;
    Real D0DotPerpD1 = D0.DotPerp(D1);
    if (Math<Real>::FAbs(D0DotPerpD1) > dotThreshold)
    {
        // Lines intersect in a single point.
        if (s)
        {
            Real invD0DotPerpD1 = ((Real)1)/D0DotPerpD1;
            Real diffDotPerpD0 = diff.DotPerp(D0);
            Real diffDotPerpD1 = diff.DotPerp(D1);
            s[0] = diffDotPerpD1*invD0DotPerpD1;
            s[1] = diffDotPerpD0*invD0DotPerpD1;
        }
        return IT_POINT;
    }

    // Lines are parallel.
    diff.Normalize();
    Real diffNDotPerpD1 = diff.DotPerp(D1);
    if (Math<Real>::FAbs(diffNDotPerpD1) <= dotThreshold)
    {
        // Lines are colinear.
        return IT_LINE;
    }

    // Lines are parallel, but distinct.
    return IT_EMPTY;
}
//----------------------------------------------------------------------------
template <typename Real>
void IntrLine2Line2<Real>::SetDotThreshold (Real dotThreshold)
{
    if (dotThreshold >= (Real)0)
    {
        mDotThreshold = dotThreshold;
        return;
    }

    assertion(false, "Dot threshold must be nonnegative.");
}
//----------------------------------------------------------------------------
template <typename Real>
Real IntrLine2Line2<Real>::GetDotThreshold () const
{
    return mDotThreshold;
}
//----------------------------------------------------------------------------
template <typename Real>
int IntrLine2Line2<Real>::GetQuantity () const
{
    return mQuantity;
}
//----------------------------------------------------------------------------
template <typename Real>
const Vector2<Real>& IntrLine2Line2<Real>::GetPoint () const
{
    return mPoint;
}
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Explicit instantiation.
//----------------------------------------------------------------------------
template WM5_MATHEMATICS_ITEM
class IntrLine2Line2<float>;

template WM5_MATHEMATICS_ITEM
class IntrLine2Line2<double>;
//----------------------------------------------------------------------------
}
