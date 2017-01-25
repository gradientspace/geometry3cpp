// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.3 (2013/11/12)

#include "Wm5MathematicsPCH.h"
#include "Wm5IntrSegment2Segment2.h"
#include "Wm5IntrLine2Line2.h"
#include "Wm5Intersector1.h"

namespace Wm5
{
//----------------------------------------------------------------------------
template <typename Real>
IntrSegment2Segment2<Real>::IntrSegment2Segment2 (
    const Segment2<Real>& segment0, const Segment2<Real>& segment1)
    :
    mSegment0(&segment0),
    mSegment1(&segment1),
    mIntervalThreshold((Real)0),
    mDotThreshold(Math<Real>::ZERO_TOLERANCE)
{
}
//----------------------------------------------------------------------------
template <typename Real>
const Segment2<Real>& IntrSegment2Segment2<Real>::GetSegment0 () const
{
    return *mSegment0;
}
//----------------------------------------------------------------------------
template <typename Real>
const Segment2<Real>& IntrSegment2Segment2<Real>::GetSegment1 () const
{
    return *mSegment1;
}
//----------------------------------------------------------------------------
template <typename Real>
bool IntrSegment2Segment2<Real>::Test ()
{
    Real s[2];
    mIntersectionType = IntrLine2Line2<Real>::Classify(mSegment0->Center,
        mSegment0->Direction, mSegment1->Center, mSegment1->Direction,
        mDotThreshold, s);

    if (mIntersectionType == IT_POINT)
    {
        // Test whether the line-line intersection is on the segments.
        if (Math<Real>::FAbs(s[0]) <= mSegment0->Extent + mIntervalThreshold
        &&  Math<Real>::FAbs(s[1]) <= mSegment1->Extent + mIntervalThreshold)
        {
            mQuantity = 1;
        }
        else
        {
            mQuantity = 0;
            mIntersectionType = IT_EMPTY;
        }
    }
    else if (mIntersectionType == IT_LINE)
    {
        // Compute the location of segment1 endpoints relative to segment0.
        Vector2<Real> diff = mSegment1->Center - mSegment0->Center;
        Real t1 = mSegment0->Direction.Dot(diff);
        Real tmin = t1 - mSegment1->Extent;
        Real tmax = t1 + mSegment1->Extent;
        Intersector1<Real> calc(-mSegment0->Extent, mSegment0->Extent,
            tmin, tmax);
        calc.Find();  // This is intentionally not calc.Test(), need mQuantity.
        mQuantity = calc.GetNumIntersections();
        if (mQuantity == 2)
        {
            mIntersectionType = IT_SEGMENT;
        }
        else if (mQuantity == 1)
        {
            mIntersectionType = IT_POINT;
        }
        else
        {
            mIntersectionType = IT_EMPTY;
        }
    }
    else
    {
        mQuantity = 0;
    }

    return mIntersectionType != IT_EMPTY;
}
//----------------------------------------------------------------------------
template <typename Real>
bool IntrSegment2Segment2<Real>::Find ()
{
    Real s[2];
    mIntersectionType = IntrLine2Line2<Real>::Classify(mSegment0->Center,
        mSegment0->Direction, mSegment1->Center, mSegment1->Direction,
        mDotThreshold, s);

    if (mIntersectionType == IT_POINT)
    {
        // Test whether the line-line intersection is on the segments.
        if (Math<Real>::FAbs(s[0]) <= mSegment0->Extent + mIntervalThreshold
        &&  Math<Real>::FAbs(s[1]) <= mSegment1->Extent + mIntervalThreshold)
        {
            mQuantity = 1;
            mPoint[0] = mSegment0->Center + s[0]*mSegment0->Direction;
        }
        else
        {
            mQuantity = 0;
            mIntersectionType = IT_EMPTY;
        }
    }
    else if (mIntersectionType == IT_LINE)
    {
        // Compute the location of segment1 endpoints relative to segment0.
        Vector2<Real> diff = mSegment1->Center - mSegment0->Center;
        Real t1 = mSegment0->Direction.Dot(diff);
        Real tmin = t1 - mSegment1->Extent;
        Real tmax = t1 + mSegment1->Extent;
        Intersector1<Real> calc(-mSegment0->Extent, mSegment0->Extent,
            tmin, tmax);
        calc.Find();
        mQuantity = calc.GetNumIntersections();
        if (mQuantity == 2)
        {
            mIntersectionType = IT_SEGMENT;
            mPoint[0] = mSegment0->Center +
                calc.GetIntersection(0)*mSegment0->Direction;
            mPoint[1] = mSegment0->Center +
                calc.GetIntersection(1)*mSegment0->Direction;
        }
        else if (mQuantity == 1)
        {
            mIntersectionType = IT_POINT;
            mPoint[0] = mSegment0->Center +
                calc.GetIntersection(0)*mSegment0->Direction;
        }
        else
        {
            mIntersectionType = IT_EMPTY;
        }
    }
    else
    {
        mQuantity = 0;
    }

    return mIntersectionType != IT_EMPTY;
}
//----------------------------------------------------------------------------
template <typename Real>
void IntrSegment2Segment2<Real>::SetIntervalThreshold (Real intervalThreshold)
{
    if (intervalThreshold >= (Real)0)
    {
        mIntervalThreshold = intervalThreshold;
        return;
    }

    assertion(false, "Interval threshold must be nonnegative.");
}
//----------------------------------------------------------------------------
template <typename Real>
Real IntrSegment2Segment2<Real>::GetIntervalThreshold () const
{
    return mIntervalThreshold;
}
//----------------------------------------------------------------------------
template <typename Real>
void IntrSegment2Segment2<Real>::SetDotThreshold (Real dotThreshold)
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
Real IntrSegment2Segment2<Real>::GetDotThreshold () const
{
    return mDotThreshold;
}
//----------------------------------------------------------------------------
template <typename Real>
int IntrSegment2Segment2<Real>::GetQuantity () const
{
    return mQuantity;
}
//----------------------------------------------------------------------------
template <typename Real>
const Vector2<Real>& IntrSegment2Segment2<Real>::GetPoint (int i) const
{
    return mPoint[i];
}
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Explicit instantiation.
//----------------------------------------------------------------------------
template WM5_MATHEMATICS_ITEM
class IntrSegment2Segment2<float>;

template WM5_MATHEMATICS_ITEM
class IntrSegment2Segment2<double>;
//----------------------------------------------------------------------------
}
