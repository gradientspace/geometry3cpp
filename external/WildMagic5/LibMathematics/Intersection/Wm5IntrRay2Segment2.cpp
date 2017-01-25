// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.2 (2013/11/12)

#include "Wm5MathematicsPCH.h"
#include "Wm5IntrRay2Segment2.h"
#include "Wm5IntrLine2Line2.h"
#include "Wm5Intersector1.h"

namespace Wm5
{
//----------------------------------------------------------------------------
template <typename Real>
IntrRay2Segment2<Real>::IntrRay2Segment2 (const Ray2<Real>& ray,
    const Segment2<Real>& segment)
    :
    mRay(&ray),
    mSegment(&segment),
    mIntervalThreshold((Real)0),
    mDotThreshold(Math<Real>::ZERO_TOLERANCE)
{
}
//----------------------------------------------------------------------------
template <typename Real>
const Ray2<Real>& IntrRay2Segment2<Real>::GetRay () const
{
    return *mRay;
}
//----------------------------------------------------------------------------
template <typename Real>
const Segment2<Real>& IntrRay2Segment2<Real>::GetSegment () const
{
    return *mSegment;
}
//----------------------------------------------------------------------------
template <typename Real>
bool IntrRay2Segment2<Real>::Test ()
{
    Real s[2];
    mIntersectionType = IntrLine2Line2<Real>::Classify(mRay->Origin,
        mRay->Direction, mSegment->Center, mSegment->Direction,
        mDotThreshold, s);

    if (mIntersectionType == IT_POINT)
    {
        // Test whether the line-line intersection is on the ray and on the
        // segment.
        if (s[0] >= (Real)0
        &&  Math<Real>::FAbs(s[1]) <= mSegment->Extent + mIntervalThreshold)
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
        // Compute the location of the segment center relative to the ray.
        Real t1 = mRay->Direction.Dot(mSegment->Center - mRay->Origin);

        // Compute the location of the right-most point of the segment
        // relative to the ray direction.
        Real tmax = t1 + mSegment->Extent;
        if (tmax > (Real)0)
        {
            mQuantity = 2;
            mIntersectionType = IT_SEGMENT;
        }
        else if (tmax < (Real)0)
        {
            mQuantity = 0;
            mIntersectionType = IT_EMPTY;
        }
        else  // tmax == 0
        {
            mQuantity = 1;
            mIntersectionType = IT_POINT;
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
bool IntrRay2Segment2<Real>::Find ()
{
    Real s[2];
    mIntersectionType = IntrLine2Line2<Real>::Classify(mRay->Origin,
        mRay->Direction, mSegment->Center, mSegment->Direction,
        mDotThreshold, s);

    if (mIntersectionType == IT_POINT)
    {
        // Test whether the line-line intersection is on the ray and on the
        // segment.
        if (s[0] >= (Real)0
        &&  Math<Real>::FAbs(s[1]) <= mSegment->Extent + mIntervalThreshold)
        {
            mQuantity = 1;
            mPoint[0] = mRay->Origin + s[0]*mRay->Direction;
        }
        else
        {
            mQuantity = 0;
            mIntersectionType = IT_EMPTY;
        }
    }
    else if (mIntersectionType == IT_LINE)
    {
        // Compute the location of the segment center relative to the ray.
        Real t1 = mRay->Direction.Dot(mSegment->Center - mRay->Origin);

        // Compute the location of the segment endpoints relative to the
        // ray direction.
        Real tmin = t1 - mSegment->Extent;
        Real tmax = t1 + mSegment->Extent;

        // Compute the intersection of [0,+infinity) and [tmin,tmax].
        Intersector1<Real> calc((Real)0, Math<Real>::MAX_REAL, tmin, tmax);
        calc.Find();
        mQuantity = calc.GetNumIntersections();
        if (mQuantity == 2)
        {
            mIntersectionType = IT_SEGMENT;
            mPoint[0] = mRay->Origin + calc.GetIntersection(0)*mRay->Direction;
            mPoint[1] = mRay->Origin + calc.GetIntersection(1)*mRay->Direction;
        }
        else if (mQuantity == 1)
        {
            mIntersectionType = IT_POINT;
            mPoint[0] = mRay->Origin;
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
void IntrRay2Segment2<Real>::SetIntervalThreshold (Real intervalThreshold)
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
Real IntrRay2Segment2<Real>::GetIntervalThreshold () const
{
    return mIntervalThreshold;
}
//----------------------------------------------------------------------------
template <typename Real>
void IntrRay2Segment2<Real>::SetDotThreshold (Real dotThreshold)
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
Real IntrRay2Segment2<Real>::GetDotThreshold () const
{
    return mDotThreshold;
}
//----------------------------------------------------------------------------
template <typename Real>
int IntrRay2Segment2<Real>::GetQuantity () const
{
    return mQuantity;
}
//----------------------------------------------------------------------------
template <typename Real>
const Vector2<Real>& IntrRay2Segment2<Real>::GetPoint (int i) const
{
    return mPoint[i];
}
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Explicit instantiation.
//----------------------------------------------------------------------------
template WM5_MATHEMATICS_ITEM
class IntrRay2Segment2<float>;

template WM5_MATHEMATICS_ITEM
class IntrRay2Segment2<double>;
//----------------------------------------------------------------------------
}
