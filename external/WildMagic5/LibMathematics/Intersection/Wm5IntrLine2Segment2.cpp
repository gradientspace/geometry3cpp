// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.2 (2013/11/12)

#include "Wm5MathematicsPCH.h"
#include "Wm5IntrLine2Segment2.h"
#include "Wm5IntrLine2Line2.h"

namespace Wm5
{
//----------------------------------------------------------------------------
template <typename Real>
IntrLine2Segment2<Real>::IntrLine2Segment2 (const Line2<Real>& line,
    const Segment2<Real>& segment)
    :
    mLine(&line),
    mSegment(&segment),
    mIntervalThreshold((Real)0),
    mDotThreshold(Math<Real>::ZERO_TOLERANCE)
{
}
//----------------------------------------------------------------------------
template <typename Real>
const Line2<Real>& IntrLine2Segment2<Real>::GetLine () const
{
    return *mLine;
}
//----------------------------------------------------------------------------
template <typename Real>
const Segment2<Real>& IntrLine2Segment2<Real>::GetSegment () const
{
    return *mSegment;
}
//----------------------------------------------------------------------------
template <typename Real>
bool IntrLine2Segment2<Real>::Test ()
{
    Real s[2];
    mIntersectionType = IntrLine2Line2<Real>::Classify(mLine->Origin,
        mLine->Direction, mSegment->Center, mSegment->Direction,
        mDotThreshold, s);

    if (mIntersectionType == IT_POINT)
    {
        // Test whether the line-line intersection is on the segment.
        if (Math<Real>::FAbs(s[1]) <= mSegment->Extent + mIntervalThreshold)
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
        mIntersectionType = IT_SEGMENT;
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
bool IntrLine2Segment2<Real>::Find ()
{
    Real s[2];
    mIntersectionType = IntrLine2Line2<Real>::Classify(mLine->Origin,
        mLine->Direction, mSegment->Center, mSegment->Direction,
        mDotThreshold, s);

    if (mIntersectionType == IT_POINT)
    {
        // Test whether the line-line intersection is on the segment.
        if (Math<Real>::FAbs(s[1]) <= mSegment->Extent + mIntervalThreshold)
        {
            mQuantity = 1;
            mPoint = mLine->Origin + s[0]*mLine->Direction;
        }
        else
        {
            mQuantity = 0;
            mIntersectionType = IT_EMPTY;
        }
    }
    else if (mIntersectionType == IT_LINE)
    {
        mIntersectionType = IT_SEGMENT;
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
void IntrLine2Segment2<Real>::SetIntervalThreshold (Real intervalThreshold)
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
Real IntrLine2Segment2<Real>::GetIntervalThreshold () const
{
    return mIntervalThreshold;
}
//----------------------------------------------------------------------------
template <typename Real>
void IntrLine2Segment2<Real>::SetDotThreshold (Real dotThreshold)
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
Real IntrLine2Segment2<Real>::GetDotThreshold () const
{
    return mDotThreshold;
}
//----------------------------------------------------------------------------
template <typename Real>
int IntrLine2Segment2<Real>::GetQuantity () const
{
    return mQuantity;
}
//----------------------------------------------------------------------------
template <typename Real>
const Vector2<Real>& IntrLine2Segment2<Real>::GetPoint () const
{
    return mPoint;
}
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Explicit instantiation.
//----------------------------------------------------------------------------
template WM5_MATHEMATICS_ITEM
class IntrLine2Segment2<float>;

template WM5_MATHEMATICS_ITEM
class IntrLine2Segment2<double>;
//----------------------------------------------------------------------------
}
