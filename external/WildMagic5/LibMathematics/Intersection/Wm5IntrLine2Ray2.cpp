// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.2 (2013/11/12)

#include "Wm5MathematicsPCH.h"
#include "Wm5IntrLine2Ray2.h"
#include "Wm5IntrLine2Line2.h"

namespace Wm5
{
//----------------------------------------------------------------------------
template <typename Real>
IntrLine2Ray2<Real>::IntrLine2Ray2 (const Line2<Real>& line,
    const Ray2<Real>& ray)
    :
    mLine(&line),
    mRay(&ray),
    mDotThreshold(Math<Real>::ZERO_TOLERANCE)
{
}
//----------------------------------------------------------------------------
template <typename Real>
const Line2<Real>& IntrLine2Ray2<Real>::GetLine () const
{
    return *mLine;
}
//----------------------------------------------------------------------------
template <typename Real>
const Ray2<Real>& IntrLine2Ray2<Real>::GetRay () const
{
    return *mRay;
}
//----------------------------------------------------------------------------
template <typename Real>
bool IntrLine2Ray2<Real>::Test ()
{
    Real s[2];
    mIntersectionType = IntrLine2Line2<Real>::Classify(mLine->Origin,
        mLine->Direction, mRay->Origin, mRay->Direction, mDotThreshold, s);

    if (mIntersectionType == IT_POINT)
    {
        // Test whether the line-line intersection is on the ray.
        if (s[1] >= (Real)0)
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
        mIntersectionType = IT_RAY;
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
bool IntrLine2Ray2<Real>::Find ()
{
    Real s[2];
    mIntersectionType = IntrLine2Line2<Real>::Classify(mLine->Origin,
        mLine->Direction, mRay->Origin, mRay->Direction, mDotThreshold, s);

    if (mIntersectionType == IT_POINT)
    {
        // Test whether the line-line intersection is on the ray.
        if (s[1] >= (Real)0)
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
        mIntersectionType = IT_RAY;
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
void IntrLine2Ray2<Real>::SetDotThreshold (Real dotThreshold)
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
Real IntrLine2Ray2<Real>::GetDotThreshold () const
{
    return mDotThreshold;
}
//----------------------------------------------------------------------------
template <typename Real>
int IntrLine2Ray2<Real>::GetQuantity () const
{
    return mQuantity;
}
//----------------------------------------------------------------------------
template <typename Real>
const Vector2<Real>& IntrLine2Ray2<Real>::GetPoint () const
{
    return mPoint;
}
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Explicit instantiation.
//----------------------------------------------------------------------------
template WM5_MATHEMATICS_ITEM
class IntrLine2Ray2<float>;

template WM5_MATHEMATICS_ITEM
class IntrLine2Ray2<double>;
//----------------------------------------------------------------------------
}
