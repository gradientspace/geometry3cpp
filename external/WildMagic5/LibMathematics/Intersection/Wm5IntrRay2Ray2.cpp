// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.3 (2013/11/19)

#include "Wm5MathematicsPCH.h"
#include "Wm5IntrRay2Ray2.h"
#include "Wm5IntrLine2Line2.h"

namespace Wm5
{
//----------------------------------------------------------------------------
template <typename Real>
IntrRay2Ray2<Real>::IntrRay2Ray2 (const Ray2<Real>& ray0,
    const Ray2<Real>& ray1)
    :
    mRay0(&ray0),
    mRay1(&ray1),
    mDotThreshold(Math<Real>::ZERO_TOLERANCE)
{
}
//----------------------------------------------------------------------------
template <typename Real>
const Ray2<Real>& IntrRay2Ray2<Real>::GetRay0 () const
{
    return *mRay0;
}
//----------------------------------------------------------------------------
template <typename Real>
const Ray2<Real>& IntrRay2Ray2<Real>::GetRay1 () const
{
    return *mRay1;
}
//----------------------------------------------------------------------------
template <typename Real>
bool IntrRay2Ray2<Real>::Test ()
{
    Real s[2];
    mIntersectionType = IntrLine2Line2<Real>::Classify(mRay0->Origin,
        mRay0->Direction, mRay1->Origin, mRay1->Direction, mDotThreshold, s);

    if (mIntersectionType == IT_POINT)
    {
        // Test whether the line-line intersection is on the rays.
        if (s[0] >= (Real)0 && s[1] >= (Real)0)
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
        if (mRay0->Direction.Dot(mRay1->Direction) > (Real)0)
        {
            // The rays are collinear and in the same direction, so they must
            // overlap.
            mQuantity = INT_MAX;
            mIntersectionType = IT_RAY;
        }
        else
        {
            // The rays are collinear but in opposite directions.  Test
            // whether they overlap.  Ray0 has interval [0,+infinity) and
            // ray1 has interval (-infinity,t1] relative to ray0.direction.
            Real t1 = mRay0->Direction.Dot(mRay1->Origin - mRay0->Origin);
            if (t1 > (Real)0)
            {
                mQuantity = 2;
                mIntersectionType = IT_SEGMENT;
            }
            else if (t1 < (Real)0)
            {
                mQuantity = 0;
                mIntersectionType = IT_EMPTY;
            }
            else  // t1 == 0
            {
                mQuantity = 1;
                mIntersectionType = IT_POINT;
            }
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
bool IntrRay2Ray2<Real>::Find ()
{
    Real s[2];
    mIntersectionType = IntrLine2Line2<Real>::Classify(mRay0->Origin,
        mRay0->Direction, mRay1->Origin, mRay1->Direction, mDotThreshold, s);

    if (mIntersectionType == IT_POINT)
    {
        // Test whether the line-line intersection is on the rays.
        if (s[0] >= (Real)0 && s[1] >= (Real)0)
        {
            mQuantity = 1;
            mPoint[0] = mRay0->Origin + s[0]*mRay0->Direction;
        }
        else
        {
            mQuantity = 0;
            mIntersectionType = IT_EMPTY;
        }
    }
    else if (mIntersectionType == IT_LINE)
    {
        // Compute t1 for which ray1.origin = ray0.origin + t1*ray0.direction.
        Real t1 = mRay0->Direction.Dot(mRay1->Origin - mRay0->Origin);
        if (mRay0->Direction.Dot(mRay1->Direction) > (Real)0)
        {
            // The rays are collinear and in the same direction, so they must
            // overlap.  Effectively, we compute the intersection of intervals
            // [0,+infinity) and [t1,+infinity).
            mQuantity = INT_MAX;
            mIntersectionType = IT_RAY;
            mPoint[0] = (t1 > (Real)0 ? mRay1->Origin : mRay0->Origin);
        }
        else
        {
            // The rays are collinear but in opposite directions.  Test
            // whether they overlap.  Ray0 has interval [0,+infinity) and
            // ray1 has interval (-infinity,t1].
            if (t1 > (Real)0)
            {
                mQuantity = 2;
                mIntersectionType = IT_SEGMENT;
                mPoint[0] = mRay0->Origin;
                mPoint[1] = mRay1->Origin;
            }
            else if (t1 < (Real)0)
            {
                mQuantity = 0;
                mIntersectionType = IT_EMPTY;
            }
            else  // t1 == 0
            {
                mQuantity = 1;
                mIntersectionType = IT_POINT;
                mPoint[0] = mRay0->Origin;
            }
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
void IntrRay2Ray2<Real>::SetDotThreshold (Real dotThreshold)
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
Real IntrRay2Ray2<Real>::GetDotThreshold () const
{
    return mDotThreshold;
}
//----------------------------------------------------------------------------
template <typename Real>
int IntrRay2Ray2<Real>::GetQuantity () const
{
    return mQuantity;
}
//----------------------------------------------------------------------------
template <typename Real>
const Vector2<Real>& IntrRay2Ray2<Real>::GetPoint (int i) const
{
    return mPoint[i];
}
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Explicit instantiation.
//----------------------------------------------------------------------------
template WM5_MATHEMATICS_ITEM
class IntrRay2Ray2<float>;

template WM5_MATHEMATICS_ITEM
class IntrRay2Ray2<double>;
//----------------------------------------------------------------------------
}
