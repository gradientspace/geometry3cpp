// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.3 (2013/11/12)

#ifndef WM5INTRSEGMENT2SEGMENT2_H
#define WM5INTRSEGMENT2SEGMENT2_H

#include "Wm5MathematicsLIB.h"
#include "Wm5Intersector.h"
#include "Wm5Segment2.h"

namespace Wm5
{

template <typename Real>
class WM5_MATHEMATICS_ITEM IntrSegment2Segment2
    : public Intersector<Real,Vector2<Real> >
{
public:
    IntrSegment2Segment2 (const Segment2<Real>& segment0,
        const Segment2<Real>& segment1);

    // Object access.
    const Segment2<Real>& GetSegment0 () const;
    const Segment2<Real>& GetSegment1 () const;

    // Static intersection query.
    virtual bool Test ();
    virtual bool Find ();

    // The intersection testing uses the center-extent form for line segments.
    // If you start with endpoints (Vector2<Real>) and create Segment2<Real>
    // objects, the conversion to center-extent form can contain small
    // numerical round-off errors.  Testing for the intersection of two
    // segments that share an endpoint might lead to a failure due to the
    // round-off errors.  To allow for this, you may specify a small positive
    // threshold that slightly enlarges the intervals for the segments.  The
    // default value is zero.
    void SetIntervalThreshold (Real intervalThreshold);
    Real GetIntervalThreshold () const;

    // The computation for determining whether the linear components are
    // parallel might contain small floating-point round-off errors.  The
    // default threshold is Math<Real>::ZERO_TOLERANCE.  If you set the value,
    // pass in a nonnegative number.
    void SetDotThreshold (Real dotThreshold);
    Real GetDotThreshold () const;

    // The intersection set.  Let q = GetQuantity().  The cases are
    //
    //   q = 0: The segments do not intersect.  GetIntersection() returns
    //          IT_EMPTY.
    //
    //   q = 1: The segments intersect in a single point.  GetIntersection()
    //          returns IT_POINT.  For Find() queries, access the intersection
    //          point using GetPoint(0).
    //          
    //   q = 2: The segments are collinear and intersect in a segment.
    //          GetIntersection() returns IT_SEGMENT.  For Find() queries,
    //          access the endpoints of the intersection segment using
    //          GetPoint(0) and GetPoint(1).
    int GetQuantity () const;
    const Vector2<Real>& GetPoint (int i) const;

private:
    using Intersector<Real,Vector2<Real> >::IT_EMPTY;
    using Intersector<Real,Vector2<Real> >::IT_POINT;
    using Intersector<Real,Vector2<Real> >::IT_SEGMENT;
    using Intersector<Real,Vector2<Real> >::IT_LINE;
    using Intersector<Real,Vector2<Real> >::mIntersectionType;

    // The objects to intersect.
    const Segment2<Real>* mSegment0;
    const Segment2<Real>* mSegment1;

    // See the comments before {Set,Get}IntervalThreshold.
    Real mIntervalThreshold;

    // See the comments before {Set,Get}DotThreshold.
    Real mDotThreshold;

    // Information about the intersection set.
    int mQuantity;
    Vector2<Real> mPoint[2];
};

typedef IntrSegment2Segment2<float> IntrSegment2Segment2f;
typedef IntrSegment2Segment2<double> IntrSegment2Segment2d;

}

#endif
