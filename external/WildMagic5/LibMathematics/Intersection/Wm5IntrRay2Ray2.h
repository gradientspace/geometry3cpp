// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.2 (2013/11/12)

#ifndef WM5INTRRAY2RAY2_H
#define WM5INTRRAY2RAY2_H

#include "Wm5MathematicsLIB.h"
#include "Wm5Intersector.h"
#include "Wm5Ray2.h"

namespace Wm5
{

template <typename Real>
class WM5_MATHEMATICS_ITEM IntrRay2Ray2
    : public Intersector<Real,Vector2<Real> >
{
public:
    IntrRay2Ray2 (const Ray2<Real>& ray0, const Ray2<Real>& ray1);

    // Object access.  Ray0 is P0+t*D0 and Ray1 is P1+t*D1.
    const Ray2<Real>& GetRay0 () const;
    const Ray2<Real>& GetRay1 () const;

    // Static intersection query.
    virtual bool Test ();
    virtual bool Find ();

    // The computation for determining whether the linear components are
    // parallel might contain small floating-point round-off errors.  The
    // default threshold is Math<Real>::ZERO_TOLERANCE.  If you set the value,
    // pass in a nonnegative number.
    void SetDotThreshold (Real dotThreshold);
    Real GetDotThreshold () const;

    // The intersection set.  Let q = GetQuantity().  The cases are
    //
    //   q = 0: The rays do not intersect.  GetIntersection() returns IT_EMPTY.
    //
    //   q = 1: The rays intersect in a single point.  GetIntersection()
    //          returns IT_POINT.  For Find() queries, access the intersection
    //          point using GetPoint(0).
    //          
    //   q = 2: The rays are collinear and intersect in a segment.  This case
    //          happens only when D1 = -D0.  GetIntersection() returns
    //          IT_SEGMENT.  For Find() queries, access the endpoints of the
    //          intersection segment using GetPoint(0) and GetPoint(1).
    //
    //   q = INT_MAX:  The rays are collinear and intersect in a ray.  This
    //          case happens only when D1 = D0.  GetIntersection() returns
    //          IT_RAY.  For Find() queries, access the origin of the
    //          intersection ray using GetPoint(0).  The direction is
    //          D0 (= D1).
    int GetQuantity () const;
    const Vector2<Real>& GetPoint (int i) const;

private:
    using Intersector<Real,Vector2<Real> >::IT_EMPTY;
    using Intersector<Real,Vector2<Real> >::IT_POINT;
    using Intersector<Real,Vector2<Real> >::IT_SEGMENT;
    using Intersector<Real,Vector2<Real> >::IT_RAY;
    using Intersector<Real,Vector2<Real> >::IT_LINE;
    using Intersector<Real,Vector2<Real> >::mIntersectionType;

    // The objects to intersect.
    const Ray2<Real>* mRay0;
    const Ray2<Real>* mRay1;

    // See the comments before {Set,Get}DotThreshold.
    Real mDotThreshold;

    // Information about the intersection set.
    int mQuantity;
    Vector2<Real> mPoint[2];
};

typedef IntrRay2Ray2<float> IntrRay2Ray2f;
typedef IntrRay2Ray2<double> IntrRay2Ray2d;

}

#endif
