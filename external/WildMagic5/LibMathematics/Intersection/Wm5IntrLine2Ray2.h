// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.2 (2013/11/12)

#ifndef WM5INTRLINE2RAY2_H
#define WM5INTRLINE2RAY2_H

#include "Wm5MathematicsLIB.h"
#include "Wm5Intersector.h"
#include "Wm5Line2.h"
#include "Wm5Ray2.h"

namespace Wm5
{

template <typename Real>
class WM5_MATHEMATICS_ITEM IntrLine2Ray2
    : public Intersector<Real,Vector2<Real> >
{
public:
    IntrLine2Ray2 (const Line2<Real>& line, const Ray2<Real>& ray);

    // Object access.
    const Line2<Real>& GetLine () const;
    const Ray2<Real>& GetRay () const;

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
    //   q = 0: The line/ray do not intersect.  GetIntersection() returns
    //          IT_EMPTY.
    //
    //   q = 1: The line/ray intersect in a single point.  GetIntersection()
    //          returns IT_POINT.  For Find() queries, access the intersection
    //          point using GetPoint().
    //          
    //   q = INT_MAX:  The line/ray are collinear.  GetIntersection() returns
    //          IT_RAY.  GetPoint() should not be called for Find() queries.
    int GetQuantity () const;
    const Vector2<Real>& GetPoint () const;

private:
    using Intersector<Real,Vector2<Real> >::IT_EMPTY;
    using Intersector<Real,Vector2<Real> >::IT_POINT;
    using Intersector<Real,Vector2<Real> >::IT_RAY;
    using Intersector<Real,Vector2<Real> >::IT_LINE;
    using Intersector<Real,Vector2<Real> >::mIntersectionType;

    // The objects to intersect.
    const Line2<Real>* mLine;
    const Ray2<Real>* mRay;

    // See the comments before {Set,Get}DotThreshold.
    Real mDotThreshold;

    // Information about the intersection set.
    int mQuantity;
    Vector2<Real> mPoint;
};

typedef IntrLine2Ray2<float> IntrLine2Ray2f;
typedef IntrLine2Ray2<double> IntrLine2Ray2d;

}

#endif
