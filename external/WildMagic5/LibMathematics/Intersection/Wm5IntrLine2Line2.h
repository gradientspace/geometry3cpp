// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.2 (2013/11/12)

#ifndef WM5INTRLINE2LINE2_H
#define WM5INTRLINE2LINE2_H

#include "Wm5MathematicsLIB.h"
#include "Wm5Intersector.h"
#include "Wm5Line2.h"

namespace Wm5
{

template <typename Real>
class WM5_MATHEMATICS_ITEM IntrLine2Line2
    : public Intersector<Real,Vector2<Real> >
{
public:
    IntrLine2Line2 (const Line2<Real>& line0, const Line2<Real>& line1);

    // Object access.
    const Line2<Real>& GetLine0 () const;
    const Line2<Real>& GetLine1 () const;

    // Static intersection query.
    virtual bool Test ();
    virtual bool Find ();

    // Classification of linear components.  The lines are P0+s*D0 and P1+s*D1.
    // The return value is IT_EMPTY when the lines do not intersect, IT_POINT
    // when the lines intersect in a single point, or IT_LINE when the lines
    // are collinear.  When the intersection is a single point, the point is
    // P0+s[0]*D0 = P1+s[1]*D1.  The 'dotThreshold' parameter is a nonnegative
    // number used for testing for parallel or perpendicular vectors.  The
    // value used in the Test() and Find() functions is mDotThreshold.  If you
    // want to know the s[] values, pass a nonnull array of two elements for
    // 's'.  The code is shared by the intersectors for ray-ray, ray-segment,
    // and segment-segment.
    static int Classify (const Vector2<Real>& P0, const Vector2<Real>& D0,
        const Vector2<Real>& P1, const Vector2<Real>& D1, Real dotThreshold,
        Real* s = 0);

    // The computation for determining whether the linear components are
    // parallel might contain small floating-point round-off errors.  The
    // default threshold is Math<Real>::ZERO_TOLERANCE.  If you set the value,
    // pass in a nonnegative number.
    void SetDotThreshold (Real dotThreshold);
    Real GetDotThreshold () const;

    // The intersection set.  Let q = GetQuantity().  The cases are
    //
    //   q = 0: The lines do not intersect.  GetIntersection() returns
    //          IT_EMPTY.
    //
    //   q = 1: The lines intersect in a single point.  GetIntersection()
    //          returns IT_POINT.  For Find() queries, access the intersection
    //          point using GetPoint().
    //          
    //   q = INT_MAX:  The lines are collinear.  GetIntersection() returns
    //          IT_LINE.  GetPoint() should not be called for Find() queries.
    int GetQuantity () const;
    const Vector2<Real>& GetPoint () const;

private:
    using Intersector<Real,Vector2<Real> >::IT_EMPTY;
    using Intersector<Real,Vector2<Real> >::IT_POINT;
    using Intersector<Real,Vector2<Real> >::IT_LINE;
    using Intersector<Real,Vector2<Real> >::mIntersectionType;

    // The objects to intersect.
    const Line2<Real>* mLine0;
    const Line2<Real>* mLine1;

    // See the comments before {Set,Get}DotThreshold.
    Real mDotThreshold;

    // Information about the intersection set.
    int mQuantity;
    Vector2<Real> mPoint;
};

typedef IntrLine2Line2<float> IntrLine2Line2f;
typedef IntrLine2Line2<double> IntrLine2Line2d;

}

#endif
