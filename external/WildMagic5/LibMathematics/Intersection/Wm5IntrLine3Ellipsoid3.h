// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.2 (2012/11/03)

#ifndef WM5INTRLINE3ELLIPSOID3_H
#define WM5INTRLINE3ELLIPSOID3_H

#include "Wm5MathematicsLIB.h"
#include "Wm5Intersector.h"
#include "Wm5Line3.h"
#include "Wm5Ellipsoid3.h"

namespace Wm5
{

template <typename Real>
class WM5_MATHEMATICS_ITEM IntrLine3Ellipsoid3
    : public Intersector<Real,Vector3<Real> >
{
public:
    IntrLine3Ellipsoid3 (const Line3<Real>& line,
        const Ellipsoid3<Real>& ellipsoid);

    // Object access.
    const Line3<Real>& GetLine () const;
    const Ellipsoid3<Real>& GetEllipsoid () const;

    // Static intersection queries.
    virtual bool Test ();
    virtual bool Find ();

    // The intersection set.
    int GetQuantity () const;
    const Vector3<Real>& GetPoint (int i) const;

    // Small thresholds are used for testing the discriminant of the quadratic
    // equation related to the computations: Q(t) = a2*t^2 + 2*a1*t + a0.  The
    // discriminant is D = a1*a1 - a0*a2.  Q(t) has no real-valued roots when
    // D < 0, one real-valued root when D = 0, or two real-valued roots when
    // D > 0.  The code logic involves user-defined thresholds:
    //   if (D < negThreshold) { no roots (no intersections) }
    //   else if (D > posThreshold) { two roots (two intersections) }
    //   else { one root (one intersection) }
    // The default values for the thresholds are zero, but you may set them
    // to be nonzero (negThreshold <= 0 and posThreshold >= 0).  Previously,
    // the negative threshold was hard-coded as zero.  The positive threshold
    // was hard-coded to Math<Real>::ZERO_TOLERANCE, which is not suitable for
    // some data sets (i.e. when ellipsoid extents are quite large).  The
    // default is now zero, so if your application relied on the old behavior,
    // you must modify this value.
    void SetNegativeThreshold (Real negThreshold);
    Real GetNegativeThreshold () const;
    void SetPositiveThreshold (Real posThreshold);
    Real GetPositiveThreshold () const;

private:
    using Intersector<Real,Vector3<Real> >::IT_EMPTY;
    using Intersector<Real,Vector3<Real> >::IT_POINT;
    using Intersector<Real,Vector3<Real> >::IT_SEGMENT;
    using Intersector<Real,Vector3<Real> >::mIntersectionType;

    // The objects to intersect.
    const Line3<Real>* mLine;
    const Ellipsoid3<Real>* mEllipsoid;

    // Information about the intersection set.
    int mQuantity;
    Vector3<Real> mPoint[2];

    // For testing the discriminant.  The default values are zero.  You may
    // set the negative threshold to a (small) negative number and the
    // positive threshold to a  (small) positive number.
    Real mNegativeThreshold;
    Real mPositiveThreshold;
};

typedef IntrLine3Ellipsoid3<float> IntrLine3Ellipsoid3f;
typedef IntrLine3Ellipsoid3<double> IntrLine3Ellipsoid3d;

}

#endif
