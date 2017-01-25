// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.2 (2013/07/28)

#ifndef WM5DISTLINE3CIRCLE3_H
#define WM5DISTLINE3CIRCLE3_H

#include "Wm5MathematicsLIB.h"
#include "Wm5Distance.h"
#include "Wm5Line3.h"
#include "Wm5Circle3.h"

namespace Wm5
{

template <typename Real>
class WM5_MATHEMATICS_ITEM DistLine3Circle3
    : public Distance<Real,Vector3<Real> >
{
public:
    DistLine3Circle3 (const Line3<Real>& line,
        const Circle3<Real>& circle);

    // Object access.
    const Line3<Real>& GetLine () const;
    const Circle3<Real>& GetCircle () const;

    // Static distance queries.  Compute the distance from the point P to the
    // circle.  When P is on the normal line C+t*N where C is the circle
    // center and N is the normal to the plane containing the circle, then
    // all circle points are equidistant from P.  In this case the returned
    // point is C+r*U, where U is a vector perpendicular to N.
    virtual Real Get ();
    virtual Real GetSquared ();

    // Member access.  The possible combinations for number of line-circle
    // closest points is (1,1), (2,2), or (1,INT_MAX).
    int GetNumClosestLine () const;
    const Vector3<Real>& GetClosestLine (int i) const;
    int GetNumClosestCircle () const;
    const Vector3<Real>& GetClosestCircle (int i) const;

    // Function calculations for dynamic distance queries.
    virtual Real Get (Real t, const Vector3<Real>& velocity0,
        const Vector3<Real>& velocity1);
    virtual Real GetSquared (Real t, const Vector3<Real>& velocity0,
        const Vector3<Real>& velocity1);

private:
    using Distance<Real,Vector3<Real> >::mClosestPoint0;
    using Distance<Real,Vector3<Real> >::mClosestPoint1;
    using Distance<Real,Vector3<Real> >::mHasMultipleClosestPoints0;
    using Distance<Real,Vector3<Real> >::mHasMultipleClosestPoints1;

    // The mClosestLine[i] is an input.  The mClosestCircle[i] is an output.
    // The squared distance between these is returned.  The value of
    // mNumClosest is set to 1 or INT_MAX.  If the latter, all circle points
    // are equidistance from mClosestLine[i].
    Real SqrDistancePointCircle (int i);

    // Bisect the function F(s) = s + m2b2 - r*m0sqr*s/sqrt(m0sqr*s*s + b1sqr)
    // on the specified interval [smin,smax].
    static Real BisectF (Real m2b2, Real rm0sqr, Real m0sqr, Real b1sqr,
        Real smin, Real smax);

    const Line3<Real>* mLine;
    const Circle3<Real>* mCircle;
    int mNumClosestLine;
    Vector3<Real> mClosestLine[2];
    int mNumClosestCircle;
    Vector3<Real> mClosestCircle[2];
};

typedef DistLine3Circle3<float> DistLine3Circle3f;
typedef DistLine3Circle3<double> DistLine3Circle3d;

}

#endif
