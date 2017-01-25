// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.10.0 (2013/03/17)

#include "Wm5MathematicsPCH.h"
#include "Wm5DistPoint2Hyperbola2.h"

namespace Wm5
{
//----------------------------------------------------------------------------
template <typename Real>
static Real H (Real t, const Vector2<Real>& point, const Vector2<Real>& resqr)
{
    Real ratio0 = point[0]/((Real)1 + t*resqr[0]);
    Real ratio1 = point[1]/((Real)1 - t*resqr[1]);
    return ratio0*ratio0 - ratio1*ratio1 - (Real)1;
}
//----------------------------------------------------------------------------
template <typename Real>
Real ComputeDistancePointToHyperbola (const Vector2<Real>& point,
    const Vector2<Real>& extent, Vector2<Real>& closest)
{
    assertion(extent[0] > (Real)0 && extent[1] > (Real)0, "Invalid inputs");

    Vector2<Real> esqr(extent[0]*extent[0], extent[1]*extent[1]);
    Vector2<Real> resqr(((Real)1)/esqr[0], ((Real)1)/esqr[1]);

    // Initialize for bisection.  It is not relevant that H(-a^2) = +infinity
    // and H(b^2) = -infinity, so we need only initialize function values with
    // accordingly signed numbers.
    Real t0 = -esqr[0], t1 = esqr[1];
    Real troot = ((Real)0.5)*(t0 + t1);
    Real hroot = H(troot, point, resqr);

    // Iterate until H(troot) is exactly zero or until one of the
    // floating-point endpoints does not change anymore.  The latter condition
    // takes advantage of the nature of IEEE floating-point numbers, so the
    // loop must terminate in a finite number of steps.
    while (hroot != (Real)0 && troot != t0 && troot != t1)
    {
        if (hroot > (Real)0)
        {
            t0 = troot;
            troot = ((Real)0.5)*(t0 + t1);
        }
        else // hroot < (Real)0
        {
            t1 = troot;
        }

        troot = ((Real)0.5)*(t0 + t1);
        hroot = H(troot, point, resqr);
    }

    closest[0] = point[0]/((Real)1 + troot*resqr[0]);
    closest[1] = point[1]/((Real)1 - troot*resqr[1]);

    Vector2<Real> diff = point - closest;
    Real distance = diff.Length();
    return distance;
}
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Explicit instantiation.
//----------------------------------------------------------------------------
template WM5_MATHEMATICS_ITEM
float ComputeDistancePointToHyperbola<float> (const Vector2<float>&,
    const Vector2<float>&, Vector2<float>&);

template WM5_MATHEMATICS_ITEM
double ComputeDistancePointToHyperbola<double> (const Vector2<double>&,
    const Vector2<double>&, Vector2<double>&);
//----------------------------------------------------------------------------
}
