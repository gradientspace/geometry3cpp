// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.10.0 (2013/03/17)

#ifndef WM5DISTPOINT2HYPERBOLA2_H
#define WM5DISTPOINT2HYPERBOLA2_H

#include "Wm5MathematicsLIB.h"
#include "Wm5Vector2.h"

namespace Wm5
{

// NOTE: This class does not derive from Distance<Real,Vector2<Real>> as do
// other distance calculators.  Wild Magic 6 will not use this paradigm.

// Compute the distance between the (x0,y0) and the standard hyperbola
// F(x,y) = (x/a)^2 - (y/b)^2 - 1 = 0, where a > 0 and b > 0.  Consider
// the case when (x0,y0) is in the first quadrant, so x0 >= 0 and
// y0 >= 0.  The closest point (x,y) is also in the first quadrant.  From
// the geometry, the vector D = (x0-x,y0-y) must be normal to the
// hyperbola at (x,y).  A normal is N = gradient(Q)/2 = (x/a^2,-y/b^2).
//
// Solution 1 (reduction to quartic polynomial):
// Because D and N are parallel, Dot(D,Perp(N)) = 0, which leads to the
// equation G(x,y) = (x0-x)*y/b^2 + (y0-y)*x/a^2 = 0.  The variable y
// may be eliminated from the two equations F(x,y) = 0 and G(x,y) = 0 to
// obtain
//   P(x) = c0 + c1*x + c2*x^2 + c3*x^3 + c4*x^4 = 0
// where
//   c0 = -a^6*x0^2
//   c1 = 2*a^4*(a^2+b^2)*x0,
//   c2 = -a^2*((a^2+b^2)^2 - a^2*x0^2 + b^2*y0^2)
//   c3 = -2*a*(a^2+b^2)*x0
//   c4 = (a^2+b^2)^2
// You may compute the roots of P(x) = 0.  For each root x = r, solve
// F(r,y) = 0 for y = s.  Compute the distance from (x0,y0) to (r,s).
// Of all distances from such pairs (r,s), choose the minimum distance
// as the distance from (x0,y0) to the hyperbola.  Numerically, root
// finding is not robust, so be wary of implementing the algorithm this
// way.
//
// Solution2 (robust):
// Because D and N are parallel, there must exist a scalar t for which
// D = t*N; that is, (x0-x,y0-y) = t*(x/a^2,-y/b^2).  Although not needed
// in the algorithm, by the geometry: sign(t) = sign(F(x0,y0)).  Some
// algebra leads to
//   (x,y) = x0/(1 + t/a^2), y0/(1 - t/b^2))
// Knowing that (x,y) is in the first quadrant, the t-value must satisfy
// -a^2 <= t <= b^2.  Substituting the (x,y) into F(x,y) = 0 leads to
//   H(t) = x0^2/(1 + t/a^2)^2 - y0^2/(1 - t/b^2)^2 - 1 = 0
// The graph of H(t) has vertical asymptotes at t = -a^2 and t = b^2.
// We care only about the graph for -a^2 <= t <= b^2.  Notice that
// as t approaches -a^2 from the right, H(t) goes to +infinity.  As t
// approaches b^2 from the left, H(t) goes to -infinity.  The derivative
// H'(t) is zero at most once; in fact, you can solve H'(t) = 0 in
// closed form.  This condition and the behavior of H(t) at the
// asymptotes guarantees that H(t) is a decreasing function for
// -a^2 <= t <= b^2, in which case H(t) = 0 has a unique root on this
// interval. Thus, you may use bisection to compute the root robustly.
// This is the implementation we use.

// 'point' is (x0,y0), 'extent' is (a,b), 'closest' is (x,y).  The return
// value is the distance.
template <typename Real> WM5_MATHEMATICS_ITEM
Real ComputeDistancePointToHyperbola (const Vector2<Real>& point,
    const Vector2<Real>& extent, Vector2<Real>& closest);

}

#endif
