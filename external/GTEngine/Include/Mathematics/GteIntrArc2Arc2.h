// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <Mathematics/GteArc2.h>
#include <Mathematics/GteIntrCircle2Circle2.h>

namespace gte
{

template <typename Real>
class FIQuery<Real, Arc2<Real>, Arc2<Real>>
{
public:
    struct Result
    {
        bool intersect;

        // The number of intersections is 0, 1, 2, or maxInt =
        // std::numeric_limits<int>::max().  When 1, the arcs intersect in a
        // single point.  When 2, the arcs are not co-circular and intersect
        // in two points.  When maxInt, the arcs are co-circular and
        // intersect.
        int numIntersections;

        // Valid only when numIntersections = 1 or 2.
        Vector2<Real> point[2];

        // Valid only when numIntersections = maxInt.
        Arc2<Real> arc;
    };

    Result operator()(Arc2<Real> const& arc0, Arc2<Real> const& arc1);
};


template <typename Real>
typename FIQuery<Real, Arc2<Real>, Arc2<Real>>::Result
FIQuery<Real, Arc2<Real>, Arc2<Real>>::operator()(
    Arc2<Real> const& arc0, Arc2<Real> const& arc1)
{
    Result result;

    Circle2<Real> circle0(arc0.center, arc0.radius);
    Circle2<Real> circle1(arc1.center, arc1.radius);
    FIQuery<Real, Circle2<Real>, Circle2<Real>> ccQuery;
    auto ccResult = ccQuery(circle0, circle1);
    if (!ccResult.intersect)
    {
        // The arcs do not intersect.
        result.intersect = false;
        result.numIntersections = 0;
        return result;
    }

    if (ccResult.numIntersections == std::numeric_limits<int>::max())
    {
        // Arcs are cocircular.  Determine if they overlap.  Let arc0 be
        // <A0,A1> and arc1 be <B0,B1>, the points ordered counterclockwise
        // around the circle of the arc.
        if (arc1.Contains(arc0.end[0]))
        {
            if (arc1.Contains(arc0.end[1]))
            {
                // Arc0 inside arc1, <B0,A0,A1,B1>.
                result.intersect = true;
                result.numIntersections = std::numeric_limits<int>::max();
                result.arc = arc0;
            }
            else
            {
                if (arc0.end[0] != arc1.end[1])
                {
                    // Arc0 and arc1 overlap, <B0,A0,B1,A1>.
                    result.intersect = true;
                    result.numIntersections = std::numeric_limits<int>::max();
                    result.arc = Arc2<Real>(arc0.center, arc0.radius,
                        arc0.end[0], arc1.end[1]);
                }
                else
                {
                    // Arc0 and arc1 share endpoint, <B0,A0,B1,A1>, A0 = B1.
                    result.intersect = true;
                    result.numIntersections = 1;
                    result.point[0] = arc0.end[0];
                }
            }
            return result;
        }

        if (arc1.Contains(arc0.end[1]))
        {
            if (arc0.end[1] != arc1.end[0])
            {
                // Arc0 and arc1 overlap, <A0,B0,A1,B1>.
                result.intersect = true;
                result.numIntersections = std::numeric_limits<int>::max();
                result.arc = Arc2<Real>(arc0.center, arc0.radius,
                    arc1.end[0], arc0.end[1]);
            }
            else
            {
                // Arc0 and arc1 share endpoint, <A0,B0,A1,B1>, B0 = A1.
                result.intersect = true;
                result.numIntersections = 1;
                result.point[0] = arc1.end[0];
            }
            return result;
        }

        if (arc0.Contains(arc1.end[0]))
        {
            // Arc1 inside arc0, <A0,B0,B1,A1>.
            result.intersect = true;
            result.numIntersections = std::numeric_limits<int>::max();
            result.arc = arc1;
        }
        else
        {
            // Arcs do not overlap, <A0,A1,B0,B1>.
            result.intersect = false;
            result.numIntersections = 0;
        }
        return result;
    }

    // Test whether circle-circle intersection points are on the arcs.
    for (int i = 0; i < ccResult.numIntersections; ++i)
    {
        result.numIntersections = 0;
        if (arc0.Contains(ccResult.point[i])
            && arc1.Contains(ccResult.point[i]))
        {
            result.point[result.numIntersections++] = ccResult.point[i];
            result.intersect = true;
        }
    }
    return result;
}


}
