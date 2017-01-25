// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <Mathematics/GteArc2.h>
#include <Mathematics/GteIntrSegment2Circle2.h>

// The queries consider the arc to be a 1-dimensional object.

namespace gte
{

template <typename Real>
class TIQuery<Real, Segment2<Real>, Arc2<Real>>
{
public:
    struct Result
    {
        bool intersect;
    };

    Result operator()(Segment2<Real> const& segment, Arc2<Real> const& arc);
};

template <typename Real>
class FIQuery<Real, Segment2<Real>, Arc2<Real>>
{
public:
    struct Result
    {
        bool intersect;
        int numIntersections;
        std::array<Real, 2> parameter;
        std::array<Vector2<Real>, 2> point;
    };

    Result operator()(Segment2<Real> const& segment, Arc2<Real> const& arc);
};


template <typename Real>
typename TIQuery<Real, Segment2<Real>, Arc2<Real>>::Result
TIQuery<Real, Segment2<Real>, Arc2<Real>>::operator()(
    Segment2<Real> const& segment, Arc2<Real> const& arc)
{
    Result result;
    FIQuery<Real, Segment2<Real>, Arc2<Real>> saQuery;
    auto saResult = saQuery(segment, arc);
    result.intersect = saResult.intersect;
    return result;
}

template <typename Real>
typename FIQuery<Real, Segment2<Real>, Arc2<Real>>::Result
FIQuery<Real, Segment2<Real>, Arc2<Real>>::operator()(
    Segment2<Real> const& segment, Arc2<Real> const& arc)
{
    Result result;

    FIQuery<Real, Segment2<Real>, Circle2<Real>> scQuery;
    Circle2<Real> circle(arc.center, arc.radius);
    auto scResult = scQuery(segment, circle);
    if (scResult.intersect)
    {
        // Test whether line-circle intersections are on the arc.
        result.numIntersections = 0;
        for (int i = 0; i < scResult.numIntersections; ++i)
        {
            if (arc.Contains(scResult.point[i]))
            {
                result.intersect = true;
                result.parameter[result.numIntersections]
                    = scResult.parameter[i];
                result.point[result.numIntersections++]
                    = scResult.point[i];
            }
        }
    }
    else
    {
        result.intersect = false;
        result.numIntersections = 0;
    }

    return result;
}


}
