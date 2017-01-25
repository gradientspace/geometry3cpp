// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <Mathematics/GteVector2.h>
#include <Mathematics/GteSegment.h>
#include <Mathematics/GteIntrIntervals.h>
#include <Mathematics/GteIntrLine2AlignedBox2.h>

// The queries consider the box to be a solid.
//
// The test-intersection queries use the method of separating axes.  The
// find-intersection queries use parametric clipping against the four edges of
// the box.

namespace gte
{


template <typename Real>
class TIQuery<Real, Segment2<Real>, AlignedBox2<Real>>
    :
    public TIQuery<Real, Line2<Real>, AlignedBox2<Real>>
{
public:
    struct Result
        :
        public TIQuery<Real, Line2<Real>, AlignedBox2<Real>>::Result
    {
        // No additional information to compute.
    };

    Result operator()(Segment2<Real> const& segment,
        AlignedBox2<Real> const& box);

protected:
    void DoQuery(Vector2<Real> const& segOrigin,
        Vector2<Real> const& segDirection, Real segExtent,
        Vector2<Real> const& boxExtent, Result& result);
};

template <typename Real>
class FIQuery<Real, Segment2<Real>, AlignedBox2<Real>>
    :
    public FIQuery<Real, Line2<Real>, AlignedBox2<Real>>
{
public:
    struct Result
        :
        public FIQuery<Real, Line2<Real>, AlignedBox2<Real>>::Result
    {
        // No additional information to compute.
    };

    Result operator()(Segment2<Real> const& segment,
        AlignedBox2<Real> const& box);

protected:
    void DoQuery(Vector2<Real> const& segOrigin,
        Vector2<Real> const& segDirection, Real segExtent,
        Vector2<Real> const& boxExtent, Result& result);
};


template <typename Real>
typename TIQuery<Real, Segment2<Real>, AlignedBox2<Real>>::Result
TIQuery<Real, Segment2<Real>, AlignedBox2<Real>>::operator()(
    Segment2<Real> const& segment, AlignedBox2<Real> const& box)
{
    // Get the centered form of the aligned box.  The axes are implicitly
    // Axis[d] = Vector2<Real>::Unit(d).
    Vector2<Real> boxCenter, boxExtent;
    box.GetCenteredForm(boxCenter, boxExtent);

    // Transform the segment to a centered form in the aligned-box coordinate
    // system.
    Vector2<Real> transformedP0 = segment.p[0] - boxCenter;
    Vector2<Real> transformedP1 = segment.p[1] - boxCenter;
    Segment2<Real> transformedSegment(transformedP0, transformedP1);
    Vector2<Real> segOrigin, segDirection;
    Real segExtent;
    transformedSegment.GetCenteredForm(segOrigin, segDirection, segExtent);

    Result result;
    DoQuery(segOrigin, segDirection, segExtent, boxExtent, result);
    return result;
}

template <typename Real>
void TIQuery<Real, Segment2<Real>, AlignedBox2<Real>>::DoQuery(
    Vector2<Real> const& segOrigin, Vector2<Real> const& segDirection,
    Real segExtent, Vector2<Real> const& boxExtent, Result& result)
{
    for (int i = 0; i < 2; ++i)
    {
        if (std::abs(segOrigin[i]) > boxExtent[i] +
            segExtent*std::abs(segDirection[i]))
        {
            result.intersect = false;
            return;
        }
    }

    TIQuery<Real, Line2<Real>, AlignedBox2<Real>>::DoQuery(segOrigin,
        segDirection, boxExtent, result);
}

template <typename Real>
typename FIQuery<Real, Segment2<Real>, AlignedBox2<Real>>::Result
FIQuery<Real, Segment2<Real>, AlignedBox2<Real>>::operator()(
    Segment2<Real> const& segment, AlignedBox2<Real> const& box)
{
    // Get the centered form of the aligned box.  The axes are implicitly
    // Axis[d] = Vector2<Real>::Unit(d).
    Vector2<Real> boxCenter, boxExtent;
    box.GetCenteredForm(boxCenter, boxExtent);

    // Transform the segment to a centered form in the aligned-box coordinate
    // system.
    Vector2<Real> transformedP0 = segment.p[0] - boxCenter;
    Vector2<Real> transformedP1 = segment.p[1] - boxCenter;
    Segment2<Real> transformedSegment(transformedP0, transformedP1);
    Vector2<Real> segOrigin, segDirection;
    Real segExtent;
    transformedSegment.GetCenteredForm(segOrigin, segDirection, segExtent);

    Result result;
    DoQuery(segOrigin, segDirection, segExtent, boxExtent, result);
    for (int i = 0; i < result.numIntersections; ++i)
    {
        result.point[i] = segOrigin + result.parameter[i] * segDirection;
    }
    return result;
}

template <typename Real>
void FIQuery<Real, Segment2<Real>, AlignedBox2<Real>>::DoQuery(
    Vector2<Real> const& segOrigin, Vector2<Real> const& segDirection,
    Real segExtent, Vector2<Real> const& boxExtent, Result& result)
{
    FIQuery<Real, Line2<Real>, AlignedBox2<Real>>::DoQuery(segOrigin,
        segDirection, boxExtent, result);

    if (result.intersect)
    {
        // The line containing the segment intersects the box; the t-interval
        // is [t0,t1].  The segment intersects the box as long as [t0,t1]
        // overlaps the segment t-interval [-segExtent,+segExtent].
        std::array<Real, 2> segInterval = { -segExtent, segExtent };
        FIQuery<Real, std::array<Real, 2>, std::array<Real, 2>> iiQuery;
        auto iiResult = iiQuery(result.parameter, segInterval);
        result.intersect = iiResult.intersect;
        result.numIntersections = iiResult.numIntersections;
        result.parameter = iiResult.overlap;
    }
}


}
