// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.1 (2017/05/23)

#pragma once

#include <Mathematics/GteOrientedBox.h>
#include <Mathematics/GteIntrSegment3AlignedBox3.h>

// The queries consider the box to be a solid.
//
// The test-intersection queries use the method of separating axes.  The
// find-intersection queries use parametric clipping against the six faces of
// the box.

namespace gte
{

template <typename Real>
class TIQuery<Real, Segment3<Real>, OrientedBox3<Real>>
    :
    public TIQuery<Real, Segment3<Real>, AlignedBox3<Real>>
{
public:
    struct Result
        :
        public TIQuery<Real, Segment3<Real>, AlignedBox3<Real>>::Result
    {
        // No additional information to compute.
    };

    Result operator()(Segment3<Real> const& segment,
        OrientedBox3<Real> const& box);
};

template <typename Real>
class FIQuery<Real, Segment3<Real>, OrientedBox3<Real>>
    :
    public FIQuery<Real, Segment3<Real>, AlignedBox3<Real>>
{
public:
    struct Result
        :
        public FIQuery<Real, Segment3<Real>, AlignedBox3<Real>>::Result
    {
        // No additional relevant information to compute.
    };

    Result operator()(Segment3<Real> const& segment,
        OrientedBox3<Real> const& box);
};


template <typename Real>
typename TIQuery<Real, Segment3<Real>, OrientedBox3<Real>>::Result
TIQuery<Real, Segment3<Real>, OrientedBox3<Real>>::operator()(
    Segment3<Real> const& segment, OrientedBox3<Real> const& box)
{
    // Transform the segment to the oriented-box coordinate system.
    Vector3<Real> tmpOrigin, tmpDirection;
    Real segExtent;
    segment.GetCenteredForm(tmpOrigin, tmpDirection, segExtent);
    Vector3<Real> diff = tmpOrigin - box.center;
    Vector3<Real> segOrigin
    {
        Dot(diff, box.axis[0]),
        Dot(diff, box.axis[1]),
        Dot(diff, box.axis[2])
    };
    Vector3<Real> segDirection
    {
        Dot(tmpDirection, box.axis[0]),
        Dot(tmpDirection, box.axis[1]),
        Dot(tmpDirection, box.axis[2])
    };

    Result result;
    this->DoQuery(segOrigin, segDirection, segExtent, box.extent, result);
    return result;
}

template <typename Real>
typename FIQuery<Real, Segment3<Real>, OrientedBox3<Real>>::Result
FIQuery<Real, Segment3<Real>, OrientedBox3<Real>>::operator()(
    Segment3<Real> const& segment, OrientedBox3<Real> const& box)
{
    // Transform the segment to the oriented-box coordinate system.
    Vector3<Real> tmpOrigin, tmpDirection;
    Real segExtent;
    segment.GetCenteredForm(tmpOrigin, tmpDirection, segExtent);
    Vector3<Real> diff = tmpOrigin - box.center;
    Vector3<Real> segOrigin
    {
        Dot(diff, box.axis[0]),
        Dot(diff, box.axis[1]),
        Dot(diff, box.axis[2])
    };
    Vector3<Real> segDirection
    {
        Dot(tmpDirection, box.axis[0]),
        Dot(tmpDirection, box.axis[1]),
        Dot(tmpDirection, box.axis[2])
    };

    Result result;
    this->DoQuery(segOrigin, segDirection, segExtent, box.extent, result);
    for (int i = 0; i < result.numPoints; ++i)
    {
        // Compute the intersection point in the oriented-box coordinate
        // system.
        Vector3<Real> c = segOrigin + result.lineParameter[i] * segDirection;

        // Transform the intersection point to the original coordinate
        // system.
        result.point[i] = box.center;
        for (int j = 0; j < 3; ++j)
        {
            result.point[i] += c[j] * box.axis[j];
        }
    }

    return result;
}


}
