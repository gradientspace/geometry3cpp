// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <Mathematics/GteOrientedBox.h>
#include <Mathematics/GteIntrLine3AlignedBox3.h>

// The queries consider the box to be a solid.
//
// The test-intersection queries use the method of separating axes.  The
// find-intersection queries use parametric clipping against the six faces of
// the box.

namespace gte
{

template <typename Real>
class TIQuery<Real, Line3<Real>, OrientedBox3<Real>>
    :
    public TIQuery<Real, Line3<Real>, AlignedBox3<Real>>
{
public:
    struct Result
        :
        public TIQuery<Real, Line3<Real>, AlignedBox3<Real>>::Result
    {
        // No additional relevant information to compute.
    };

    Result operator()(Line3<Real> const& line, OrientedBox3<Real> const& box);
};

template <typename Real>
class FIQuery<Real, Line3<Real>, OrientedBox3<Real>>
    :
    public FIQuery<Real, Line3<Real>, AlignedBox3<Real>>
{
public:
    struct Result
        :
        public FIQuery<Real, Line3<Real>, AlignedBox3<Real>>::Result
    {
        // No additional relevant information to compute.
    };

    Result operator()(Line3<Real> const& line, OrientedBox3<Real> const& box);
};


template <typename Real>
typename TIQuery<Real, Line3<Real>, OrientedBox3<Real>>::Result
TIQuery<Real, Line3<Real>, OrientedBox3<Real>>::operator()(
    Line3<Real> const& line, OrientedBox3<Real> const& box)
{
    // Transform the line to the oriented-box coordinate system.
    Vector3<Real> diff = line.origin - box.center;
    Vector3<Real> lineOrigin
    {
        Dot(diff, box.axis[0]),
        Dot(diff, box.axis[1]),
        Dot(diff, box.axis[2])
    };
    Vector3<Real> lineDirection
    {
        Dot(line.direction, box.axis[0]),
        Dot(line.direction, box.axis[1]),
        Dot(line.direction, box.axis[2])
    };

    Result result;
    this->DoQuery(lineOrigin, lineDirection, box.extent, result);
    return result;
}

template <typename Real>
typename FIQuery<Real, Line3<Real>, OrientedBox3<Real>>::Result
FIQuery<Real, Line3<Real>, OrientedBox3<Real>>::operator()(
    Line3<Real> const& line, OrientedBox3<Real> const& box)
{
    // Transform the line to the oriented-box coordinate system.
    Vector3<Real> diff = line.origin - box.center;
    Vector3<Real> lineOrigin
    {
        Dot(diff, box.axis[0]),
        Dot(diff, box.axis[1]),
        Dot(diff, box.axis[2])
    };
    Vector3<Real> lineDirection
    {
        Dot(line.direction, box.axis[0]),
        Dot(line.direction, box.axis[1]),
        Dot(line.direction, box.axis[2])
    };

    Result result;
    this->DoQuery(lineOrigin, lineDirection, box.extent, result);
    for (int i = 0; i < result.numPoints; ++i)
    {
        result.point[i] =
            line.origin + result.lineParameter[i] * line.direction;
    }
    return result;
}


}
