// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <Mathematics/GteVector2.h>
#include <Mathematics/GteOrientedBox.h>
#include <Mathematics/GteFIQuery.h>
#include <Mathematics/GteTIQuery.h>

// The queries consider the box to be a solid.
//
// The test-intersection query uses the method of separating axes.  The set of
// potential separating directions includes the 2 edge normals of box0 and the
// 2 edge normals of box1.  The integer 'separating' identifies the axis that
// reported separation; there may be more than one but only one is reported.
// The value is 0 when box0.axis[0] separates, 1 when box0.axis[1] separates,
// 2 when box1.axis[0] separates, or 3 when box1.axis[1] separates.

namespace gte
{

template <typename Real>
class TIQuery<Real, OrientedBox2<Real>, OrientedBox2<Real>>
{
public:
    struct Result
    {
        bool intersect;
        int separating;
    };

    Result operator()(OrientedBox2<Real> const& box0,
        OrientedBox2<Real> const& box1);
};


template <typename Real>
typename TIQuery<Real, OrientedBox2<Real>, OrientedBox2<Real>>::Result
TIQuery<Real, OrientedBox2<Real>, OrientedBox2<Real>>::operator()(
    OrientedBox2<Real> const& box0, OrientedBox2<Real> const& box1)
{
    Result result;

    // Convenience variables.
    Vector2<Real> const* A0 = &box0.axis[0];
    Vector2<Real> const* A1 = &box1.axis[0];
    Vector2<Real> const& E0 = box0.extent;
    Vector2<Real> const& E1 = box1.extent;

    // Compute difference of box centers, D = C1-C0.
    Vector2<Real> D = box1.center - box0.center;

    Real absA0dA1[2][2], rSum;

    // Test box0.axis[0].
    absA0dA1[0][0] = std::abs(Dot(A0[0], A1[0]));
    absA0dA1[0][1] = std::abs(Dot(A0[0], A1[1]));
    rSum = E0[0] + E1[0] * absA0dA1[0][0] + E1[1] * absA0dA1[0][1];
    if (std::abs(Dot(A0[0], D)) > rSum)
    {
        result.intersect = false;
        result.separating = 0;
        return result;
    }

    // Test axis box0.axis[1].
    absA0dA1[1][0] = std::abs(Dot(A0[1], A1[0]));
    absA0dA1[1][1] = std::abs(Dot(A0[1], A1[1]));
    rSum = E0[1] + E1[0] * absA0dA1[1][0] + E1[1] * absA0dA1[1][1];
    if (std::abs(Dot(A0[1], D)) > rSum)
    {
        result.intersect = false;
        result.separating = 1;
        return result;
    }

    // Test axis box1.axis[0].
    rSum = E1[0] + E0[0] * absA0dA1[0][0] + E0[1] * absA0dA1[1][0];
    if (std::abs(Dot(A1[0], D)) > rSum)
    {
        result.intersect = false;
        result.separating = 2;
        return result;
    }

    // Test axis box1.axis[1].
    rSum = E1[1] + E0[0] * absA0dA1[0][1] + E0[1] * absA0dA1[1][1];
    if (std::abs(Dot(A1[1], D)) > rSum)
    {
        result.intersect = false;
        result.separating = 3;
        return result;
    }

    result.intersect = true;
    return result;
}


}
