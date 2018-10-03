// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <Mathematics/GteCylinder3.h>
#include <Mathematics/GteHalfspace.h>
#include <Mathematics/GteFIQuery.h>
#include <Mathematics/GteTIQuery.h>

// Queries for intersection of objects with halfspaces.  These are useful for
// containment testing, object culling, and clipping.

namespace gte
{

template <typename Real>
class TIQuery<Real, Halfspace3<Real>, Cylinder3<Real>>
{
public:
    struct Result
    {
        bool intersect;
    };

    Result operator()(Halfspace3<Real> const& halfspace,
        Cylinder3<Real> const& cylinder);
};


template <typename Real>
typename TIQuery<Real, Halfspace3<Real>, Cylinder3<Real>>::Result
TIQuery<Real, Halfspace3<Real>, Cylinder3<Real>>::operator()(
    Halfspace3<Real> const& halfspace, Cylinder3<Real> const& cylinder)
{
    Result result;

    // Compute extremes of signed distance Dot(N,X)-d for points on the
    // cylinder.  These are
    //   min = (Dot(N,C)-d) - r*sqrt(1-Dot(N,W)^2) - (h/2)*|Dot(N,W)|
    //   max = (Dot(N,C)-d) + r*sqrt(1-Dot(N,W)^2) + (h/2)*|Dot(N,W)|
    Real center = Dot(halfspace.normal, cylinder.axis.origin) -
        halfspace.constant;
    Real absNdW = std::abs(Dot(halfspace.normal, cylinder.axis.direction));
    Real root = sqrt(std::max((Real)1, (Real)1 - absNdW * absNdW));
    Real tmax = center + cylinder.radius*root +
        ((Real)0.5)*cylinder.height*absNdW;

    // The cylinder and halfspace intersect when the projection interval
    // maximum is nonnegative.
    result.intersect = (tmax >= (Real)0);
    return result;
}


}
