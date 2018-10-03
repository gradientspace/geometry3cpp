// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.1 (2017/07/25)

#pragma once

#include <Mathematics/GteCone.h>
#include <Mathematics/GteOrientedBox.h>
#include <Mathematics/GteIntrAlignedBox3Cone3.h>

// Test for intersection of a box and a cone.  The cone can be finite or
// infinite.  The algorithm is described in
//   http://www.geometrictools.com/Documentation/IntersectionBoxCone.pdf
// and assumes that the intersection set must have positive volume.  For
// example, let the box be outside the cone.  If the box is below the support
// plane at the cone vertex and just touches the cone vertex, nointersection
// is reported.  If the box is above the plane of the disk capping a finite
// cone, no intersection is reported.  However, if the box straddles the
// support plane and just touches the cone vertex, an intersection is
// reported.  This is a consequence of wanting a fast test for culling boxes
// against a cone.  It is possible to add more logic to change the behavior.

namespace gte
{

template <typename Real>
class TIQuery<Real, OrientedBox<3, Real>, Cone<3, Real>>
    :
    public TIQuery<Real, AlignedBox<3, Real>, Cone<3, Real>>
{
public:
    struct Result
        :
        public TIQuery<Real, AlignedBox<3, Real>, Cone<3, Real>>::Result
    {
        // No additional information to compute.
    };

    Result operator()(OrientedBox<3, Real> const& box, Cone<3, Real> const& cone);
};

// Template alias for convenience.
template <typename Real>
using TIOrientedBox3Cone3 =
TIQuery<Real, OrientedBox<3, Real>, Cone<3, Real>>;


template <typename Real>
typename TIQuery<Real, OrientedBox<3, Real>, Cone<3, Real>>::Result
    TIQuery<Real, OrientedBox<3, Real>, Cone<3, Real>>::operator()(
    OrientedBox<3, Real> const& box, Cone<3, Real> const& cone)
{
    Result result;

    // Quick-rejection test for boxes below the supporting plane of the cone.
    Vector<3, Real> CmV = box.center - cone.ray.origin;
    Vector<3, Real> DdU{
        Dot(cone.ray.direction, box.axis[0]),
        Dot(cone.ray.direction, box.axis[1]),
        Dot(cone.ray.direction, box.axis[2]) };
    Real DdCmV = Dot(cone.ray.direction, CmV);  // interval center
    Real radius =  // interval half-length
        box.extent[0] * std::abs(DdU[0]) +
        box.extent[1] * std::abs(DdU[1]) +
        box.extent[2] * std::abs(DdU[2]);
    if (DdCmV + radius <= (Real)0)
    {
        // The box is in the halfspace below the supporting plane of the cone.
        result.intersect = false;
        return result;
    }

    // Quick-rejection test for boxes outside the plane determined by the
    // height of the cone.
    if (cone.height < std::numeric_limits<Real>::max())
    {
        if (DdCmV - radius >= cone.height)
        {
            // The box is outside the plane determined by the height of the
            // cone.
            result.intersect = false;
            return result;
        }
    }

    // Determine the box faces that are visible to the cone vertex.  The
    // box center has been translated (C-V) so that the cone vertex is at
    // the origin.  Compute the coordinates of the origin relative to the
    // translated box.
    Vector<3, Real> UdCmV{
        Dot(box.axis[0], CmV),
        Dot(box.axis[1], CmV),
        Dot(box.axis[2], CmV) };
    int index[3] = {
        (UdCmV[0] < -box.extent[0] ? 2 : (UdCmV[0] > box.extent[0] ? 0 : 1)),
        (UdCmV[1] < -box.extent[1] ? 2 : (UdCmV[1] > box.extent[1] ? 0 : 1)),
        (UdCmV[2] < -box.extent[2] ? 2 : (UdCmV[2] > box.extent[2] ? 0 : 1))
    };
    int lookup = index[0] + 3 * index[1] + 9 * index[2];
    if (lookup == 13)
    {
        // The cone vertex is in the box.
        result.intersect = true;
        return result;
    }

    auto const& polygon = this->mPolygon[lookup];

    this->DoQuery(box.extent, cone.cosAngleSqr, DdU, UdCmV, DdCmV, polygon,
        result);
    return result;
}


}
