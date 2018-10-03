// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <Mathematics/GteIntrAlignedBox3Sphere3.h>
#include <Mathematics/GteOrientedBox.h>

namespace gte
{

template <typename Real>
class TIQuery<Real, OrientedBox3<Real>, Sphere3<Real>>
    :
    public TIQuery<Real, AlignedBox3<Real>, Sphere3<Real>>
{
public:
    // The intersection query considers the box and sphere to be solids.
    // For example, if the sphere is strictly inside the box (does not touch
    // the box faces), the objects intersect.
    struct Result
        :
        public TIQuery<Real, AlignedBox3<Real>, Sphere3<Real>>::Result
    {
        // No additional information to compute.
    };

    Result operator()(OrientedBox3<Real> const& box, Sphere3<Real> const& sphere);
};

template <typename Real>
class FIQuery<Real, OrientedBox3<Real>, Sphere3<Real>>
    :
    public FIQuery<Real, AlignedBox3<Real>, Sphere3<Real>>
{
public:
    // Currently, only a dynamic query is supported.  The static query must
    // compute the intersection set of (solid) box and sphere.
    struct Result
        :
        public FIQuery<Real, AlignedBox3<Real>, Sphere3<Real>>::Result
    {
        // No additional information to compute.
    };

    Result operator()(Real maxTime, OrientedBox3<Real> const& box,
        Vector3<Real> const& boxVelocity, Sphere3<Real> const& sphere,
        Vector3<Real> const& sphereVelocity);
};


template <typename Real>
typename TIQuery<Real, OrientedBox3<Real>, Sphere3<Real>>::Result
TIQuery<Real, OrientedBox3<Real>, Sphere3<Real>>::operator()(
    OrientedBox3<Real> const& box, Sphere3<Real> const& sphere)
{
    // Test for intersection in the coordinate system of the box by
    // transforming the sphere into that coordinate system.
    Vector3<Real> temp = sphere.center - box.center;
    Vector3<Real> cdiff{ Dot(temp, box.axis[0]), Dot(temp, box.axis[1]), Dot(temp, box.axis[2]) };

    Result result;
    this->DoQuery(box.extent, cdiff, sphere.radius, result);
    return result;
}


template <typename Real>
typename FIQuery<Real, OrientedBox3<Real>, Sphere3<Real>>::Result
FIQuery<Real, OrientedBox3<Real>, Sphere3<Real>>::operator()(Real maxTime,
    OrientedBox3<Real> const& box, Vector3<Real> const& boxVelocity,
    Sphere3<Real> const& sphere, Vector3<Real> const& sphereVelocity)
{
    // Find intersections relative to the coordinate system of the box.
    // The sphere is transformed to the box coordinates and the velocity of
    // the sphere is relative to the box.
    Vector3<Real> temp = sphere.center - box.center;
    Vector3<Real> cdiff{ Dot(temp, box.axis[0]), Dot(temp, box.axis[1]), Dot(temp, box.axis[2]) };
    temp = sphereVelocity - boxVelocity;
    Vector3<Real> relvel{ Dot(temp, box.axis[0]), Dot(temp, box.axis[1]), Dot(temp, box.axis[2]) };

    Result result;
    this->DoQuery(maxTime, box.center, box.extent, sphere, cdiff, relvel, result);
    return result;
}

}
