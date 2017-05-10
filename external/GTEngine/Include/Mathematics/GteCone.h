// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <Mathematics/GteConstants.h>
#include <Mathematics/GteRay.h>
#include <limits>

// An acute cone is Dot(A,X-V) = |X-V| cos(t) where V is the vertex, A is the
// unit-length direction of the axis of the cone, and T is the cone angle with
// 0 < t < pi/2.  The cone interior is defined by the inequality
// Dot(A,X-V) >= |X-V| cos(t).  Since cos(t) > 0, we can avoid computing
// square roots.  The solid cone is defined by the inequality
// Dot(A,X-V)^2 >= Dot(X-V,X-V) cos(t)^2.  This is an infinite, single-sided
// cone.
//
// The cone may be truncated by a plane perpendicular to its axis at a height
// h from the vertex (distance from the vertex to the intersection of the
// plane and the axis).  The infinite cone has h = infinity.  The finite cone
// has a disk of intersection between the plane and infinite cone.  The radius
// r of the disk is r = h*tan(t).

namespace gte
{

template <int N, typename Real>
class Cone
{
public:
    // Construction and destruction.  The default constructor sets center to
    // (0,...,0), axis to (0,...,0,1), cosAngle and sinAngle to 1/sqrt(2)
    // [angle is pi/4], and height to 1.
    Cone();

    // The axis direction must be unit-length and the angle must be in
    // (0,pi/2).  The height is set to std::numeric_limits<float>::max().
    Cone(Ray<N, Real> const& inRay, Real inAngle);

    // The axis direction must be unit-length and the angle must be in
    // (0,pi/2).  The height must be positive.
    // std::numeric_limits<float>::max().
    Cone(Ray<N, Real> const& inRay, Real inAngle, Real inHeight);

    // The angle must be in (0,pi/2).  The function sets 'angle' and computes
    // 'cosAngle' and 'sinAngle'.
    void SetAngle(Real inAngle);

    // The cone vertex is the ray origin and the cone axis direction is the
    // ray direction.  The direction must be unit length.  The angle must be
    // in (0,pi/2).  The height must be in (0,+infinity), where +infinity is
    // std::numeric_limits<Real>::max().
    Ray<N, Real> ray;
    Real angle;
    Real height;

    // Members derived from 'angle', to avoid calling trigonometric functions
    // in geometric queries (for speed).  You may set 'angle' and compute
    // these by calling SetAngle(inAngle).
    Real cosAngle, sinAngle, cosAngleSqr;

public:
    // Comparisons to support sorted containers.  These based only on 'ray',
    // 'angle', and 'height'.
    bool operator==(Cone const& cone) const;
    bool operator!=(Cone const& cone) const;
    bool operator< (Cone const& cone) const;
    bool operator<=(Cone const& cone) const;
    bool operator> (Cone const& cone) const;
    bool operator>=(Cone const& cone) const;
};

// Template alias for convenience.
template <typename Real>
using Cone3 = Cone<3, Real>;


template <int N, typename Real>
Cone<N, Real>::Cone()
    :
    angle((Real)GTE_C_QUARTER_PI),
    height((Real)1),
    cosAngle(angle),
    sinAngle(angle),
    cosAngleSqr(cosAngle * cosAngle)
{
    ray.origin.MakeZero();
    ray.direction.MakeUnit(N - 1);
}

template <int N, typename Real>
Cone<N, Real>::Cone(Ray<N, Real> const& inRay, Real inAngle)
    :
    ray(inRay),
    angle(inAngle),
    height(std::numeric_limits<Real>::max()),
    cosAngle(cos(angle)),
    sinAngle(sin(angle)),
    cosAngleSqr(cosAngle * cosAngle)
{
}

template <int N, typename Real>
Cone<N, Real>::Cone(Ray<N, Real> const& inRay, Real inAngle, Real inHeight)
    :
    ray(inRay),
    angle(inAngle),
    height(inHeight),
    cosAngle(cos(angle)),
    sinAngle(sin(angle)),
    cosAngleSqr(cosAngle * cosAngle)
{
}

template <int N, typename Real>
void Cone<N, Real>::SetAngle(Real inAngle)
{
    angle = inAngle;
    cosAngle = cos(angle);
    sinAngle = sin(angle);
    cosAngleSqr = cosAngle * cosAngle;
}

template <int N, typename Real>
bool Cone<N, Real>::operator==(Cone const& cone) const
{
    return ray == cone.ray && angle == cone.angle && height == cone.height;
}

template <int N, typename Real>
bool Cone<N, Real>::operator!=(Cone const& cone) const
{
    return !operator==(cone);
}

template <int N, typename Real>
bool Cone<N, Real>::operator<(Cone const& cone) const
{
    if (ray < cone.ray)
    {
        return true;
    }

    if (ray > cone.ray)
    {
        return false;
    }

    if (angle < cone.angle)
    {
        return true;
    }

    if (angle > cone.angle)
    {
        return false;
    }

    return height < cone.height;
}

template <int N, typename Real>
bool Cone<N, Real>::operator<=(Cone const& cone) const
{
    return operator<(cone) || operator==(cone);
}

template <int N, typename Real>
bool Cone<N, Real>::operator>(Cone const& cone) const
{
    return !operator<=(cone);
}

template <int N, typename Real>
bool Cone<N, Real>::operator>=(Cone const& cone) const
{
    return !operator<(cone);
}


}
