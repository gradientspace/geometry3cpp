// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <Mathematics/GteVector3.h>
#include <Mathematics/GteHypersphere.h>
#include <Mathematics/GteCone.h>
#include <Mathematics/GteFIQuery.h>
#include <Mathematics/GteTIQuery.h>

namespace gte
{

template <typename Real>
class TIQuery<Real, Sphere3<Real>, Cone3<Real>>
{
public:
    struct Result
    {
        bool intersect;
    };

    Result operator()(Sphere3<Real> const& sphere, Cone3<Real> const& cone);
};

template <typename Real>
class FIQuery<Real, Sphere3<Real>, Cone3<Real>>
{
public:
    struct Result
    {
        // If an intersection occurs, it is potentially an infinite set.  If
        // the cone vertex is inside the sphere, 'point' is set to the cone
        // vertex; else if the sphere center is inside the cone, 'point' is
        // set to the sphere center; else 'point' is set to the cone point
        // closest to the cone vertex.
        bool intersect;
        Vector3<Real> point;
    };

    Result operator()(Sphere3<Real> const& sphere, Cone3<Real> const& cone);
};


template <typename Real>
typename TIQuery<Real, Sphere3<Real>, Cone3<Real>>::Result
TIQuery<Real, Sphere3<Real>, Cone3<Real>>::operator()(
    Sphere3<Real> const& sphere, Cone3<Real> const& cone)
{
    Result result;

    Real invSin = ((Real)1) / cone.sinAngle;
    Real cosSqr = cone.cosAngle * cone.cosAngle;
    Vector3<Real> CmV = sphere.center - cone.ray.origin;
    Vector3<Real> D = CmV + (sphere.radius * invSin) * cone.ray.direction;
    Real lenSqr = Dot(D, D);
    Real e = Dot(D, cone.ray.direction);
    if (e > (Real)0 && e*e >= lenSqr*cosSqr)
    {
        Real sinSqr = cone.sinAngle * cone.sinAngle;
        lenSqr = Dot(CmV, CmV);
        e = -Dot(CmV, cone.ray.direction);
        if (e > (Real)0 && e*e >= lenSqr*sinSqr)
        {
            Real rSqr = sphere.radius * sphere.radius;
            result.intersect = (lenSqr <= rSqr);
        }
        else
        {
            result.intersect = true;
        }
    }
    else
    {
        result.intersect = false;
    }

    return result;
}

template <typename Real>
typename FIQuery<Real, Sphere3<Real>, Cone3<Real>>::Result
FIQuery<Real, Sphere3<Real>, Cone3<Real>>::operator()(
    Sphere3<Real> const& sphere, Cone3<Real> const& cone)
{
    Result result;

    // Test whether the cone vertex is inside the sphere.
    Vector3<Real> diff = sphere.center - cone.ray.origin;
    Real rSqr = sphere.radius * sphere.radius;
    Real lenSqr = Dot(diff, diff);
    if (lenSqr <= rSqr)
    {
        // The cone vertex is inside the sphere, so the sphere and cone
        // intersect.
        result.intersect = true;
        result.point = cone.ray.origin;
        return result;
    }

    // Test whether the sphere center is inside the cone.
    Real dot = Dot(diff, cone.ray.direction);
    Real dotSqr = dot*dot;
    Real cosSqr = cone.cosAngle * cone.cosAngle;
    if (dotSqr >= lenSqr*cosSqr && dot > (Real)0)
    {
        // The sphere center is inside cone, so the sphere and cone
        // intersect.
        result.intersect = true;
        result.point = sphere.center;
        return result;
    }

    // The sphere center is outside the cone.  The problem now reduces to
    // computing an intersection between the circle and the ray in the plane
    // containing the cone vertex and spanned by the cone axis and vector
    // from the cone vertex to the sphere center.

    // The ray is t*D+V (t >= 0) where |D| = 1 and dot(A,D) = cos(angle).
    // Also, D = e*A+f*(C-V).  Plugging the ray equation into sphere equation
    // yields R^2 = |t*D+V-C|^2, so the quadratic for intersections is
    // t^2 - 2*dot(D,C-V)*t + |C-V|^2 - R^2 = 0.  An intersection occurs
    // if and only if the discriminant is nonnegative.  This test becomes
    //
    //     dot(D,C-V)^2 >= dot(C-V,C-V) - R^2
    //
    // Note that if the right-hand side is nonpositive, then the inequality
    // is true (the sphere contains V).  This is already ruled out in the
    // first block of code in this function.

    Real uLen = sqrt(std::max(lenSqr - dotSqr, (Real)0));
    Real test = cone.cosAngle * dot + cone.sinAngle * uLen;
    Real discr = test * test - lenSqr + rSqr;

    if (discr >= (Real)0 && test >= (Real)0)
    {
        // Compute the point of intersection closest to the cone vertex.
        result.intersect = true;
        Real t = test - sqrt(std::max(discr, (Real)0));
        Vector3<Real> B = diff - dot * cone.ray.direction;
        Real tmp = cone.sinAngle / uLen;
        result.point = t * (cone.cosAngle * cone.ray.direction + tmp * B);
    }
    else
    {
        result.intersect = false;
    }

    return result;
}


}
