// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <Mathematics/GteVector3.h>

// The torus has center C with plane of symmetry containing C and having
// directions D0 and D1.  The axis of symmetry is the line containing C
// and having direction N (the plane normal).  The radius from the center
// of the torus is r0 (outer radius) and the radius of the tube of the
// torus is r1 (inner radius).  It must be that r0 >= r1.  A point X may
// be written as X = C + y0*D0 + y1*D1 + y2*N, where matrix [U V N] is
// orthonormal and has determinant 1.  Thus, y0 = Dot(D0,X-C),
// y1 = Dot(D1,X-C), and y2 = Dot(N,X-C).  The implicit form is
//      [|X-C|^2 - (r0^2 + r1^2)]^2 + 4*r0^2*((Dot(N,X-C))^2 - r1^2) = 0
// Note that D0 and D1 are not present in the equation, which is to be
// expected by the symmetry.  The parametric form is
//      X(u,v) = (r0 + r1*cos(v))*(cos(u)*D0 + sin(u)*D1) + r1*sin(v)*N
// for -pi <= u < pi, -pi <= v < pi.  The member 'center' is C, 'direction0'
// is D0, 'direction1' is D1, 'normal' is N, 'radius0' is r0, and 'radius1'
// is r1.

namespace gte
{

template <typename Real>
class Torus3
{
public:
    // Construction and destruction.  The default constructor sets center to
    // (0,0,0), direction0 to (1,0,0), direction1 to (0,1,0), normal to
    // (0,0,1), radius0 to 2, and radius1 to 1.
    Torus3();
    Torus3(Vector3<Real> const& inCenter, Vector3<Real> const& inDirection0,
        Vector3<Real> const& inDirection1, Vector3<Real> const& inNormal,
        Real inRadius0, Real inRadius1);

    // Evaluation of the surface.  The function supports derivative
    // calculation through order 2; that is, maxOrder <= 2 is required.  If
    // you want only the position, pass in maxOrder of 0.  If you want the
    // position and first-order derivatives, pass in maxOrder of 1, and so on.
    // The output 'values' are ordered as: position X; first-order derivatives
    // dX/du, dX/dv; second-order derivatives d2X/du2, d2X/dudv, d2X/dv2.
    void Evaluate(Real u, Real v, unsigned int maxOrder,
        Vector3<Real> values[6]) const;

    // Reverse lookup of parameters from position.
    void GetParameters(Vector3<Real> const& X, Real& u, Real& v) const;

    Vector3<float> center, direction0, direction1, normal;
    Real radius0, radius1;

public:
    // Comparisons to support sorted containers.
    bool operator==(Torus3 const& torus) const;
    bool operator!=(Torus3 const& torus) const;
    bool operator< (Torus3 const& torus) const;
    bool operator<=(Torus3 const& torus) const;
    bool operator> (Torus3 const& torus) const;
    bool operator>=(Torus3 const& torus) const;
};


template <typename Real>
Torus3<Real>::Torus3()
    :
    center(Vector3<Real>::Zero()),
    direction0(Vector3<Real>::Unit(0)),
    direction1(Vector3<Real>::Unit(1)),
    normal(Vector3<Real>::Unit(2)),
    radius0((Real)2),
    radius1((Real)1)
{
}

template <typename Real>
Torus3<Real>::Torus3(Vector3<Real> const& inCenter,
    Vector3<Real> const& inDirection0, Vector3<Real> const& inDirection1,
    Vector3<Real> const& inNormal, Real inRadius0, Real inRadius1)
    :
    center(inCenter),
    direction0(inDirection0),
    direction1(inDirection1),
    normal(inNormal),
    radius0(inRadius0),
    radius1(inRadius1)
{
}

template <typename Real>
void Torus3<Real>::Evaluate(Real u, Real v, unsigned int maxOrder,
    Vector3<Real> values[6]) const
{
    // Compute position.
    Real csu = cos(u), snu = sin(u), csv = cos(v), snv = sin(v);
    Real r1csv = radius1 * csv;
    Real r1snv = radius1 * snv;
    Real r0pr1csv = radius0 + r1csv;
    Vector3<Real> combo0 = csu * direction0 + snu * direction1;
    Vector3<Real> r0pr1csvcombo0 = r0pr1csv * combo0;
    Vector3<Real> r1snvnormal = r1snv * normal;
    values[0] = center + r0pr1csvcombo0 + r1snvnormal;

    if (maxOrder >= 1)
    {
        // Compute first-order derivatives.
        Vector3<Real> combo1 = -snu * direction0 + csu * direction1;
        values[1] = r0pr1csv * combo1;
        values[2] = -r1snv * combo0 + r1csv * normal;

        if (maxOrder >= 2)
        {
            // Compute second-order derivatives.
            values[3] = -r0pr1csvcombo0;
            values[4] = -r1snv * combo1;
            values[5] = -r1csv * combo0 - r1snvnormal;

            if (maxOrder >= 3)
            {
                // These orders are not supported.
                for (int i = 0; i < 6; ++i)
                {
                    values[i] = Vector3<float>::Zero();
                }
            }
        }
    }
}

template <typename Real>
void Torus3<Real>::GetParameters(Vector3<Real> const& X, Real& u, Real& v)
const
{
    Vector3<Real> delta = X - center;
    Real dot0 = Dot(direction0, delta);  // (r0 + r1*cos(v))*cos(u)
    Real dot1 = Dot(direction1, delta);  // (r0 + r1*cos(v))*sin(u)
    Real dot2 = Dot(normal, delta);      // r1*sin(v)
    Real r1csv = sqrt(dot0 * dot0 + dot1 * dot1) - radius0;  // r1*cos(v)
    u = atan2(dot1, dot0);
    v = atan2(dot2, r1csv);
}

template <typename Real>
bool Torus3<Real>::operator==(Torus3 const& torus) const
{
    return center == torus.center
        && direction0 == torus.direction0
        && direction1 == torus.direction1
        && normal == torus.normal
        && radius0 == torus.radius0
        && radius1 == torus.radius1;
}

template <typename Real>
bool Torus3<Real>::operator!=(Torus3 const& torus) const
{
    return !operator==(torus);
}

template <typename Real>
bool Torus3<Real>::operator<(Torus3 const& torus) const
{
    if (center < torus.center)
    {
        return true;
    }

    if (center > torus.center)
    {
        return false;
    }

    if (direction0 < torus.direction0)
    {
        return true;
    }

    if (direction0 > torus.direction0)
    {
        return false;
    }

    if (direction1 < torus.direction1)
    {
        return true;
    }

    if (direction1 > torus.direction1)
    {
        return false;
    }

    if (normal < torus.normal)
    {
        return true;
    }

    if (normal > torus.normal)
    {
        return false;
    }

    if (radius0 < torus.radius0)
    {
        return true;
    }

    if (radius0 > torus.radius0)
    {
        return false;
    }

    return radius1 < torus.radius1;
}

template <typename Real>
bool Torus3<Real>::operator<=(Torus3 const& torus) const
{
    return operator<(torus) || operator==(torus);
}

template <typename Real>
bool Torus3<Real>::operator>(Torus3 const& torus) const
{
    return !operator<=(torus);
}

template <typename Real>
bool Torus3<Real>::operator>=(Torus3 const& torus) const
{
    return !operator<(torus);
}


}
