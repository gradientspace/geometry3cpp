// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <Mathematics/GteVector2.h>
#include <Mathematics/GteDistPointOrientedBox.h>
#include <Mathematics/GteHypersphere.h>
#include <Mathematics/GteFIQuery.h>
#include <Mathematics/GteTIQuery.h>

namespace gte
{

template <typename Real>
class TIQuery<Real, OrientedBox2<Real>, Circle2<Real>>
{
public:
    // The intersection query considers the box and circle to be solids.
    // For example, if the circle is strictly inside the box (does not touch
    // the box edges), the objects intersect.
    struct Result
    {
        bool intersect;
    };

    Result operator()(OrientedBox2<Real> const& box,
        Circle2<Real> const& circle);
};

template <typename Real>
class FIQuery<Real, OrientedBox2<Real>, Circle2<Real>>
{
public:
    // Currently, only a dynamic query is supported.  The static query must
    // compute the intersection set of (solid) box and circle.
    struct Result
    {
        bool intersect;
        Real contactTime;
        Vector2<Real> contactPoint;
    };

    Result operator()(Real maxTime, OrientedBox2<Real> const& box,
        Vector2<Real> const& boxVelocity, Circle2<Real> const& circle,
        Vector2<Real> const& circleVelocity);

private:
    // Support for dynamic query.  Both functions return -1 if the objects are
    // initially intersecting, 0 if no intersection, or +1 if they intersect
    // at some positive time.
    int TestVertexRegion(Real cx, Real cy, Real vx, Real vy, Real ex,
        Real ey, Real& ix, Real& iy, Real radius, Real& contactTime);

    int TestEdgeRegion(Real cx, Real cy, Real vx, Real vy, Real ex, Real ey,
        Real& ix, Real& iy, Real radius, Real& contactTime);
};


template <typename Real>
typename TIQuery<Real, OrientedBox2<Real>, Circle2<Real>>::Result
TIQuery<Real, OrientedBox2<Real>, Circle2<Real>>::operator()(
    OrientedBox2<Real> const& box, Circle2<Real> const& circle)
{
    DCPQuery<Real, Vector2<Real>, OrientedBox2<Real>> pbQuery;
    auto pbResult = pbQuery(circle.center, box);
    Result result;
    result.intersect = (pbResult.distance <= circle.radius);
    return result;
}



template <typename Real>
typename FIQuery<Real, OrientedBox2<Real>, Circle2<Real>>::Result
FIQuery<Real, OrientedBox2<Real>, Circle2<Real>>::operator()(Real maxTime,
    OrientedBox2<Real> const& box, Vector2<Real> const& boxVelocity,
    Circle2<Real> const& circle, Vector2<Real> const& circleVelocity)
{
    Result result;

    // Convert circle center to box coordinates.
    Vector2<Real> diff = circle.center - box.center;
    Vector2<Real> relativeVelocity = circleVelocity - boxVelocity;
    Real cx = Dot(diff, box.axis[0]);
    Real cy = Dot(diff, box.axis[1]);
    Real vx = Dot(relativeVelocity, box.axis[0]);
    Real vy = Dot(relativeVelocity, box.axis[1]);
    Real ex = box.extent[0];
    Real ey = box.extent[1];
    Real ix, iy;

    int type = 0;

    if (cx < -ex)
    {
        if (cy < -ey)
        {
            // region Rmm
            type = TestVertexRegion(cx, cy, vx, vy, ex, ey, ix, iy,
                circle.radius, result.contactTime);
        }
        else if (cy <= ey)
        {
            // region Rmz
            type = TestEdgeRegion(cx, cy, vx, vy, ex, ey, ix, iy,
                circle.radius, result.contactTime);
        }
        else
        {
            // region Rmp
            type = TestVertexRegion(cx, -cy, vx, -vy, ex, ey, ix, iy,
                circle.radius, result.contactTime);
            iy = -iy;
        }
    }
    else if (cx <= ex)
    {
        if (cy < -ey)
        {
            // region Rzm
            type = TestEdgeRegion(cy, cx, vy, vx, ey, ex, iy, ix,
                circle.radius, result.contactTime);
        }
        else if (cy <= ey)
        {
            // region Rzz: The circle is already intersecting the box.  Use
            // it as the intersection point.
            result.intersect = true;
            result.contactTime = (Real)0;
            result.contactPoint = circle.center;
            return result;
        }
        else
        {
            // region Rzp
            type = TestEdgeRegion(-cy, cx, -vy, vx, ey, ex, iy, ix,
                circle.radius, result.contactTime);
            iy = -iy;
        }
    }
    else
    {
        if (cy < -ey)
        {
            // region Rpm
            type = TestVertexRegion(-cx, cy, -vx, vy, ex, ey, ix, iy,
                circle.radius, result.contactTime);
            ix = -ix;
        }
        else if (cy <= ey)
        {
            // region Rpz
            type = TestEdgeRegion(-cx, cy, -vx, vy, ex, ey, ix, iy,
                circle.radius, result.contactTime);
            ix = -ix;
        }
        else
        {
            // region Rpp
            type = TestVertexRegion(-cx, -cy, -vx, -vy, ex, ey, ix, iy,
                circle.radius, result.contactTime);
            ix = -ix;
            iy = -iy;
        }
    }

    if (type != 1 || result.contactTime > maxTime)
    {
        result.intersect = false;
        return result;
    }

    result.intersect = true;
    result.contactPoint = box.center + ix*box.axis[0] + iy*box.axis[1];
    return result;
}

template <typename Real>
int FIQuery<Real, OrientedBox2<Real>, Circle2<Real>>::TestVertexRegion(
    Real cx, Real cy, Real vx, Real vy, Real ex, Real ey, Real& ix, Real& iy,
    Real radius, Real& contactTime)
{
    Real dx = cx + ex;
    Real dy = cy + ey;
    Real rsqr = radius*radius;
    Real diff = dx*dx + dy*dy - rsqr;
    if (diff <= (Real)0)
    {
        // Circle is already intersecting the box.
        contactTime = (Real)0;
        return -1;
    }

    Real dot = vx*dx + vy*dy;
    if (dot >= (Real)0)
    {
        // Circle not moving towards box.
        return 0;
    }

    Real dotPerp = vx*dy - vy*dx;
    Real vsqr, inv;

    if (dotPerp >= (Real)0)
    {
        // Potential contact on left edge.
        if (dotPerp <= radius*vy)
        {
            // Lower left corner is first point of contact.
            ix = -ex;
            iy = -ey;
            vsqr = vx*vx + vy*vy;
            inv = ((Real)1) / sqrt(std::abs(dot*dot - vsqr*diff));
            contactTime = diff*inv / ((Real)1 - dot*inv);
            return 1;
        }

        if (vx <= (Real)0)
        {
            // Passed corner, moving away from box.
            return 0;
        }

        vsqr = vx*vx + vy*vy;
        dy = cy - ey;
        dotPerp = vx*dy - vy*dx;
        if (dotPerp >= (Real)0 && dotPerp*dotPerp > rsqr*vsqr)
        {
            // Circle misses box.
            return 0;
        }

        // Circle will intersect box.  Determine first time and place of
        // contact with x = xmin.
        ix = -ex;

        if (dotPerp <= radius*vy)
        {
            // First contact on left edge of box.
            contactTime = -(dx + radius) / vx;
            iy = cy + contactTime*vy;
        }
        else
        {
            // First contact at upper left corner of box.
            dot = vx*dx + vy*dy;
            diff = dx*dx + dy*dy - rsqr;
            inv = ((Real)1) / sqrt(std::abs(dot*dot - vsqr*diff));
            contactTime = diff*inv / ((Real)1 - dot*inv);
            iy = ey;
        }
    }
    else
    {
        // Potential contact on bottom edge.
        if (-dotPerp <= radius*vx)
        {
            // Lower left corner is first point of contact.
            ix = -ex;
            iy = -ey;
            vsqr = vx*vx + vy*vy;
            inv = ((Real)1) / sqrt(std::abs(dot*dot - vsqr*diff));
            contactTime = diff*inv / ((Real)1 - dot*inv);
            return 1;
        }

        if (vy <= (Real)0)
        {
            // Passed corner, moving away from box.
            return 0;
        }

        vsqr = vx*vx + vy*vy;
        dx = cx - ex;
        dotPerp = vx*dy - vy*dx;
        if (-dotPerp >= (Real)0 && dotPerp*dotPerp > rsqr*vsqr)
        {
            // Circle misses box.
            return 0;
        }

        // Circle will intersect box.  Determine first time and place of
        // contact with y = ymin.
        iy = -ey;

        if (-dotPerp <= radius*vx)
        {
            // First contact on bottom edge of box.
            contactTime = -(dy + radius) / vy;
            ix = cx + contactTime*vx;
        }
        else
        {
            // First contact at lower right corner of box.
            dot = vx*dx + vy*dy;
            diff = dx*dx + dy*dy - rsqr;
            inv = ((Real)1) / sqrt(std::abs(dot*dot - vsqr*diff));
            contactTime = diff*inv / ((Real)1 - dot*inv);
            ix = ex;
        }
    }

    return 1;
}

template <typename Real>
int FIQuery<Real, OrientedBox2<Real>, Circle2<Real>>::TestEdgeRegion(Real cx,
    Real cy, Real vx, Real vy, Real ex, Real ey, Real& ix, Real& iy,
    Real radius, Real& contactTime)
{
    Real dx = cx + ex;
    Real xSignedDist = dx + radius;
    if (xSignedDist >= (Real)0)
    {
        // Circle is already intersecting the box.
        contactTime = (Real)0;
        return -1;
    }

    if (vx <= (Real)0)
    {
        // Circle not moving towards box.
        return 0;
    }

    Real rsqr = radius*radius;
    Real vsqr = vx*vx + vy*vy;
    Real dy, dot, dotPerp, diff, inv;

    if (vy >= (Real)0)
    {
        dy = cy - ey;
        dotPerp = vx*dy - vy*dx;
        if (dotPerp >= (Real)0 && dotPerp*dotPerp > rsqr*vsqr)
        {
            // Circle misses box.
            return 0;
        }

        // Circle will intersect box.  Determine first time and place of
        // contact with x = xmin.
        ix = -ex;

        if (dotPerp <= radius*vy)
        {
            // First contact on left edge of box.
            contactTime = -xSignedDist / vx;
            iy = cy + contactTime*vy;
        }
        else
        {
            // First contact at corner of box.
            dot = vx*dx + vy*dy;
            diff = dx*dx + dy*dy - rsqr;
            inv = ((Real)1) / sqrt(std::abs(dot*dot - vsqr*diff));
            contactTime = diff*inv / ((Real)1 - dot*inv);
            iy = ey;
        }
    }
    else
    {
        dy = cy + ey;
        dotPerp = vx*dy - vy*dx;
        if (dotPerp <= (Real)0 && dotPerp*dotPerp > rsqr*vsqr)
        {
            // Circle misses box.
            return 0;
        }

        // Circle will intersect box.  Determine first time and place of
        // contact with x = xmin.
        ix = -ex;

        if (dotPerp >= radius*vy)
        {
            // First contact on left edge of box.
            contactTime = -xSignedDist / vx;
            iy = cy + contactTime*vy;
        }
        else
        {
            // First contact at corner of box.
            dot = vx*dx + vy*dy;
            diff = dx*dx + dy*dy - rsqr;
            inv = ((Real)1) / sqrt(std::abs(dot*dot - vsqr*diff));
            contactTime = diff*inv / ((Real)1 - dot*inv);
            iy = -ey;
        }
    }

    return 1;
}


}
