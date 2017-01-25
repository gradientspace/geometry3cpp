// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <Mathematics/GteVector3.h>
#include <Mathematics/GteHypersphere.h>
#include <Mathematics/GteOrientedBox.h>
#include <Mathematics/GteFIQuery.h>
#include <Mathematics/GteTIQuery.h>

namespace gte
{

template <typename Real>
class TIQuery<Real, OrientedBox3<Real>, Sphere3<Real>>
{
public:
    // The intersection query considers the box and sphere to be solids.
    // For example, if the sphere is strictly inside the box (does not touch
    // the box faces), the objects intersect.
    struct Result
    {
        bool intersect;
    };

    Result operator()(OrientedBox3<Real> const& box,
        Sphere3<Real> const& sphere);
};

template <typename Real>
class FIQuery<Real, OrientedBox3<Real>, Sphere3<Real>>
{
public:
    // Currently, only a dynamic query is supported.  The static query must
    // compute the intersection set of (solid) box and sphere.
    struct Result
    {
        bool intersect;
        Real contactTime;
        Vector3<Real> contactPoint;
    };

    Result operator()(Real maxTime, OrientedBox3<Real> const& box,
        Vector3<Real> const& boxVelocity, Sphere3<Real> const& sphere,
        Vector3<Real> const& sphereVelocity);

private:
    // Support for dynamic query.
    Real GetVertexIntersection(Real dx, Real dy, Real dz, Real vx,
        Real vy, Real vz, Real rsqr);

    Real GetEdgeIntersection(Real dx, Real dz, Real vx, Real vz,
        Real vsqr, Real rsqr);

    int FindFaceRegionIntersection(Real ex, Real ey, Real ez, Real cx,
        Real cy, Real cz, Real vx, Real vy, Real vz, Real& ix, Real& iy,
        Real& iz, bool aboveFace, Real radius, Real& contactTime);

    int FindJustEdgeIntersection(Real cy, Real ex, Real ey, Real ez, Real dx,
        Real dz, Real vx, Real vy, Real vz, Real& ix, Real& iy, Real& iz,
        Real radius, Real& contactTime);

    int FindEdgeRegionIntersection(Real ex, Real ey, Real ez, Real cx,
        Real cy, Real cz, Real vx, Real vy, Real vz, Real& ix, Real& iy,
        Real& iz, bool aboveEdge, Real radius, Real& contactTime);

    int FindVertexRegionIntersection(Real ex, Real ey, Real ez, Real cx,
        Real cy, Real cz, Real vx, Real vy, Real vz, Real& ix, Real& iy,
        Real& iz, Real radius, Real& contactTime);
};


template <typename Real>
typename TIQuery<Real, OrientedBox3<Real>, Sphere3<Real>>::Result
TIQuery<Real, OrientedBox3<Real>, Sphere3<Real>>::operator()(
    OrientedBox3<Real> const& box, Sphere3<Real> const& sphere)
{
    Result result;

    // Test for intersection in the coordinate system of the box by
    // transforming the sphere into that coordinate system.
    Vector3<Real> cdiff = sphere.center - box.center;
    Real ax = std::abs(Dot(cdiff, box.axis[0]));
    Real ay = std::abs(Dot(cdiff, box.axis[1]));
    Real az = std::abs(Dot(cdiff, box.axis[2]));
    Real dx = ax - box.extent[0];
    Real dy = ay - box.extent[1];
    Real dz = az - box.extent[2];

    if (ax <= box.extent[0])
    {
        if (ay <= box.extent[1])
        {
            if (az <= box.extent[2])
            {
                // Sphere center inside box.
                result.intersect = true;
                return result;
            }
            else
            {
                // Potential sphere-face intersection with face z.
                result.intersect = (dz <= sphere.radius);
                return result;
            }
        }
        else
        {
            if (az <= box.extent[2])
            {
                // Potential sphere-face intersection with face y.
                result.intersect = (dy <= sphere.radius);
                return result;
            }
            else
            {
                // Potential sphere-edge intersection with edge formed
                // by faces y and z.
                Real rsqr = sphere.radius*sphere.radius;
                result.intersect = (dy*dy + dz*dz <= rsqr);
                return result;
            }
        }
    }
    else
    {
        if (ay <= box.extent[1])
        {
            if (az <= box.extent[2])
            {
                // Potential sphere-face intersection with face x.
                result.intersect = (dx <= sphere.radius);
                return result;
            }
            else
            {
                // Potential sphere-edge intersection with edge formed
                // by faces x and z.
                Real rsqr = sphere.radius*sphere.radius;
                result.intersect = (dx*dx + dz*dz <= rsqr);
                return result;
            }
        }
        else
        {
            if (az <= box.extent[2])
            {
                // Potential sphere-edge intersection with edge formed
                // by faces x and y.
                Real rsqr = sphere.radius*sphere.radius;
                result.intersect = (dx*dx + dy*dy <= rsqr);
                return result;
            }
            else
            {
                // Potential sphere-vertex intersection at corner formed
                // by faces x,y,z.
                Real rsqr = sphere.radius*sphere.radius;
                result.intersect = (dx*dx + dy*dy + dz*dz <= rsqr);
                return result;
            }
        }
    }
}



template <typename Real>
typename FIQuery<Real, OrientedBox3<Real>, Sphere3<Real>>::Result
FIQuery<Real, OrientedBox3<Real>, Sphere3<Real>>::operator()(Real maxTime,
    OrientedBox3<Real> const& box, Vector3<Real> const& boxVelocity,
    Sphere3<Real> const& sphere, Vector3<Real> const& sphereVelocity)
{
    Result result;

    // Find intersections relative to the coordinate system of the box.
    // The sphere is transformed to the box coordinates and the velocity of
    // the sphere is relative to the box.
    Vector3<Real> cdiff = sphere.center - box.center;
    Vector3<Real> relativeVelocity = sphereVelocity - boxVelocity;
    Real ax = Dot(cdiff, box.axis[0]);
    Real ay = Dot(cdiff, box.axis[1]);
    Real az = Dot(cdiff, box.axis[2]);
    Real vx = Dot(relativeVelocity, box.axis[0]);
    Real vy = Dot(relativeVelocity, box.axis[1]);
    Real vz = Dot(relativeVelocity, box.axis[2]);

    // Flip coordinate frame into the first octant.
    int signX = 1;
    if (ax < (Real)0)
    {
        ax = -ax;
        vx = -vx;
        signX = -1;
    }

    int signY = 1;
    if (ay < (Real)0)
    {
        ay = -ay;
        vy = -vy;
        signY = -1;
    }

    int signZ = 1;
    if (az < (Real)0)
    {
        az = -az;
        vz = -vz;
        signZ = -1;
    }

    // Intersection coordinates.
    Real ix, iy, iz;
    int retVal;

    if (ax <= box.extent[0])
    {
        if (ay <= box.extent[1])
        {
            if (az <= box.extent[2])
            {
                // The sphere center is inside box.  Use it as the contact
                // point.
                result.intersect = true;
                result.contactTime = (Real)0;
                result.contactPoint = sphere.center;
                return result;
            }
            else
            {
                // Sphere above face on axis Z.
                retVal = FindFaceRegionIntersection(box.extent[0],
                    box.extent[1], box.extent[2], ax, ay, az, vx, vy,
                    vz, ix, iy, iz, true, sphere.radius, result.contactTime);
            }
        }
        else
        {
            if (az <= box.extent[2])
            {
                // Sphere above face on axis Y.
                retVal = FindFaceRegionIntersection(box.extent[0],
                    box.extent[2], box.extent[1], ax, az, ay, vx, vz,
                    vy, ix, iz, iy, true, sphere.radius, result.contactTime);
            }
            else
            {
                // Sphere is above the edge formed by faces y and z.
                retVal = FindEdgeRegionIntersection(box.extent[1],
                    box.extent[0], box.extent[2], ay, ax, az, vy, vx,
                    vz, iy, ix, iz, true, sphere.radius, result.contactTime);
            }
        }
    }
    else
    {
        if (ay <= box.extent[1])
        {
            if (az <= box.extent[2])
            {
                // Sphere above face on axis X.
                retVal = FindFaceRegionIntersection(box.extent[1],
                    box.extent[2], box.extent[0], ay, az, ax, vy, vz,
                    vx, iy, iz, ix, true, sphere.radius, result.contactTime);
            }
            else
            {
                // Sphere is above the edge formed by faces x and z.
                retVal = FindEdgeRegionIntersection(box.extent[0],
                    box.extent[1], box.extent[2], ax, ay, az, vx, vy,
                    vz, ix, iy, iz, true, sphere.radius, result.contactTime);
            }
        }
        else
        {
            if (az <= box.extent[2])
            {
                // Sphere is above the edge formed by faces x and y.
                retVal = FindEdgeRegionIntersection(box.extent[0],
                    box.extent[2], box.extent[1], ax, az, ay, vx, vz,
                    vy, ix, iz, iy, true, sphere.radius, result.contactTime);
            }
            else
            {
                // sphere is above the corner formed by faces x,y,z
                retVal = FindVertexRegionIntersection(box.extent[0],
                    box.extent[1], box.extent[2], ax, ay, az, vx, vy,
                    vz, ix, iy, iz, sphere.radius, result.contactTime);
            }
        }
    }

    if (retVal == 0 || result.contactTime > maxTime)
    {
        result.intersect = false;
        return result;
    }

    // Calculate actual intersection (move point back into world coordinates).
    result.intersect = true;
    result.contactPoint = box.center + (signX*ix)*box.axis[0] +
        (signY*iy)*box.axis[1] + (signZ*iz)*box.axis[2];
    return result;
}

template <typename Real>
Real FIQuery<Real, OrientedBox3<Real>, Sphere3<Real>>::GetVertexIntersection(
    Real dx, Real dy, Real dz, Real vx, Real vy, Real vz, Real rsqr)
{
    // Find the time of a 3D line-sphere intersection between a line P = Dt,
    // where P = (dx, dy, dz) and D = (vx, vy, vz) and a sphere of radius^2
    // rsqr.  This is valid only if there is, in fact, an intersection.

    Real vsqr = vx*vx + vy*vy + vz*vz;
    Real dot = dx*vx + dy*vy + dz*vz;
    Real diff = dx*dx + dy*dy + dz*dz - rsqr;
    Real inv = ((Real)1) / sqrt(std::abs(dot*dot - vsqr*diff));
    return diff*inv / ((Real)1 - dot*inv);
}

template <typename Real>
Real FIQuery<Real, OrientedBox3<Real>, Sphere3<Real>>::GetEdgeIntersection(
    Real dx, Real dz, Real vx, Real vz, Real vsqr, Real rsqr)
{
    // Find the time of a 2D line-circle intersection between a line
    // P = Dt where P = (dx,dz) and D = (vx, vz) and a circle of radius^2
    // rsqr.  This is valid only if there is, in fact, an intersection.

    Real dot = vx*dx + vz*dz;
    Real diff = dx*dx + dz*dz - rsqr;
    Real inv = ((Real)1) / sqrt(std::abs(dot*dot - vsqr*diff));
    return diff*inv / ((Real)1 - dot*inv);
}

template <typename Real>
int FIQuery<Real, OrientedBox3<Real>, Sphere3<Real>>::
FindFaceRegionIntersection(Real ex, Real ey, Real ez, Real cx, Real cy,
    Real cz, Real vx, Real vy, Real vz, Real& ix, Real& iy, Real& iz,
    bool aboveFace, Real radius, Real& contactTime)
{
    // Returns when and whether a sphere in the region above face +Z
    // intersects face +Z or any of its vertices or edges.  The input
    // aboveFace is true when the x and y coordinates are within the x and y
    // extents.  The function will still work if they are not, but it needs
    // to be false then, to avoid some checks that assume that x and y are
    // within the extents.  This function checks face z, and the vertex and
    // two edges that the velocity is headed towards on the face.

    // Check for already intersecting if above face.
    if (cz <= ez + radius && aboveFace)
    {
        contactTime = (Real)0;
        return -1;
    }

    // Check for easy out (moving away on Z axis).
    if (vz >= (Real)0)
    {
        return 0;
    }

    Real rsqr = radius*radius;

    Real vsqrX = vz*vz + vx*vx;
    Real vsqrY = vz*vz + vy*vy;
    Real dx, dy, dz = cz - ez;
    Real crossX, crossY;
    int signX, signY;

    // This determines which way the box is heading and finds the values of
    // CrossX and CrossY which are positive if the sphere center will not
    // pass through the box.  Then it is only necessary to check two edges,
    // the face and the vertex for intersection.

    if (vx >= (Real)0)
    {
        signX = 1;
        dx = cx - ex;
        crossX = vx*dz - vz*dx;
    }
    else
    {
        signX = -1;
        dx = cx + ex;
        crossX = vz*dx - vx*dz;
    }

    if (vy >= (Real)0)
    {
        signY = 1;
        dy = cy - ey;
        crossY = vy*dz - vz*dy;
    }
    else
    {
        signY = -1;
        dy = cy + ey;
        crossY = vz*dy - vy*dz;
    }

    // Does the circle intersect along the x edge?
    if (crossX > radius*vx*signX)
    {
        if (crossX*crossX > rsqr*vsqrX)
        {
            // Sphere overshoots box on the x-axis (either side).
            return 0;
        }

        // Does the circle hit the y edge?
        if (crossY > radius*vy*signY)
        {
            // Potential vertex intersection.
            if (crossY*crossY > rsqr*vsqrY)
            {
                // Sphere overshoots box on the y-axis (either side).
                return 0;
            }

            Vector3<Real> relVelocity{ vx, vy, vz };
            Vector3<Real> D{ dx, dy, dz };
            Vector3<Real> cross = Cross(D, relVelocity);
            Real crossSqrLen = Dot(cross, cross);
            Real relvelSqrLen = Dot(relVelocity, relVelocity);
            if (crossSqrLen > rsqr*relvelSqrLen)
            {
                // Sphere overshoots the box on the corner.
                return 0;
            }

            contactTime = GetVertexIntersection(dx, dy, dz, vx, vy, vz, rsqr);
            ix = ex*signX;
            iy = ey*signY;
        }
        else
        {
            // x-edge intersection
            contactTime = GetEdgeIntersection(dx, dz, vx, vz, vsqrX, rsqr);
            ix = ex*signX;
            iy = cy + vy*contactTime;
        }
    }
    else
    {
        // Does the circle hit the y edge?
        if (crossY > radius*vy*signY)
        {
            // Potential y-edge intersection.
            if (crossY*crossY > rsqr*vsqrY)
            {
                // Sphere overshoots box on the y-axis (either side).
                return 0;
            }

            contactTime = GetEdgeIntersection(dy, dz, vy, vz, vsqrY, rsqr);
            ix = cx + vx*contactTime;
            iy = ey*signY;
        }
        else
        {
            // Face intersection (easy).
            contactTime = (-dz + radius) / vz;
            ix = contactTime*vx + cx;
            iy = contactTime*vy + cy;
        }
    }

    // z coordinate of any intersection must be the face of z.
    iz = ez;
    return 1;
}

template <typename Real>
int FIQuery<Real, OrientedBox3<Real>, Sphere3<Real>>::
FindJustEdgeIntersection(Real cy, Real ex, Real ey, Real ez, Real dx, Real dz,
    Real vx, Real vy, Real vz, Real& ix, Real& iy, Real& iz, Real radius,
    Real& contactTime)
{
    // Finds the intersection of a point dx and dz away from an edge with
    // direction y.  The sphere is at a point cy, and the edge is at the
    // point ex.  Checks the edge and the vertex the velocity is heading
    // towards.

    Real rsqr = radius*radius;
    Real dy, crossZ, crossX;  // possible edge/vertex intersection
    int signY;

    // Depending on the sign of Vy, pick the vertex that the velocity is
    // heading towards on the edge, as well as creating crossX and crossZ
    // such that their sign will always be positive if the sphere center goes
    // over that edge.

    if (vy >= (Real)0)
    {
        signY = 1;
        dy = cy - ey;
        crossZ = dx*vy - dy*vx;
        crossX = dz*vy - dy*vz;
    }
    else
    {
        signY = -1;
        dy = cy + ey;
        crossZ = dy*vx - dx*vy;
        crossX = dy*vz - dz*vy;
    }

    // Check where on edge this intersection will occur.
    if (crossZ >= (Real)0 && crossX >= (Real)0
        && crossX*crossX + crossZ*crossZ > rsqr*vy*vy)
    {
        // Sphere potentially intersects with vertex.
        Vector3<Real> relVelocity{ vx, vy, vz };
        Vector3<Real> D{ dx, dy, dz };
        Vector3<Real> cross = Cross(D, relVelocity);
        Real crossSqrLen = Dot(cross, cross);
        Real relvelSqrLen = Dot(relVelocity, relVelocity);
        if (crossSqrLen > rsqr*relvelSqrLen)
        {
            // Sphere overshoots the box on the vertex.
            return 0;
        }

        // Sphere actually does intersect the vertex.
        contactTime = GetVertexIntersection(dx, dy, dz, vx, vy, vz, rsqr);
        ix = ex;
        iy = signY*ey;
        iz = ez;
    }
    else
    {
        // Sphere intersects with edge.
        Real vsqrX = vz*vz + vx*vx;
        contactTime = GetEdgeIntersection(dx, dz, vx, vz, vsqrX, rsqr);
        ix = ex;
        iy = cy + contactTime*vy;
        iz = ez;
    }
    return 1;
}

template <typename Real>
int FIQuery<Real, OrientedBox3<Real>, Sphere3<Real>>::
FindEdgeRegionIntersection(Real ex, Real ey, Real ez, Real cx, Real cy,
    Real cz, Real vx, Real vy, Real vz, Real& ix, Real& iy, Real& iz,
    bool aboveEdge, Real radius, Real& contactTime)
{
    // Assumes the sphere center is in the region above the x and z planes.
    // The input aboveEdge is true when the y coordinate is within the y
    // extents.  The function will still work if it is not, but it needs to be
    // false then, to avoid some checks that assume that y is within the
    // extent.  This function checks the edge that the region is above, as
    // well as does a "face region" check on the face it is heading towards.

    Real dx = cx - ex;
    Real dz = cz - ez;
    Real rsqr = radius*radius;

    if (aboveEdge)
    {
        Real diff = dx*dx + dz*dz - rsqr;
        if (diff <= (Real)0)
        {
            // Circle is already intersecting the box.
            contactTime = (Real)0;
            return -1;
        }
    }

    Real dot = vx*dx + vz*dz;
    if (dot >= (Real)0)
    {
        // Circle not moving towards box.
        return 0;
    }

    // The value dotPerp splits the region of interest along the edge in the
    // middle of that region.
    Real dotPerp = vz*dx - vx*dz;
    if (dotPerp >= (Real)0)
    {
        // Sphere moving towards +z face.
        if (vx >= (Real)0)
        {
            // Passed corner, moving away from box.
            return 0;
        }

        // Intersection with x-z edge.  If there is trouble with objects that
        // "scrape" the surface (velocity perpendicular to face normal, and
        // point of contact with a radius direction parallel to the face
        // normal), this check may need to be a little more inclusive (small
        // tolerance due to floating point errors) as the edge check needs
        // to get "scraping" objects (as they hit the edge with the point)
        // and the face region check will not catch it because the object is
        // not moving towards the face.
        if (dotPerp <= -radius*vx)
        {
            return FindJustEdgeIntersection(cy, ez, ey, ex, dz, dx, vz, vy,
                vx, iz, iy, ix, radius, contactTime);
        }

        // Now, check the face of z for intersections.
        return FindFaceRegionIntersection(ex, ey, ez, cx, cy, cz, vx, vy,
            vz, ix, iy, iz, false, radius, contactTime);
    }
    else
    {
        // Sphere moving towards +x face.
        if (vz >= (Real)0)
        {
            // Passed corner, moving away from box.
            return 0;
        }

        // Check intersection with x-z edge.  See the note above about
        // "scraping" objects.
        if (dotPerp >= radius*vz)
        {
            // Possible edge/vertex intersection.
            return FindJustEdgeIntersection(cy, ex, ey, ez, dx, dz, vx, vy,
                vz, ix, iy, iz, radius, contactTime);
        }

        // Now, check the face of x for intersections.
        return FindFaceRegionIntersection(ez, ey, ex, cz, cy, cx, vz, vy,
            vx, iz, iy, ix, false, radius, contactTime);
    }
}

template <typename Real>
int FIQuery<Real, OrientedBox3<Real>, Sphere3<Real>>::
FindVertexRegionIntersection(Real ex, Real ey, Real ez, Real cx, Real cy,
    Real cz, Real vx, Real vy, Real vz, Real& ix, Real& iy, Real& iz,
    Real radius, Real& contactTime)
{
    // Assumes the sphere is above the vertex +ex, +ey, +ez.

    Real dx = cx - ex;
    Real dy = cy - ey;
    Real dz = cz - ez;
    Real rsqr = radius*radius;
    Real diff = dx*dx + dy*dy + dz*dz - rsqr;
    if (diff <= (Real)0)
    {
        // Sphere is already intersecting the box.
        contactTime = (Real)0;
        return -1;
    }

    if (vx*dx + vy*dy + vz*dz >= (Real)0)
    {
        // Sphere not moving towards box.
        return 0;
    }

    // The box can be divided up into 3 regions, which simplifies checking to
    // see what the sphere hits.  The regions are divided by which edge
    // (formed by the vertex and some axis) is closest to the velocity
    // vector.
    //
    // To check if it hits the vertex, look at the edge (E) it is going
    // towards.  Create a plane formed by the other two edges (with E as the
    // plane normal) with the vertex at the origin.  The intercept along an
    // axis, in that plane, of the line formed by the sphere center and the
    // velocity as its direction, will be crossAxis/vEdge.  Thus, the
    // distance from the origin to the point in the plane that that line from
    // the sphere in the velocity direction crosses will be the squared sum
    // of those two intercepts.  If that sum is less than the radius squared,
    // then the sphere will hit the vertex otherwise, it will continue on
    // past the vertex.  If it misses, since it is known which edge the box
    // is near, the find edge region test can be used.
    //
    // Also, due to the constrained conditions, only those six cases (two for
    // each region, since crossEdge can be + or -) of signs of cross values
    // can occur.
    //
    // The third case will also pick up cases where D = V, causing a zero
    // cross.

    Real crossX = vy*dz - vz*dy;
    Real crossY = vx*dz - vz*dx;
    Real crossZ = vy*dx - vx*dy;
    Real crX2 = crossX*crossX;
    Real crY2 = crossY*crossY;
    Real crZ2 = crossZ*crossZ;
    Real vx2 = vx*vx;
    Real vy2 = vy*vy;
    Real vz2 = vz*vz;

    // Intersection with the vertex?
    if ((crossY <  (Real)0 && crossZ >= (Real)0 && crY2 + crZ2 <= rsqr*vx2)
        || (crossZ <  (Real)0 && crossX <  (Real)0 && crX2 + crZ2 <= rsqr*vy2)
        || (crossY >= (Real)0 && crossX >= (Real)0 && crX2 + crY2 <= rsqr*vz2))
    {
        // Standard line-sphere intersection.
        contactTime = GetVertexIntersection(dx, dy, dz, vx, vy, vz, rsqr);
        ix = contactTime*vx + cx;
        iy = contactTime*vy + cy;
        iz = contactTime*vz + cz;
        return 1;
    }
    else if (crossY < (Real)0 && crossZ >= (Real)0)
    {
        // x edge region, check y,z planes.
        return FindEdgeRegionIntersection(ey, ex, ez, cy, cx, cz, vy, vx,
            vz, iy, ix, iz, false, radius, contactTime);
    }
    else if (crossZ < (Real)0 && crossX < (Real)0)
    {
        // y edge region, check x,z planes.
        return FindEdgeRegionIntersection(ex, ey, ez, cx, cy, cz, vx, vy,
            vz, ix, iy, iz, false, radius, contactTime);
    }
    else // crossY >= 0 && crossX >= 0
    {
        // z edge region, check x,y planes.
        return FindEdgeRegionIntersection(ex, ez, ey, cx, cz, cy, vx, vz,
            vy, ix, iz, iy, false, radius, contactTime);
    }
}


}
