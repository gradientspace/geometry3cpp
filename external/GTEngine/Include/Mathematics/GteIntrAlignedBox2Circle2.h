// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.15.1 (2018/08/01)

#pragma once

#include <Mathematics/GteFunctions.h>
#include <Mathematics/GteVector2.h>
#include <Mathematics/GteHypersphere.h>
#include <Mathematics/GteDistPointAlignedBox.h>
#include <Mathematics/GteFIQuery.h>
#include <Mathematics/GteTIQuery.h>

// The find-intersection query is based on the document
// https://www.geometrictools.com/Documentation/IntersectionMovingCircleRectangle.pdf

namespace gte
{
    template <typename Real>
    class TIQuery<Real, AlignedBox2<Real>, Circle2<Real>>
    {
    public:
        // The intersection query considers the box and circle to be solids;
        // that is, the circle object includes the region inside the circular
        // boundary and the box object includes the region inside the
        // rectangular boundary.  If the circle object and rectangle object
        // overlap, the objects intersect.
        struct Result
        {
            bool intersect;
        };

        Result operator()(AlignedBox2<Real> const& box, Circle2<Real> const& circle)
        {
            DCPQuery<Real, Vector2<Real>, AlignedBox2<Real>> pbQuery;
            auto pbResult = pbQuery(circle.center, box);
            Result result;
            result.intersect = (pbResult.sqrDistance <= circle.radius * circle.radius);
            return result;
        }
    };

    template <typename Real>
    class FIQuery<Real, AlignedBox2<Real>, Circle2<Real>>
    {
    public:
        // Currently, only a dynamic query is supported.  A static query will
        // need to compute the intersection set of (solid) box and circle.
        struct Result
        {
            // The cases are
            // 1. Objects initially overlapping.  The contactPoint is invalid.
            //      intersectionType = -1
            //      contactTime = 0
            //      contactPoint = (0,0) (TODO: Set to point in intersection?)
            // 2. Objects initially separated but do not intersect later.  The
            //      contactTime and contactPoint are invalid.
            //      intersectionType = 0
            //      contactTime = 0
            //      contactPoint = (0,0)
            // 3. Objects initially separated but intersect later.
            //      intersectionType = +1
            //      contactTime = first time T > 0
            //      contactPoint = corresponding first contact
            int intersectionType;
            Real contactTime;
            Vector2<Real> contactPoint;

            // TODO: To support arbitrary precision for the contactTime,
            // return q0, q1 and q2 where contactTime = (q0 - sqrt(q1)) / q2.
            // The caller can compute contactTime to desired number of digits
            // of precision.  These are valid in the case when intersectionType
            // is +1 but are set to zero (invalid) in the other cases.  Do the
            // same for the contactPoint.
        };

        Result operator()(AlignedBox2<Real> const& box, Vector2<Real> const& boxVelocity,
            Circle2<Real> const& circle, Vector2<Real> const& circleVelocity)
        {
            // Translate the circle and box so that the box center becomes
            // the origin.  Compute the velocity of the circle relative to
            // the box.
            Real const zero(0), one(1), minusOne(-1), half(0.5);
            Vector2<Real> boxCenter = (box.max + box.min) * half;
            Vector2<Real> extent = (box.max - box.min) * half;
            Vector2<Real> C = circle.center - boxCenter;
            Vector2<Real> V = circleVelocity - boxVelocity;

            // Change signs on components, if necessary, to transform C to the
            // first quadrant.  Adjust the velocity accordingly.
            Real sign[2];
            for (int i = 0; i < 2; ++i)
            {
                if (C[i] >= zero)
                {
                    sign[i] = one;
                }
                else
                {
                    C[i] = -C[i];
                    V[i] = -V[i];
                    sign[i] = minusOne;
                }
            }

            Result result = { 0, zero, { zero, zero } };
            DoQuery(extent, C, circle.radius, V, result);

            if (result.intersectionType != 0)
            {
                // Translate back to the original coordinate system.
                for (int i = 0; i < 2; ++i)
                {
                    if (sign[i] < zero)
                    {
                        result.contactPoint[i] = -result.contactPoint[i];
                    }
                }

                result.contactPoint += boxCenter;
            }
            return result;
        }

    protected:
        void DoQuery(Vector2<Real> const& K, Vector2<Real> const& C,
            Real radius, Vector2<Real> const& V, Result& result)
        {
            Real const zero(0);
            Vector2<Real> delta = C - K;
            if (delta[1] <= radius)
            {
                if (delta[0] <= radius)
                {
                    if (delta[1] <= zero)
                    {
                        if (delta[0] <= zero)
                        {
                            InRegionInteriorOverlap(C, result);
                        }
                        else
                        {
                            InRegionEdgeOverlap(0, 1, C, delta, radius, result);
                        }
                    }
                    else
                    {
                        if (delta[0] <= zero)
                        {
                            InRegionEdgeOverlap(1, 0, C, delta, radius, result);
                        }
                        else
                        {
                            Real sqrRadius = radius * radius;
                            Real sqrDistance = delta[0] * delta[0] + delta[1] * delta[1];
                            if (sqrDistance <= sqrRadius)
                            {
                                InRegionVertexOverlap(K, delta, radius, result);
                            }
                            else
                            {
                                InRegionVertexSeparated(K, delta, V, radius, result);
                            }
                        }

                    }
                }
                else
                {
                    InRegionEdgeUnbounded(0, 1, K, C, radius, delta, V, result);
                }
            }
            else
            {
                if (delta[0] <= radius)
                {
                    InRegionEdgeUnbounded(1, 0, K, C, radius, delta, V, result);
                }
                else
                {
                    InRegionVertexUnbounded(K, C, radius, delta, V, result);
                }
            }
        }

    private:
        void InRegionInteriorOverlap(Vector2<Real> const& C, Result& result)
        {
            Real const zero(0);
            result.intersectionType = -1;
            result.contactTime = zero;
            result.contactPoint = C;
        }

        void InRegionEdgeOverlap(int i0, int i1, Vector2<Real> const& C,
            Vector2<Real> const& delta, Real radius, Result& result)
        {
            Real const zero(0);
            result.intersectionType = (delta[i0] < radius ? -1 : 1);
            result.contactTime = zero;
            result.contactPoint[i0] = C[i0] - delta[i0];
            result.contactPoint[i1] = C[i1];
        }

        void InRegionVertexOverlap(Vector2<Real> const& K0, Vector2<Real> const& delta,
            Real radius, Result& result)
        {
            Real const zero(0);
            Real sqrDistance = delta[0] * delta[0] + delta[1] * delta[1];
            Real sqrRadius = radius * radius;
            result.intersectionType = (sqrDistance < sqrRadius ? -1 : 1);
            result.contactTime = zero;
            result.contactPoint = K0;
        }

        void InRegionVertexSeparated(Vector2<Real> const& K0, Vector2<Real> const& delta0,
            Vector2<Real> const& V, Real radius, Result& result)
        {
            Real const zero(0);
            Real q0 = -Dot(V, delta0);
            if (q0 > zero)
            {
                Real dotVPerpD0 = Dot(V, Perp(delta0));
                Real q2 = Dot(V, V);
                Real q1 = radius * radius * q2 - dotVPerpD0 * dotVPerpD0;
                if (q1 >= zero)
                {
                    IntersectsVertex(0, 1, K0, q0, q1, q2, result);
                }
            }
        }

        void InRegionEdgeUnbounded(int i0, int i1, Vector2<Real> const& K0, Vector2<Real> const& C,
            Real radius, Vector2<Real> const& delta0, Vector2<Real> const& V, Result& result)
        {
            const Real zero(0);
            if (V[i0] < zero)
            {
                Real dotVPerpD0 = V[i0] * delta0[i1] - V[i1] * delta0[i0];
                if (radius * V[1] + dotVPerpD0 >= zero)
                {
                    Vector2<Real> K1, delta1;
                    K1[i0] = K0[i0];
                    K1[i1] = -K0[i1];
                    delta1[i0] = C[i0] - K1[i0];
                    delta1[i1] = C[i1] - K1[i1];
                    Real dotVPerpD1 = V[i0] * delta1[i1] - V[i1] * delta1[i0];
                    if (radius * V[1] + dotVPerpD1 <= zero)
                    {
                        IntersectsEdge(i0, i1, K0, C, radius, V, result);
                    }
                    else
                    {
                        Real q2 = Dot(V, V);
                        Real q1 = radius * radius * q2 - dotVPerpD1 * dotVPerpD1;
                        if (q1 >= zero)
                        {
                            Real q0 = -(V[i0] * delta1[i0] + V[i1] * delta1[i1]);
                            IntersectsVertex(i0, i1, K1, q0, q1, q2, result);
                        }
                    }
                }
                else
                {
                    Real q2 = Dot(V, V);
                    Real q1 = radius * radius * q2 - dotVPerpD0 * dotVPerpD0;
                    if (q1 >= zero)
                    {
                        Real q0 = -(V[i0] * delta0[i0] + V[i1] * delta0[i1]);
                        IntersectsVertex(i0, i1, K0, q0, q1, q2, result);
                    }
                }
            }
        }

        void InRegionVertexUnbounded(Vector2<Real> const& K0, Vector2<Real> const& C, Real radius,
            Vector2<Real> const& delta0, Vector2<Real> const& V, Result& result)
        {
            Real const zero(0);
            if (V[0] < zero && V[1] < zero)
            {
                Real dotVPerpD0 = Dot(V, Perp(delta0));
                if (radius * V[0] - dotVPerpD0 <= zero)
                {
                    if (-radius * V[1] - dotVPerpD0 >= zero)
                    {
                        Real q2 = Dot(V, V);
                        Real q1 = radius * radius * q2 - dotVPerpD0 * dotVPerpD0;
                        Real q0 = -Dot(V, delta0);
                        IntersectsVertex(0, 1, K0, q0, q1, q2, result);
                    }
                    else
                    {
                        Vector2<Real> K1{ K0[0], -K0[1] };
                        Vector2<Real> delta1 = C - K1;
                        Real dotVPerpD1 = Dot(V, Perp(delta1));
                        if (-radius * V[1] - dotVPerpD1 >= zero)
                        {
                            IntersectsEdge(0, 1, K0, C, radius, V, result);
                        }
                        else
                        {
                            Real q2 = Dot(V, V);
                            Real q1 = radius * radius * q2 - dotVPerpD1 * dotVPerpD1;
                            if (q1 >= zero)
                            {
                                Real q0 = -Dot(V, delta1);
                                IntersectsVertex(0, 1, K1, q0, q1, q2, result);
                            }
                        }
                    }
                }
                else
                {
                    Vector2<Real> K2{ -K0[0], K0[1] };
                    Vector2<Real> delta2 = C - K2;
                    Real dotVPerpD2 = Dot(V, Perp(delta2));
                    if (radius * V[0] - dotVPerpD2 <= zero)
                    {
                        IntersectsEdge(1, 0, K0, C, radius, V, result);
                    }
                    else
                    {
                        Real q2 = Dot(V, V);
                        Real q1 = radius * radius * q2 - dotVPerpD2 * dotVPerpD2;
                        if (q1 >= zero)
                        {
                            Real q0 = -Dot(V, delta2);
                            IntersectsVertex(1, 0, K2, q0, q1, q2, result);
                        }
                    }
                }
            }
        }

        void IntersectsVertex(int i0, int i1, Vector2<Real> const& K,
            Real q0, Real q1, Real q2, Result& result)
        {
            result.intersectionType = +1;
            result.contactTime = (q0 - Function<Real>::Sqrt(q1)) / q2;
            result.contactPoint[i0] = K[i0];
            result.contactPoint[i1] = K[i1];
        }

        void IntersectsEdge(int i0, int i1, Vector2<Real> const& K0, Vector2<Real> const& C,
            Real radius, Vector2<Real> const& V, Result& result)
        {
            result.intersectionType = +1;
            result.contactTime = (K0[i0] + radius - C[i0]) / V[i0];
            result.contactPoint[i0] = K0[i0];
            result.contactPoint[i1] = C[i1] + result.contactTime * V[i1];
        }
    };
}
