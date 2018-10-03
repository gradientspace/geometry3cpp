// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.1 (2017/07/25)

#pragma once

#include <Mathematics/GteAlignedBox.h>
#include <Mathematics/GteCone.h>
#include <Mathematics/GteTIQuery.h>
#include <cstdlib>

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
class TIQuery<Real, AlignedBox<3, Real>, Cone<3, Real>>
{
public:
    TIQuery();

    struct Result
    {
        bool intersect;
    };

    Result operator()(AlignedBox<3, Real> const& box, Cone<3, Real> const& cone);

protected:
    // The spherical polygons have vertices stored in counterclockwise order
    // when viewed from outside the sphere.  The 'indices' are lookups into
    // {0..26}, where there are 27 possible spherical polygon configurations
    // based on the location of the cone vertex related to the box.
    struct Polygon
    {
        int numPoints;
        std::array<int, 6> indices;
    };

    // The inputs D and CmV have components relative to the basis of box axes.
    void DoQuery(Vector<3, Real> const& boxExtent, Real cosSqr,
        Vector<3, Real> const& D, Vector<3, Real> const& CmV, Real DdCmV,
        Polygon const& polygon, Result& result);

    std::array<Polygon, 27> mPolygon;
    std::array<int, 4> mMod3;
};


// Template alias for convenience.
template <typename Real>
using TIAlignedBox3Cone3 =
TIQuery<Real, AlignedBox<3, Real>, Cone<3, Real>>;


template <typename Real>
TIQuery<Real, AlignedBox<3, Real>, Cone<3, Real>>::TIQuery()
{
    // The vector z[] stores the box coordinates of the cone vertex, z[i] =
    // Dot(box.axis[i],cone.vertex-box.center).  Each mPolygon[] has a comment
    // with signs:  s[0] s[1] s[2].  If the sign is '-', z[i] < -e[i].  If the
    // sign is '+', z[i] > e[i].  If the sign is '0', |z[i]| <= e[i].
    mPolygon[ 0] = { 6, { 1, 5, 4, 6, 2, 3 } };  // ---
    mPolygon[ 1] = { 6, { 0, 2, 3, 1, 5, 4 } };  // 0--
    mPolygon[ 2] = { 6, { 0, 2, 3, 7, 5, 4 } };  // +--
    mPolygon[ 3] = { 6, { 0, 4, 6, 2, 3, 1 } };  // -0-
    mPolygon[ 4] = { 4, { 0, 2, 3, 1 }       };  // 00-
    mPolygon[ 5] = { 6, { 0, 2, 3, 7, 5, 1 } };  // +0-
    mPolygon[ 6] = { 6, { 0, 4, 6, 7, 3, 1 } };  // -+-
    mPolygon[ 7] = { 6, { 0, 2, 6, 7, 3, 1 } };  // 0+-
    mPolygon[ 8] = { 6, { 0, 2, 6, 7, 5, 1 } };  // ++-
    mPolygon[ 9] = { 6, { 0, 1, 5, 4, 6, 2 } };  // --0
    mPolygon[10] = { 4, { 0, 1, 5, 4 }       };  // 0-0
    mPolygon[11] = { 6, { 0, 1, 3, 7, 5, 4 } };  // +-0
    mPolygon[12] = { 4, { 0, 4, 6, 2 }       };  // -00
    mPolygon[13] = { 0, { 0 }                };  // 000
    mPolygon[14] = { 4, { 1, 3, 7, 5 }       };  // +00
    mPolygon[15] = { 6, { 0, 4, 6, 7, 3, 2 } };  // -+0
    mPolygon[16] = { 4, { 2, 6, 7, 3 }       };  // 0+0
    mPolygon[17] = { 6, { 1, 3, 2, 6, 7, 5 } };  // ++0
    mPolygon[18] = { 6, { 0, 1, 5, 7, 6, 2 } };  // --+
    mPolygon[19] = { 6, { 0, 1, 5, 7, 6, 4 } };  // 0-+
    mPolygon[20] = { 6, { 0, 1, 3, 7, 6, 4 } };  // +-+
    mPolygon[21] = { 6, { 0, 4, 5, 7, 6, 2 } };  // -0+
    mPolygon[22] = { 4, { 4, 5, 7, 6 }       };  // 00+
    mPolygon[23] = { 6, { 1, 3, 7, 6, 4, 5 } };  // +0+
    mPolygon[24] = { 6, { 0, 4, 5, 7, 3, 2 } };  // -++
    mPolygon[25] = { 6, { 2, 6, 4, 5, 7, 3 } };  // 0++
    mPolygon[26] = { 6, { 1, 3, 2, 6, 4, 5 } };  // +++

    // To avoid the modulus operator (%).  For k0 in {0,1,2}, k1 = mod3[k0+1]
    // is (k0+1)%3 and k2 = mod3[k1+1] is (k1+1)%3.
    mMod3 = { 0, 1, 2, 0 };
}

template <typename Real>
typename TIQuery<Real, AlignedBox<3, Real>, Cone<3, Real>>::Result
    TIQuery<Real, AlignedBox<3, Real>, Cone<3, Real>>::operator()(
    AlignedBox<3, Real> const& box, Cone<3, Real> const& cone)
{
    Result result;

    // Get the centered form for the box.
    Vector<3, Real> boxCenter, boxExtent;
    box.GetCenteredForm(boxCenter, boxExtent);

    // Quick-rejection test for boxes below the supporting plane of the cone.
    Vector<3, Real> CmV = boxCenter - cone.ray.origin;
    Real DdCmV = Dot(cone.ray.direction, CmV);  // interval center
    Real radius =  // interval half-length
        boxExtent[0] * std::abs(cone.ray.direction[0]) +
        boxExtent[1] * std::abs(cone.ray.direction[1]) +
        boxExtent[2] * std::abs(cone.ray.direction[2]);
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
    int index[3] = {
        (CmV[0] < -boxExtent[0] ? 2 : (CmV[0] > boxExtent[0] ? 0 : 1)),
        (CmV[1] < -boxExtent[1] ? 2 : (CmV[1] > boxExtent[1] ? 0 : 1)),
        (CmV[2] < -boxExtent[2] ? 2 : (CmV[2] > boxExtent[2] ? 0 : 1))
    };
    int lookup = index[0] + 3 * index[1] + 9 * index[2];
    if (lookup == 13)
    {
        // The cone vertex is in the box.
        result.intersect = true;
        return result;
    }

    Polygon const& polygon = mPolygon[lookup];

    DoQuery(boxExtent, cone.cosAngleSqr, cone.ray.direction, CmV, DdCmV,
        polygon, result);
    return result;
}

template <typename Real>
void TIQuery<Real, AlignedBox<3, Real>, Cone<3, Real>>::DoQuery(
    Vector<3, Real> const& boxExtent, Real cosSqr, Vector<3, Real> const& D,
    Vector<3, Real> const& CmV, Real DdCmV, Polygon const& polygon,
    Result& result)
{
    // Test polygon points.
    Vector<3, Real> X[8], PmV[8];
    Real DdPmV[8], sqrDdPmV[8], sqrLenPmV[8], q;
    int iMax = -1, jMax = -1;
    for (int i = 0; i < polygon.numPoints; ++i)
    {
        int j = polygon.indices[i];
        X[j][0] = (j & 1 ? boxExtent[0] : -boxExtent[0]);
        X[j][1] = (j & 2 ? boxExtent[1] : -boxExtent[1]);
        X[j][2] = (j & 4 ? boxExtent[2] : -boxExtent[2]);
        DdPmV[j] = Dot(D, X[j]) + DdCmV;
        if (DdPmV[j] > (Real)0)
        {
            PmV[j] = X[j] + CmV;
            sqrDdPmV[j] = DdPmV[j] * DdPmV[j];
            sqrLenPmV[j] = Dot(PmV[j], PmV[j]);
            q = sqrDdPmV[j] - cosSqr * sqrLenPmV[j];
            if (q > (Real)0)
            {
                result.intersect = true;
                return;
            }

            // Keep track of the maximum in case we must process box edges.
            // This supports the gradient ascent search.
            if (iMax == -1 ||
                sqrDdPmV[j] * sqrLenPmV[jMax] > sqrDdPmV[jMax] * sqrLenPmV[j])
            {
                iMax = i;
                jMax = j;
            }
        }
    }

    // Theoretically, DoQuery is called when the box has at least one corner
    // above the supporting plane, in which case DdPmV[j] > 0 for at least one
    // j and consequently iMax should not be -1.  But in case of numerical
    // rounding errors, return a no-intersection result if iMax is -1: the
    // box is below the supporting plane within numerical rounding errors.
    if (iMax == -1)
    {
        result.intersect = false;
        return;
    }

    // Start the gradient ascent search at index jMax.
    Real maxSqrLenPmV = sqrLenPmV[jMax];
    Real maxDdPmV = DdPmV[jMax];
    Vector<3, Real>& maxX = X[jMax];
    Vector<3, Real>& maxPmV = PmV[jMax];
    int k0, k1, k2, jDiff;
    Real s, fder, numer, denom, DdMmV, det;
    Vector<3, Real> MmV;

    // Search the counterclockwise edge <corner[jMax],corner[jNext]>.
    int iNext = (iMax < polygon.numPoints - 1 ? iMax + 1 : 0);
    int jNext = polygon.indices[iNext];
    jDiff = jNext - jMax;
    s = (jDiff > 0 ? (Real)1 : (Real)-1);
    k0 = std::abs(jDiff) >> 1;
    fder = s * (D[k0] * maxSqrLenPmV - maxDdPmV * maxPmV[k0]);
    if (fder > (Real)0)
    {
        // The edge has an interior local maximum in F because
        // F(K[j0]) >= F(K[j1]) and the directional derivative of F at K0
        // is positive.  Compute the local maximum point.
        k1 = mMod3[k0 + 1];
        k2 = mMod3[k1 + 1];
        numer = maxPmV[k1] * maxPmV[k1] + maxPmV[k2] * maxPmV[k2];
        denom = D[k1] * maxPmV[k1] + D[k2] * maxPmV[k2];
        MmV[k0] = numer * D[k0];
        MmV[k1] = denom * (maxX[k1] + CmV[k1]);
        MmV[k2] = denom * (maxX[k2] + CmV[k2]);

        // Theoretically, DdMmV > 0, so there is no need to test positivity.
        DdMmV = Dot(D, MmV);
        q = DdMmV * DdMmV - cosSqr * Dot(MmV, MmV);
        if (q > (Real)0)
        {
            result.intersect = true;
            return;
        }

        // Determine on which side of the spherical arc D lives on.  If the
        // polygon side, then the cone ray intersects the polygon and the cone
        // and box intersect.  Otherwise, the D is outside the polygon and the
        // cone and box do not intersect.
        det = s * (D[k1] * maxPmV[k2] - D[k2] * maxPmV[k1]);
        result.intersect = (det <= (Real)0);
        return;
    }

    // Search the clockwise edge <corner[jMax],corner[jPrev]>.
    int iPrev = (iMax > 0 ? iMax - 1 : polygon.numPoints - 1);
    int jPrev = polygon.indices[iPrev];
    jDiff = jMax - jPrev;
    s = (jDiff > 0 ? (Real)1 : (Real)-1);
    k0 = std::abs(jDiff) >> 1;
    fder = -s * (D[k0] * maxSqrLenPmV - maxDdPmV * maxPmV[k0]);
    if (fder > (Real)0)
    {
        // The edge has an interior local maximum in F because
        // F(K[j0]) >= F(K[j1]) and the directional derivative of F at K0
        // is positive.  Compute the local maximum point.
        k1 = mMod3[k0 + 1];
        k2 = mMod3[k1 + 1];
        numer = maxPmV[k1] * maxPmV[k1] + maxPmV[k2] * maxPmV[k2];
        denom = D[k1] * maxPmV[k1] + D[k2] * maxPmV[k2];
        MmV[k0] = numer * D[k0];
        MmV[k1] = denom * (maxX[k1] + CmV[k1]);
        MmV[k2] = denom * (maxX[k2] + CmV[k2]);

        // Theoretically, DdMmV > 0, so there is no need to test positivity.
        DdMmV = Dot(D, MmV);
        q = DdMmV * DdMmV - cosSqr * Dot(MmV, MmV);
        if (q > (Real)0)
        {
            result.intersect = true;
            return;
        }

        // Determine on which side of the spherical arc D lives on.  If the
        // polygon side, then the cone ray intersects the polygon and the cone
        // and box intersect.  Otherwise, the D is outside the polygon and the
        // cone and box do not intersect.
        det = s * (D[k1] * maxPmV[k2] - D[k2] * maxPmV[k1]);
        result.intersect = (det <= (Real)0);
        return;
    }

    result.intersect = false;
}


}
