// Geometric Tools, LLC
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.16.0 (2017/08/24)

#ifndef WM5THPLANE_H
#define WM5THPLANE_H

#include "Wm5MathematicsLIB.h"
#include "Wm5TAPoint.h"

namespace Wm5
{

template <typename Real>
class THPlane
{
public:
    // The plane is represented as Dot(N,X) - c = 0, where N = (n0,n1,n2,0)
    // is a unit-length normal vector, c is the plane constant, and
    // X = (x0,x1,x2,1) is any point on the plane.  The user must ensure
    // that the normal vector is unit length.  The storage as a 4-tuple is
    // (n0,n1,n2,-c).

    // Construction and destruction.
    THPlane()
    {
        // uninitialized members
    }

    THPlane(const THPlane& plane)
        :
        mTuple(plane.mTuple)
    {
    }

    ~THPlane()
    {
    }

    // Specify N and c directly.
    THPlane(Real normal0, Real normal1, Real normal2, Real constant)
    {
        mTuple[0] = normal0;
        mTuple[1] = normal1;
        mTuple[2] = normal2;
        mTuple[3] = -constant;
    }

    THPlane(const TAVector<Real>& normal, Real constant)
    {
        mTuple[0] = normal[0];
        mTuple[1] = normal[1];
        mTuple[2] = normal[2];
        mTuple[3] = -constant;
    }

    // N is specified, c = Dot(N,P) where P = (p0,p1,p2,1) is a point on the
    // plane.
    THPlane(const TAVector<Real>& normal, const TAPoint<Real>& p)
    {
        mTuple[0] = normal[0];
        mTuple[1] = normal[1];
        mTuple[2] = normal[2];
        mTuple[3] = -p.Dot(normal);
    }

    // N = Cross(P1-P0,P2-P0)/Length(Cross(P1-P0,P2-P0)), c = Dot(N,P0) where
    // P0, P1, P2 are points on the plane.
    THPlane(const TAPoint<Real>& p0, const TAPoint<Real>& p1, const TAPoint<Real>& p2)
    {
        TAVector<Real> edge1 = p1 - p0;
        TAVector<Real> edge2 = p2 - p0;
        TAVector<Real> normal = edge1.UnitCross(edge2);
        mTuple[0] = normal[0];
        mTuple[1] = normal[1];
        mTuple[2] = normal[2];
        mTuple[3] = -p0.Dot(normal);
    }

    // Specify the entire (n0,n1,n2,-c) tuple.
    THPlane(const THPoint<Real>& tuple)
        :
        mTuple(tuple)
    {
    }

    // Implicit conversion to THPoint<Real>.
    inline operator THPoint<Real>()
    {
        return mTuple;
    }

    inline operator THPoint<Real>() const
    {
        return mTuple;
    }

    // Coordinate access.
    inline operator const Real* () const
    {
        return (const Real*)mTuple;
    }

    inline operator Real* ()
    {
        return (Real*)mTuple;
    }

    inline const Real& operator[] (int i) const
    {
        return mTuple[i];
    }

    inline Real& operator[] (int i)
    {
        return mTuple[i];
    }

    // Assignment.
    THPlane& operator= (const THPlane& plane)
    {
        mTuple = plane.mTuple;
        return *this;
    }

    // Comparison (for use by STL containers).
    bool operator== (const THPlane& plane) const
    {
        for (int i = 0; i < 4; ++i)
        {
            if (mTuple[i] != plane.mTuple[i])
            {
                return false;
            }
        }
        return true;
    }

    bool operator!= (const THPlane& plane) const
    {
        return !operator==(plane);
    }

    bool operator< (const THPlane& plane) const
    {
        // lexicographical ordering
        for (int i = 0; i < 4; ++i)
        {
            if (mTuple[i] < plane.mTuple[i])
            {
                return true;
            }
            if (mTuple[i] > plane.mTuple[i])
            {
                return false;
            }
        }
        return false;
    }

    bool operator<= (const THPlane& plane) const
    {
        // (x <= y) <=> !(y < x)
        return !(plane.operator<(*this));
    }

    bool operator>  (const THPlane& plane) const
    {
        // (x > y) <=> (y < x)
        return plane.operator<(*this);
    }

    bool operator>= (const THPlane& plane) const
    {
        // (x >= y) <=> !(x < y)
        return !operator<(plane);
    }

    // Access to individual components.
    inline void SetNormal(const TAVector<Real>& normal)
    {
        mTuple[0] = normal[0];
        mTuple[1] = normal[1];
        mTuple[2] = normal[2];
    }

    inline void SetConstant(Real constant)
    {
        mTuple[3] = -constant;
    }

    inline TAVector<Real> GetNormal() const
    {
        return TAVector<Real>(mTuple[0], mTuple[1], mTuple[2]);
    }

    inline Real GetConstant() const
    {
        return -mTuple[3];
    }

    // Compute L = Length(n0,n1,n2) and set the plane to (n0,n1,n2,-c)/L.
    // This is useful when transforming planes by homogeneous matrices, where
    // the unit-length normal is not guaranteed.  The function returns L.
    Real Normalize(const Real epsilon = (Real)0)
    {
        Real length = sqrt(mTuple[0] * mTuple[0] + mTuple[1] * mTuple[1] +
            mTuple[2] * mTuple[2]);

        if (length > epsilon)
        {
            Real invLength = (Real)1 / length;
            mTuple[0] *= invLength;
            mTuple[1] *= invLength;
            mTuple[2] *= invLength;
            mTuple[3] *= invLength;
        }

        return length;
    }

    // Compute d = Dot(N,P)-c where N is the plane normal and c is the plane
    // constant.  This is a signed distance.  The sign of the return value is
    // positive if the point is on the positive side of the plane, negative if
    // the point is on the negative side, and zero if the point is on the
    // plane.
    Real DistanceTo(const TAPoint<Real>& p) const
    {
        return mTuple[0] * p[0] + mTuple[1] * p[1] + mTuple[2] * p[2] + mTuple[3];
    }

    // The "positive side" of the plane is the half space to which the plane
    // normal points.  The "negative side" is the other half space.  The
    // function returns +1 when P is on the positive side, -1 when P is on the
    // the negative side, or 0 when P is on the plane.
    int WhichSide(const TAPoint<Real>& p) const
    {
        Real distance = DistanceTo(p);

        if (distance < (Real)0)
        {
            return -1;
        }
        else if (distance >(Real)0)
        {
            return +1;
        }
        else
        {
            return 0;
        }
    }

private:
    // Storage is (n0,n1,n2,-c).
    THPoint<Real> mTuple;
};

}

#endif
