// Geometric Tools, LLC
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.16.0 (2017/08/24)

#ifndef WM5THPOINT_H
#define WM5THPOINT_H

#include "Wm5MathematicsLIB.h"

namespace Wm5
{

template <typename Real>
class THPoint
{
public:
    // Construction and destruction.  THPoint represents a homogeneous point
    // of the form (x,y,z,w).  Affine points are characterized by w = 1
    // (see class TAPoint) and affine vectors are characterized by w = 0
    // (see class TAVector).  The default constructor does not initialize
    // its members.
    THPoint()
    {
        // uninitialized members
    }

    THPoint(const THPoint& pnt)
    {
        mTuple[0] = pnt.mTuple[0];
        mTuple[1] = pnt.mTuple[1];
        mTuple[2] = pnt.mTuple[2];
        mTuple[3] = pnt.mTuple[3];
    }

    THPoint(Real x, Real y, Real z, Real w)
    {
        mTuple[0] = x;
        mTuple[1] = y;
        mTuple[2] = z;
        mTuple[3] = w;
    }

    ~THPoint()
    {
    }

    // Coordinate access.
    inline operator const Real*() const
    {
        return mTuple;
    }

    inline operator Real*()
    {
        return mTuple;
    }

    inline const Real& operator[] (int i) const
    {
        return mTuple[i];
    }

    inline Real& operator[] (int i)
    {
        return mTuple[i];
    }

    inline Real X() const
    {
        return mTuple[0];
    }

    inline Real& X()
    {
        return mTuple[0];
    }

    inline Real Y() const
    {
        return mTuple[1];
    }

    inline Real& Y()
    {
        return mTuple[1];
    }

    inline Real Z() const
    {
        return mTuple[2];
    }

    inline Real& Z()
    {
        return mTuple[2];
    }

    inline Real W() const
    {
        return mTuple[3];
    }

    inline Real& W()
    {
        return mTuple[3];
    }

    // Assignment.
    THPoint& operator= (const THPoint& pnt)
    {
        mTuple[0] = pnt.mTuple[0];
        mTuple[1] = pnt.mTuple[1];
        mTuple[2] = pnt.mTuple[2];
        mTuple[3] = pnt.mTuple[3];
        return *this;
    }

    // Comparison (for use by STL containers).
    bool operator== (const THPoint& pnt) const
    {
        for (int i = 0; i < 4; ++i)
        {
            if (mTuple[i] != pnt.mTuple[i])
            {
                return false;
            }
        }
        return true;
    }

    bool operator!= (const THPoint& pnt) const
    {
        return !operator==(pnt);
    }

    bool operator< (const THPoint& pnt) const
    {
        // lexicographical ordering
        for (int i = 0; i < 4; ++i)
        {
            if (mTuple[i] < pnt.mTuple[i])
            {
                return true;
            }
            if (mTuple[i] > pnt.mTuple[i])
            {
                return false;
            }
        }
        return false;
    }

    bool operator<= (const THPoint& pnt) const
    {
        // (x <= y) <=> !(y < x)
        return !(pnt.operator<(*this));
    }

    bool operator> (const THPoint& pnt) const
    {
        // (x > y) <=> (y < x)
        return pnt.operator<(*this);
    }

    bool operator>= (const THPoint& pnt) const
    {
        // (x >= y) <=> !(x < y)
        return !operator<(pnt);
    }

protected:
    Real mTuple[4];
};

}

#endif
