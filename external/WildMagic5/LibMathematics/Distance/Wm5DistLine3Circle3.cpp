// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.2 (2013/07/28)

#include "Wm5MathematicsPCH.h"
#include "Wm5DistLine3Circle3.h"
#include "Wm5PolynomialRoots.h"
#include "Wm5DistPoint3Circle3.h"

namespace Wm5
{
//----------------------------------------------------------------------------
template <typename Real>
DistLine3Circle3<Real>::DistLine3Circle3 (const Line3<Real>& rkLine,
    const Circle3<Real>& rkCircle)
    :
    mLine(&rkLine),
    mCircle(&rkCircle),
    mNumClosestLine(0),
    mNumClosestCircle(0)
{
}
//----------------------------------------------------------------------------
template <typename Real>
const Line3<Real>& DistLine3Circle3<Real>::GetLine () const
{
    return *mLine;
}
//----------------------------------------------------------------------------
template <typename Real>
const Circle3<Real>& DistLine3Circle3<Real>::GetCircle () const
{
    return *mCircle;
}
//----------------------------------------------------------------------------
template <typename Real>
Real DistLine3Circle3<Real>::Get ()
{
    return Math<Real>::Sqrt(GetSquared());
}
//----------------------------------------------------------------------------
template <typename Real>
Real DistLine3Circle3<Real>::GetSquared ()
{
    Vector3<Real> D = mLine->Origin - mCircle->Center;
    Real s, sqrDistance, temp;

    Vector3<Real> MxN = mLine->Direction.Cross(mCircle->Normal);
    Real m0sqr = MxN.Dot(MxN);
    if (m0sqr > (Real)0)
    {
        Real m0 = sqrt(m0sqr);
        Real rm0 = mCircle->Radius*m0;
        Vector3<Real> DxN = D.Cross(mCircle->Normal);
        Real lambda = -MxN.Dot(DxN)/m0sqr;
        D += lambda*mLine->Direction;
        DxN += lambda*MxN;
        Real m2b2 = mLine->Direction.Dot(D);
        Real b1sqr = DxN.Dot(DxN);
        if (b1sqr > (Real)0)
        {
            Real b1 = sqrt(b1sqr);
            Real rm0sqr = mCircle->Radius*m0sqr;
            if (rm0sqr > b1)
            {
                const Real twoThirds = (Real)2/(Real)3;
                Real sHat = sqrt(pow(rm0sqr*b1sqr, twoThirds) - b1sqr)/m0;
                Real gHat = rm0sqr*sHat/sqrt(m0sqr*sHat*sHat + b1sqr);
                Real cutoff = gHat - sHat;
                if (m2b2 <= -cutoff)
                {
                    s = BisectF(m2b2, rm0sqr, m0sqr, b1sqr, -m2b2,
                        -m2b2 + rm0);
                    mNumClosestLine = 1;
                    mClosestLine[0] =
                        mLine->Origin + (s + lambda)*mLine->Direction;
                    sqrDistance = SqrDistancePointCircle(0);
                    if (m2b2 == -cutoff)
                    {
                        mClosestLine[1] =
                            mLine->Origin + (-sHat + lambda)*mLine->Direction;
                        temp = SqrDistancePointCircle(1);
                        if (temp < sqrDistance)
                        {
                            sqrDistance = temp;
                            std::swap(mClosestLine[0], mClosestLine[1]);
                            std::swap(mClosestCircle[0], mClosestCircle[1]);
                        }
                        else if (temp == sqrDistance)
                        {
                            mNumClosestLine = 2;
                            mNumClosestCircle = 2;
                        }
                    }
                }
                else if (m2b2 >= cutoff)
                {
                    s = BisectF(m2b2, rm0sqr, m0sqr, b1sqr, -m2b2 - rm0,
                        -m2b2);
                    mNumClosestLine = 1;
                    mClosestLine[0] =
                        mLine->Origin + (s + lambda)*mLine->Direction;
                    sqrDistance = SqrDistancePointCircle(0);
                    if (m2b2 == cutoff)
                    {
                        mClosestLine[1] =
                            mLine->Origin + (+sHat + lambda)*mLine->Direction;
                        temp = SqrDistancePointCircle(1);
                        if (temp < sqrDistance)
                        {
                            sqrDistance = temp;
                            std::swap(mClosestLine[0], mClosestLine[1]);
                            std::swap(mClosestCircle[0], mClosestCircle[1]);
                        }
                        else if (temp == sqrDistance)
                        {
                            mNumClosestLine = 2;
                            mNumClosestCircle = 2;
                        }
                    }
                }
                else
                {
                    if (m2b2 <= (Real)0)
                    {
                        s = BisectF(m2b2, rm0sqr, m0sqr, b1sqr, -m2b2,
                            -m2b2 + rm0);
                        mNumClosestLine = 1;
                        mClosestLine[0] =
                            mLine->Origin + (s + lambda)*mLine->Direction;
                        sqrDistance = SqrDistancePointCircle(0);
                        s = BisectF(m2b2, rm0sqr, m0sqr, b1sqr, -m2b2 - rm0,
                            -sHat);
                        mClosestLine[1] =
                            mLine->Origin + (s + lambda)*mLine->Direction;
                        temp = SqrDistancePointCircle(1);
                        if (temp < sqrDistance)
                        {
                            sqrDistance = temp;
                            std::swap(mClosestLine[0], mClosestLine[1]);
                            std::swap(mClosestCircle[0], mClosestCircle[1]);
                        }
                        else if (temp == sqrDistance)
                        {
                            mNumClosestLine = 2;
                            mNumClosestCircle = 2;
                        }
                    }
                    else
                    {
                        s = BisectF(m2b2, rm0sqr, m0sqr, b1sqr, -m2b2 - rm0,
                            -m2b2);
                        mNumClosestLine = 1;
                        mClosestLine[0] =
                            mLine->Origin + (s + lambda)*mLine->Direction;
                        sqrDistance = SqrDistancePointCircle(0);
                        s = BisectF(m2b2, rm0sqr, m0sqr, b1sqr, sHat,
                            -m2b2 + rm0);
                        mClosestLine[1] =
                            mLine->Origin + (s + lambda)*mLine->Direction;
                        temp = SqrDistancePointCircle(1);
                        if (temp < sqrDistance)
                        {
                            sqrDistance = temp;
                            std::swap(mClosestLine[0], mClosestLine[1]);
                            std::swap(mClosestCircle[0], mClosestCircle[1]);
                        }
                        else if (temp == sqrDistance)
                        {
                            mNumClosestLine = 2;
                            mNumClosestCircle = 2;
                        }
                    }
                }
            }
            else
            {
                if (m2b2 < (Real)0)
                {
                    s = BisectF(m2b2, rm0sqr, m0sqr, b1sqr, -m2b2,
                        -m2b2 + rm0);
                    mNumClosestLine = 1;
                    mClosestLine[0] =
                        mLine->Origin + (s + lambda)*mLine->Direction;
                    sqrDistance = SqrDistancePointCircle(0);
                }
                else if (m2b2 > (Real)0)
                {
                    s = BisectF(m2b2, rm0sqr, m0sqr, b1sqr, -m2b2 - rm0,
                        -m2b2);
                    mNumClosestLine = 1;
                    mClosestLine[0] =
                        mLine->Origin + (s + lambda)*mLine->Direction;
                    sqrDistance = SqrDistancePointCircle(0);
                }
                else
                {
                    mNumClosestLine = 1;
                    mClosestLine[0] = mLine->Origin + lambda*mLine->Direction;
                    sqrDistance = SqrDistancePointCircle(0);
                }
            }
        }
        else
        {
            if (m2b2 < (Real)0)
            {
                s = -m2b2 + rm0;
                mNumClosestLine = 1;
                mClosestLine[0] =
                    mLine->Origin + (s + lambda)*mLine->Direction;
                sqrDistance = SqrDistancePointCircle(0);
            }
            else if (m2b2 > (Real)0)
            {
                s = -m2b2 - rm0;
                mNumClosestLine = 1;
                mClosestLine[0] =
                    mLine->Origin + (s + lambda)*mLine->Direction;
                sqrDistance = SqrDistancePointCircle(0);
            }
            else
            {
                s = -m2b2 + rm0;
                mClosestLine[0] =
                    mLine->Origin + (s + lambda)*mLine->Direction;
                sqrDistance = SqrDistancePointCircle(0);
                s = -m2b2 - rm0;
                mClosestLine[1] =
                    mLine->Origin + (s + lambda)*mLine->Direction;
                sqrDistance = SqrDistancePointCircle(1);
                mNumClosestLine = 2;
                mNumClosestCircle = 2;
            }
        }
    }
    else
    {
        mNumClosestLine = 1;
        mClosestLine[0] =
            mLine->Origin - mCircle->Normal.Dot(D)*mCircle->Normal;
        sqrDistance = SqrDistancePointCircle(0);
    }

    mClosestPoint0 = mClosestLine[0];
    mClosestPoint1 = mClosestCircle[0];
    mHasMultipleClosestPoints0 = (mNumClosestLine != 1);
    mHasMultipleClosestPoints1 = (mNumClosestCircle != 1);

    return sqrDistance;
}
//----------------------------------------------------------------------------
template <typename Real>
int DistLine3Circle3<Real>::GetNumClosestLine () const
{
    return mNumClosestLine;
}
//----------------------------------------------------------------------------
template <typename Real>
const Vector3<Real>& DistLine3Circle3<Real>::GetClosestLine (int i) const
{
    return mClosestLine[i];
}
//----------------------------------------------------------------------------
template <typename Real>
int DistLine3Circle3<Real>::GetNumClosestCircle () const
{
    return mNumClosestCircle;
}
//----------------------------------------------------------------------------
template <typename Real>
const Vector3<Real>& DistLine3Circle3<Real>::GetClosestCircle (int i) const
{
    return mClosestCircle[i];
}
//----------------------------------------------------------------------------
template <typename Real>
Real DistLine3Circle3<Real>::Get (Real t,
    const Vector3<Real>& velocity0, const Vector3<Real>& velocity1)
{
    Vector3<Real> movedOrigin = mLine->Origin + t*velocity0;
    Vector3<Real> movedCenter = mCircle->Center + t*velocity1;
    Line3<Real> movedLine(movedOrigin, mLine->Direction);
    Circle3<Real> movedCircle(movedCenter, mCircle->Direction0,
        mCircle->Direction1, mCircle->Normal, mCircle->Radius);
    return DistLine3Circle3<Real>(movedLine, movedCircle).Get();
}
//----------------------------------------------------------------------------
template <typename Real>
Real DistLine3Circle3<Real>::GetSquared (Real t,
    const Vector3<Real>& velocity0, const Vector3<Real>& velocity1)
{
    Vector3<Real> movedOrigin = mLine->Origin + t*velocity0;
    Vector3<Real> movedCenter = mCircle->Center + t*velocity1;
    Line3<Real> movedLine(movedOrigin, mLine->Direction);
    Circle3<Real> movedCircle(movedCenter, mCircle->Direction0,
        mCircle->Direction1, mCircle->Normal, mCircle->Radius);
    return DistLine3Circle3<Real>(movedLine, movedCircle).GetSquared();
}
//----------------------------------------------------------------------------
template <typename Real>
Real DistLine3Circle3<Real>::SqrDistancePointCircle (int i)
{
    Vector3<Real> PmC = mClosestLine[i] - mCircle->Center;
    Vector3<Real> QmC = PmC - mCircle->Normal.Dot(PmC)*mCircle->Normal;
    Real lengthQmC = QmC.Length();
    if (lengthQmC > (Real)0)
    {
        mNumClosestCircle = 1;
        mClosestCircle[i] = mCircle->Center + mCircle->Radius*QmC/lengthQmC;
    }
    else
    {
        // All circle points are equidistant from P.  Return one of them.
        mNumClosestCircle = INT_MAX;
        mClosestCircle[i] =
            mCircle->Center + mCircle->Radius*mCircle->Direction0;
    }
    Vector3<Real> diff = mClosestLine[i] - mClosestCircle[i];
    return diff.Dot(diff);
}
//----------------------------------------------------------------------------
template <typename Real>
Real DistLine3Circle3<Real>::BisectF (Real m2b2, Real rm0sqr, Real m0sqr,
    Real b1sqr, Real smin, Real smax)
{
    Real s = (Real)0, f = (Real)0;
    for (int i = 0; i < 1024; ++i)
    {
        s = ((Real)0.5)*(smin + smax);
        f = s + m2b2 - rm0sqr*s/sqrt(m0sqr*s*s + b1sqr);
        if (f == (Real)0 || s == smin || s == smax)
        {
            return s;
        }
        if (f > (Real)0)
        {
            smax = s;
        }
        else
        {
            smin = s;
        }
    }

    assertion(false, "Exceeded maximum iterations.");
    return s;
}
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Explicit instantiation.
//----------------------------------------------------------------------------
template WM5_MATHEMATICS_ITEM
class DistLine3Circle3<float>;

template WM5_MATHEMATICS_ITEM
class DistLine3Circle3<double>;
//----------------------------------------------------------------------------
}
