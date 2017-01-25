// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.2 (2013/07/28)

#include "Wm5MathematicsPCH.h"
#include "Wm5DistPoint3Circle3.h"

namespace Wm5
{
//----------------------------------------------------------------------------
template <typename Real>
DistPoint3Circle3<Real>::DistPoint3Circle3 (const Vector3<Real>& point,
    const Circle3<Real>& circle)
    :
    mPoint(&point),
    mCircle(&circle)
{
}
//----------------------------------------------------------------------------
template <typename Real>
const Vector3<Real>& DistPoint3Circle3<Real>::GetPoint () const
{
    return *mPoint;
}
//----------------------------------------------------------------------------
template <typename Real>
const Circle3<Real>& DistPoint3Circle3<Real>::GetCircle () const
{
    return *mCircle;
}
//----------------------------------------------------------------------------
template <typename Real>
Real DistPoint3Circle3<Real>::Get ()
{
    return Math<Real>::Sqrt(GetSquared());
}
//----------------------------------------------------------------------------
template <typename Real>
Real DistPoint3Circle3<Real>::GetSquared ()
{
    // Projection of P-C onto plane is Q-C = P-C - Dot(N,P-C)*N.
    Vector3<Real> PmC = *mPoint - mCircle->Center;
    Vector3<Real> QmC = PmC - mCircle->Normal.Dot(PmC)*mCircle->Normal;
    Real lengthQmC = QmC.Length();
    if (lengthQmC > (Real)0)
    {
        mClosestPoint1 = mCircle->Center + mCircle->Radius*QmC/lengthQmC;
        mHasMultipleClosestPoints1 = false;
    }
    else
    {
        // All circle points are equidistant from P.  Return one of them.
        mClosestPoint1 = mCircle->Center + mCircle->Radius*mCircle->Direction0;
        mHasMultipleClosestPoints1 = true;
    }

    mClosestPoint0 = *mPoint;
    mHasMultipleClosestPoints0 = false;

    Vector3<Real> diff = mClosestPoint0 - mClosestPoint1;
    return diff.Dot(diff);
}
//----------------------------------------------------------------------------
template <typename Real>
Real DistPoint3Circle3<Real>::Get (Real t,
    const Vector3<Real>& velocity0, const Vector3<Real>& velocity1)
{
    Vector3<Real> movedPoint = *mPoint + t*velocity0;
    Vector3<Real> movedCenter = mCircle->Center + t*velocity1;
    Circle3<Real> movedCircle(movedCenter, mCircle->Direction0,
        mCircle->Direction1, mCircle->Normal, mCircle->Radius);
    return DistPoint3Circle3<Real>(movedPoint, movedCircle).Get();
}
//----------------------------------------------------------------------------
template <typename Real>
Real DistPoint3Circle3<Real>::GetSquared (Real t,
    const Vector3<Real>& velocity0, const Vector3<Real>& velocity1)
{
    Vector3<Real> movedPoint = *mPoint + t*velocity0;
    Vector3<Real> movedCenter = mCircle->Center + t*velocity1;
    Circle3<Real> movedCircle(movedCenter, mCircle->Direction0,
        mCircle->Direction1, mCircle->Normal, mCircle->Radius);
    return DistPoint3Circle3<Real>(movedPoint, movedCircle).GetSquared();
}
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Explicit instantiation.
//----------------------------------------------------------------------------
template WM5_MATHEMATICS_ITEM
class DistPoint3Circle3<float>;

template WM5_MATHEMATICS_ITEM
class DistPoint3Circle3<double>;
//----------------------------------------------------------------------------
}
