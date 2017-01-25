// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.0 (2010/01/01)

//----------------------------------------------------------------------------
template <typename Real>
AxisAlignedBox3<Real>::AxisAlignedBox3 ()
{
}
//----------------------------------------------------------------------------
template <typename Real>
AxisAlignedBox3<Real>::~AxisAlignedBox3 ()
{
}
//----------------------------------------------------------------------------
template <typename Real>
AxisAlignedBox3<Real>::AxisAlignedBox3 (Real xmin, Real xmax, Real ymin,
    Real ymax, Real zmin, Real zmax)
{
    Min[0] = xmin;
    Max[0] = xmax;
    Min[1] = ymin;
    Max[1] = ymax;
    Min[2] = zmin;
    Max[2] = zmax;
}
//
template <typename Real>
AxisAlignedBox3<Real>::AxisAlignedBox3(const Real * vMin, const Real * vMax)
{
	Min[0] = vMin[0]; Min[1] = vMin[1]; Min[2] = vMin[2];
	Max[0] = vMax[0]; Max[1] = vMax[1]; Max[2] = vMax[2];
}
//
template <typename Real>
AxisAlignedBox3<Real>::AxisAlignedBox3(const Vector3<Real>& vMin, const Vector3<Real>& vMax)
{
	Min[0] = vMin[0]; Min[1] = vMin[1]; Min[2] = vMin[2];
	Max[0] = vMax[0]; Max[1] = vMax[1]; Max[2] = vMax[2];
}
//
template <typename Real>
AxisAlignedBox3<Real>::AxisAlignedBox3(const Vector3<Real>& vCenter, Real fRadius)
{
	Min[0] = vCenter.X() - fRadius;
	Max[0] = vCenter.X() + fRadius;
	Min[1] = vCenter.Y() - fRadius;
	Max[1] = vCenter.Y() + fRadius;
	Min[2] = vCenter.Z() - fRadius;
	Max[2] = vCenter.Z() + fRadius;
}
//----------------------------------------------------------------------------
template <typename Real>
void AxisAlignedBox3<Real>::GetCenterExtents (Vector3<Real>& center,
    Real extent[3])
{
    center[0] = ((Real)0.5)*(Max[0] + Min[0]);
    center[1] = ((Real)0.5)*(Max[1] + Min[1]);
    center[2] = ((Real)0.5)*(Max[2] + Min[2]);
    extent[0] = ((Real)0.5)*(Max[0] - Min[0]);
    extent[1] = ((Real)0.5)*(Max[1] - Min[1]);
    extent[2] = ((Real)0.5)*(Max[2] - Min[2]);
}
//----------------------------------------------------------------------------
template <typename Real>
bool AxisAlignedBox3<Real>::HasXOverlap (const AxisAlignedBox3& box) const
{
    return (Max[0] >= box.Min[0] && Min[0] <= box.Max[0]);
}
//----------------------------------------------------------------------------
template <typename Real>
bool AxisAlignedBox3<Real>::HasYOverlap (const AxisAlignedBox3& box) const
{
    return (Max[1] >= box.Min[1] && Min[1] <= box.Max[1]);
}
//----------------------------------------------------------------------------
template <typename Real>
bool AxisAlignedBox3<Real>::HasZOverlap (const AxisAlignedBox3& box) const
{
    return (Max[2] >= box.Min[2] && Min[2] <= box.Max[2]);
}
//----------------------------------------------------------------------------
template <typename Real>
bool AxisAlignedBox3<Real>::TestIntersection (const AxisAlignedBox3& box)
    const
{
    for (int i = 0; i < 3; i++)
    {
        if (Max[i] < box.Min[i] || Min[i] > box.Max[i])
        {
            return false;
        }
    }
    return true;
}
//----------------------------------------------------------------------------
template <typename Real>
bool AxisAlignedBox3<Real>::FindIntersection (const AxisAlignedBox3& box,
    AxisAlignedBox3& intersection) const
{
    int i;
    for (i = 0; i < 3; i++)
    {
        if (Max[i] < box.Min[i] || Min[i] > box.Max[i])
        {
            return false;
        }
    }

    for (i = 0; i < 3; i++)
    {
        if (Max[i] <= box.Max[i])
        {
            intersection.Max[i] = Max[i];
        }
        else
        {
            intersection.Max[i] = box.Max[i];
        }

        if (Min[i] <= box.Min[i])
        {
            intersection.Min[i] = box.Min[i];
        }
        else
        {
            intersection.Min[i] = Min[i];
        }
    }
    return true;
}
//----------------------------------------------------------------------------



template <typename Real>
Real AxisAlignedBox3<Real>::Dimension(int i) const
{
	return Max[i] - Min[i];
}
template <typename Real>
Real AxisAlignedBox3<Real>::MaxDimension() const
{
	return std::max(Max[0] - Min[0], std::max(Max[1] - Min[1], Max[2] - Min[2]));
}
template <typename Real>
Real AxisAlignedBox3<Real>::MinDimension() const
{
	return std::min(Max[0] - Min[0], std::min(Max[1] - Min[1], Max[2] - Min[2]));
}

template <typename Real>
Vector3<Real> AxisAlignedBox3<Real>::Center() const
{
	return Vector3<Real>((Max[0] + Min[0])*(Real)0.5, (Max[1] + Min[1])*(Real)0.5, (Max[2] + Min[2])*(Real)0.5);
}

template <typename Real>
Vector3<Real> AxisAlignedBox3<Real>::Diagonal() const
{
	return Vector3<Real>(Max[0] - Min[0], Max[1] - Min[1], Max[2] - Min[2]);
}

template <typename Real>
Real AxisAlignedBox3<Real>::Volume() const
{
	return (Max[0] - Min[0]) * (Max[1] - Min[1]) * (Max[2] - Min[2]);
}

template <typename Real>
void AxisAlignedBox3<Real>::Contain(const Vector3<Real> & v)
{
	for (int k = 0; k < 3; ++k) {
		if (v[k] < Min[k])
			Min[k] = v[k];
		if (v[k] > Max[k])
			Max[k] = v[k];
	}
}

template <typename Real>
void AxisAlignedBox3<Real>::Contain(const AxisAlignedBox3<Real> & o)
{
	for (int k = 0; k < 3; ++k) {
		if (o.Min[k] < Min[k])
			Min[k] = o.Min[k];
		if (o.Max[k] > Max[k])
			Max[k] = o.Max[k];
	}
}


template <typename Real>
bool AxisAlignedBox3<Real>::Contained(const Vector3<Real> & v) const
{
	return (v[0] >= Min[0] && v[0] <= Max[0]
		 && v[1] >= Min[1] && v[1] <= Max[1]
		 && v[2] >= Min[2] && v[2] <= Max[2]);
}


template <typename Real>
void AxisAlignedBox3<Real>::Expand(Real f)
{
	Min[0] -= f; Min[1] -= f; Min[2] -= f;
	Max[0] += f; Max[1] += f; Max[2] += f;
}

template <typename Real>
void AxisAlignedBox3<Real>::Translate(const Vector3<Real> & v)
{
	Min[0] += v[0]; Min[1] += v[1]; Min[2] += v[2];
	Max[0] += v[0]; Max[1] += v[1]; Max[2] += v[2];
}