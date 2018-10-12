// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.0 (2010/01/01)

//----------------------------------------------------------------------------
template <typename Real>
Box3<Real>::Box3 ()
{
}
//----------------------------------------------------------------------------
template <typename Real>
Box3<Real>::~Box3 ()
{
}
//----------------------------------------------------------------------------
template <typename Real>
Box3<Real>::Box3 (const Vector3<Real>& center)
	:
	Center(center)
{
	Axis[0] = Vector3<Real>::UNIT_X;
	Axis[1] = Vector3<Real>::UNIT_Y;
	Axis[2] = Vector3<Real>::UNIT_Z;
	Extent[0] = Extent[1] = Extent[2] = (Real)0;
}
//----------------------------------------------------------------------------
template <typename Real>
Box3<Real>::Box3 (const Vector3<Real>& center, const Vector3<Real> axis[3],
    const Real extent[3])
    :
    Center(center)
{
    Axis[0] = axis[0];
    Axis[1] = axis[1];
    Axis[2] = axis[2];
    Extent[0] = extent[0];
    Extent[1] = extent[1];
    Extent[2] = extent[2];
}
//----------------------------------------------------------------------------
template <typename Real>
Box3<Real>::Box3 (const Vector3<Real>& center, const Vector3<Real>& axis0,
    const Vector3<Real>& axis1, const Vector3<Real>& axis2,
    const Real extent0, const Real extent1, const Real extent2)
    :
    Center(center)
{
    Axis[0] = axis0;
    Axis[1] = axis1;
    Axis[2] = axis2;
    Extent[0] = extent0;
    Extent[1] = extent1;
    Extent[2] = extent2;
}
//====
template <typename Real>
Box3<Real>::Box3( const AxisAlignedBox3<Real>& aaBox )
{
	Center = aaBox.Center();
	Axis[0] = Vector3<Real>::UNIT_X;
	Axis[1] = Vector3<Real>::UNIT_Y;
	Axis[2] = Vector3<Real>::UNIT_Z;
	Extent[0] = aaBox.Dimension(0) * (Real)0.5;
	Extent[1] = aaBox.Dimension(1) * (Real)0.5;
	Extent[2] = aaBox.Dimension(2) * (Real)0.5;
}
//----------------------------------------------------------------------------
template <typename Real>
void Box3<Real>::ComputeVertices (Vector3<Real> vertex[8]) const
{
    Vector3<Real> extAxis0 = Extent[0]*Axis[0];
    Vector3<Real> extAxis1 = Extent[1]*Axis[1];
    Vector3<Real> extAxis2 = Extent[2]*Axis[2];

    vertex[0] = Center - extAxis0 - extAxis1 - extAxis2;
    vertex[1] = Center + extAxis0 - extAxis1 - extAxis2;
    vertex[2] = Center + extAxis0 + extAxis1 - extAxis2;
    vertex[3] = Center - extAxis0 + extAxis1 - extAxis2;
    vertex[4] = Center - extAxis0 - extAxis1 + extAxis2;
    vertex[5] = Center + extAxis0 - extAxis1 + extAxis2;
    vertex[6] = Center + extAxis0 + extAxis1 + extAxis2;
    vertex[7] = Center - extAxis0 + extAxis1 + extAxis2;
}
//----------------------------------------------------------------------------





template <typename Real>
Real Box3<Real>::MaxExtent() const
{
	return std::max(Extent[0], std::max(Extent[1], Extent[2]));
}
template <typename Real>
Real Box3<Real>::MinExtent() const
{
	return std::min(Extent[0], std::min(Extent[1], Extent[2]));
}

template <typename Real>
Vector3<Real> Box3<Real>::Diagonal() const
{
	return 
		(Extent[0]*Axis[0] + Extent[1]*Axis[1] + Extent[2]*Axis[2]) - 
		(-Extent[0]*Axis[0] - Extent[1]*Axis[1] - Extent[2]*Axis[2]);
}

template <typename Real>
Real Box3<Real>::Volume() const
{
	return ((Real)2*Extent[0]) + ((Real)2*Extent[1]) + ((Real)2*Extent[2]);
}

template <typename Real>
void Box3<Real>::Contain(const Vector3<Real> & v)
{
	Vector3<Real> lv = v - Center;
	for (int k = 0; k < 3; ++k) {
		Real t = lv.Dot(Axis[k]);
		if ( fabs(t) > Extent[k]) {
			Real min = -Extent[k], max = Extent[k];
			if ( t < min )
				min = t;
			else if ( t > max )
				max = t;
			Extent[k] = (max-min) * (Real)0.5;
			Center = Center + ((max+min) * (Real)0.5) * Axis[k];
		}
	}
}

template <typename Real>
void Box3<Real>::Contain(const Box3<Real> & o)
{
	Vector3<Real> v[8];
	o.ComputeVertices(v);
	for (int k = 0; k < 8; ++k) 
		Contain(v[k]);
}


template <typename Real>
bool Box3<Real>::Contains(const Vector3<Real> & v) const
{
	Vector3<Real> lv = v - Center;
	return (fabs(lv.Dot(Axis[0])) <= Extent[0]) &&
		   (fabs(lv.Dot(Axis[1])) <= Extent[1]) &&
		   (fabs(lv.Dot(Axis[2])) <= Extent[2]);
}


template <typename Real>
void Box3<Real>::Expand(Real f)
{
	Extent[0] += f; Extent[1] += f; Extent[2] += f;
}

template <typename Real>
void Box3<Real>::Translate(const Vector3<Real> & v)
{
	Center += v;
}
