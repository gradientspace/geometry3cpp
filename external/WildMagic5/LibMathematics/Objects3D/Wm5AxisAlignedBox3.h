// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.0 (2010/01/01)

#ifndef WM5AXISALIGNEDBOX3_H
#define WM5AXISALIGNEDBOX3_H

#include "Wm5MathematicsLIB.h"
#include "Wm5Vector3.h"

namespace Wm5
{

template <typename Real>
class AxisAlignedBox3
{
public:
    // Construction and destruction.
    AxisAlignedBox3 ();  // uninitialized
    ~AxisAlignedBox3 ();

    // The caller must ensure that xmin <= xmax, ymin <= ymax, and
    // zmin <= zmax.
private:
	// [RMSG3] disabling this constructor because order of args is different than
	// in g3Sharp. Use Vector3 version instead
    AxisAlignedBox3 (Real xmin, Real xmax, Real ymin, Real ymax,
        Real zmin, Real zmax);
public:

	AxisAlignedBox3(const Real pMin[3], const Real pMax[3]);

	AxisAlignedBox3(const Vector3<Real>& vMin, const Vector3<Real>& vMax );

	AxisAlignedBox3 ( const Vector3<Real>& center, Real extent);

    // Compute the center of the box and the extents (half-lengths)
    // of the box edges.
    void GetCenterExtents (Vector3<Real>& center, Real extent[3]);

    // Overlap testing is in the strict sense.  If the two boxes are just
    // touching along a common edge or a common face, the boxes are reported
    // as overlapping.
    bool HasXOverlap (const AxisAlignedBox3& box) const;
    bool HasYOverlap (const AxisAlignedBox3& box) const;
    bool HasZOverlap (const AxisAlignedBox3& box) const;
    bool TestIntersection (const AxisAlignedBox3& box) const;

    // The return value is 'true' if there is overlap.  In this case the
    // intersection is stored in 'intersection'.  If the return value is
    // 'false', there is no overlap.  In this case 'intersection' is
    // undefined.
    bool FindIntersection (const AxisAlignedBox3& box,
        AxisAlignedBox3& intersection) const;

    Real Min[3], Max[3];


	// g3 extensions
	Real Dimension(int i) const;
	Real MaxDimension() const;
	Real MinDimension() const;
	Vector3<Real> Center() const;
	Vector3<Real> Diagonal() const;
	Vector3<Real> Extents() const;
	Real Volume() const;
	void Contain(const Vector3<Real> & v);
	void Contain(const AxisAlignedBox3<Real> & o);
	bool Contains(const Vector3<Real> & v) const;
	void Expand(Real f);
	void Translate(const Vector3<Real> & v);


	WM5_MATHEMATICS_ITEM static const AxisAlignedBox3 EMPTY;    // [-inf,inf]
};

#include "Wm5AxisAlignedBox3.inl"

typedef AxisAlignedBox3<float> AxisAlignedBox3f;
typedef AxisAlignedBox3<double> AxisAlignedBox3d;

}

#endif
