// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.0 (2010/01/01)
//
// [RMS] 10/3/2018 - added this

#ifndef WM5INTERVAL1_H
#define WM5INTERVAL1_H

#include "Wm5MathematicsLIB.h"

namespace Wm5
{

template <typename Real>
class Interval1
{
public:
    // Construction and destruction.
    Interval1 ();  // uninitialized
    ~Interval1 ();

    // The caller must ensure that xmin <= xmax.
    Interval1 (Real xmin, Real xmax);

    // Overlap testing is in the strict sense.  If the two boxes are just
    // touching along a common edge, the boxes are reported as overlapping.
    bool Overlaps (const Interval1& other) const;

    // The return value is 'true' if there is overlap.  In this case the
    // intersection is stored in 'intersection'.  If the return value is
    // 'false', there is no overlap.  In this case 'intersection' is
    // undefined.
    Interval1 IntersectionWith (const Interval1& other) const;

	Real Min, Max;

	// g3 extensions
	Real Center() const;
	Real Length() const;
	Real LengthSquared() const;
	void Contain(const Real & v);
	void Contain(const Interval1<Real> & o);
	bool Contains(const Real & v) const;
	void Expand(Real f);
	void Translate(const Real & v);

	WM5_MATHEMATICS_ITEM static const Interval1 EMPTY;    // [inf,-inf]
};

#include "Wm5Interval1.inl"

typedef Interval1<float> Interval1f;
typedef Interval1<double> Interval1d;

}

#endif
