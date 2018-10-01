// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.1 (2010/10/01)

#ifndef WM5IVECTOR2_H
#define WM5IVECTOR2_H

#include "Wm5MathematicsLIB.h"
#include "Wm5IVector.h"

namespace Wm5
{

class WM5_MATHEMATICS_ITEM IVector2 : public IVector<2>
{
public:
    // Construction.
    IVector2 ();
    IVector2 (const IVector2& vec);
    IVector2 (const IVector<2>& vec);
    IVector2 (const int& x, const int& y);

    // Member access.
    inline int X () const;
    inline int& X ();
    inline int Y () const;
    inline int& Y ();

    // Assignment,
    IVector2& operator= (const IVector2& vec);
    IVector2& operator= (const IVector<2>& vec);

    // Returns Dot(this,V).
    int Dot (const IVector2& vec) const;

    // Returns (y,-x).
    IVector2 Perp () const;

    // Returns Cross((x,y,0),(V.x,V.y,0)) = x*V.y - y*V.x.
    int DotPerp (const IVector2& vec) const;

protected:
    using IVector<2>::mTuple;


// [geometry3]
public:
	using EVector2 = Eigen::Matrix<int, 2, 1>;
	operator EVector2() const {
		return Eigen::Map<EVector2>((int *)this);
	}
	IVector2(const EVector2 & vec) {
		const int * p = vec.data();
		mTuple[0] = p[0];
		mTuple[1] = p[1];
	}
// [geometry3]

};

#include "Wm5IVector2.inl"

}

#endif
