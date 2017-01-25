// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.1 (2011/03/27)

#ifndef WM5BVECTOR_H
#define WM5BVECTOR_H

#include "Wm5MathematicsLIB.h"
#include "Wm5Assert.h"

namespace Wm5
{

template <int VSIZE>
class BVector
{
public:
    // Construction.
    BVector ();
    BVector (const BVector& vec);

    // Coordinate access.
    inline operator const unsigned char* () const;
    inline operator unsigned char* ();
    inline const unsigned char& operator[] (int i) const;
    inline unsigned char& operator[] (int i);

    // Assignment.
    BVector& operator= (const BVector& vec);

    // Comparison.
    bool operator== (const BVector& vec) const;
    bool operator!= (const BVector& vec) const;
    bool operator<  (const BVector& vec) const;
    bool operator<= (const BVector& vec) const;
    bool operator>  (const BVector& vec) const;
    bool operator>= (const BVector& vec) const;

    // Arithmetic operations.
    BVector operator+ (const BVector& vec) const;
    BVector operator- (const BVector& vec) const;
    BVector operator* (const int& scalar) const;
    BVector operator/ (const int& scalar) const;

    // Arithmetic updates.
    BVector& operator+= (const BVector& vec);
    BVector& operator-= (const BVector& vec);
    BVector& operator*= (const int& scalar);
    BVector& operator/= (const int& scalar);

    // Vector operations.
	int SquaredLength () const;
	int Dot (const BVector& vec) const;

protected:
    // Support for comparisons.
    int CompareArrays (const BVector& vec) const;
	static unsigned char clamp(int n);

	unsigned char mTuple[VSIZE];
};

template <int VSIZE>
BVector<VSIZE> operator* (const int & scalar, const BVector<VSIZE>& vec);

#include "Wm5BVector.inl"

}

#endif
