// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.1 (2011/03/27)

#ifndef WM5IVECTOR_H
#define WM5IVECTOR_H

#include "Wm5MathematicsLIB.h"
#include "Wm5Assert.h"

namespace Wm5
{

template <int VSIZE>
class IVector
{
public:
    // Construction.
    IVector ();
    IVector (const IVector& vec);

    // Coordinate access.
    inline operator const int* () const;
    inline operator int* ();
    inline const int& operator[] (int i) const;
    inline int& operator[] (int i);

    // Assignment.
    IVector& operator= (const IVector& vec);

    // Comparison.
    bool operator== (const IVector& vec) const;
    bool operator!= (const IVector& vec) const;
    bool operator<  (const IVector& vec) const;
    bool operator<= (const IVector& vec) const;
    bool operator>  (const IVector& vec) const;
    bool operator>= (const IVector& vec) const;

    // Arithmetic operations.
    IVector operator+ (const IVector& vec) const;
    IVector operator- (const IVector& vec) const;
    IVector operator* (const int& scalar) const;
    IVector operator/ (const int& scalar) const;
    IVector operator- () const;

    // Arithmetic updates.
    IVector& operator+= (const IVector& vec);
    IVector& operator-= (const IVector& vec);
    IVector& operator*= (const int& scalar);
    IVector& operator/= (const int& scalar);

    // Vector operations.
	int SquaredLength () const;
	int Dot (const IVector& vec) const;

protected:
    // Support for comparisons.
    int CompareArrays (const IVector& vec) const;

	int mTuple[VSIZE];
};

template <int VSIZE>
IVector<VSIZE> operator* (const int & scalar, const IVector<VSIZE>& vec);

#include "Wm5IVector.inl"

}

#endif
