// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.1 (2010/10/01)

#ifndef WM5COLOR4B_H
#define WM5COLOR4B_H

#include "Wm5MathematicsLIB.h"
#include "Wm5BVector.h"
#include "Wm5Vector4.h"

namespace Wm5
{

// like #AABBCCFF, but without the leading "#"
struct WM5_MATHEMATICS_ITEM HexColor
{
	char s[9];
};


class WM5_MATHEMATICS_ITEM Color4b : public BVector<4>
{
public:
	// [TODO] support color maps as arguments to various functions...eg RGBA, BGRA, ARGB
	//    could declare these as int[4] static members, then we can use directly as swizzlers

	// [TODO] move .cpp functions to an .inl file so they can be inlined...

    // Construction.
    Color4b ();
    Color4b (const Color4b& vec);
    Color4b (const BVector<4>& vec);
	Color4b (const unsigned char & r, const unsigned char & g, const unsigned char & b, const unsigned char &a = 255);
	Color4b (const unsigned char * pValues);
	Color4b (unsigned char grey, unsigned char a = 255);
	Color4b (const int& r, const int& g, const int& b, const int &a);
	Color4b (const int * pValues);
	Color4b (const float * pValues);
	Color4b (const double * pValues);
	Color4b (float r, float g, float b, float a = 1.0f);
	Color4b (double r, double g, double b, double a = 1.0);
	Color4b (float grey, float a = 1.0f);
	Color4b (double grey, double a = 1.0);

    // Assignment.
    Color4b& operator= (const Color4b& vec);
    Color4b& operator= (const BVector<4>& vec);
	Color4b& operator= (const Vector4f& vec);
	Color4b& operator= (const Vector4d& vec);

    // returns Dot(this,V)
    int Dot (const Color4b& vec) const;

	Vector4f ToFloat() const;
	Vector4d ToDouble() const;
	Vector3f ToFloat3() const;
	Vector3d ToDouble3() const;
	HexColor ToHexadecimal() const;

	struct ColorMap {
		int m[4];
	};
	static const ColorMap RGBA;
	static const ColorMap BGRA;
	static const ColorMap ARGB;
	Color4b Remap( const ColorMap & from, const ColorMap & to );

	static Color4b Lerp(const Color4b & c0, const Color4b & c1, float f);

protected:
    using BVector<4>::mTuple;
};

}

#endif
