// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.0 (2010/01/01)

#include "Wm5MathematicsPCH.h"
#include "Wm5Color4b.h"
using namespace Wm5;

//----------------------------------------------------------------------------
Color4b::Color4b ()
{
    mTuple[0] = mTuple[1] = mTuple[2] = mTuple[3] = 255;
}
//----------------------------------------------------------------------------
Color4b::Color4b (const Color4b& vec)
{
    mTuple[0] = vec.mTuple[0];
    mTuple[1] = vec.mTuple[1];
    mTuple[2] = vec.mTuple[2];
	mTuple[3] = vec.mTuple[3];
}
//----------------------------------------------------------------------------
Color4b::Color4b (const BVector<4>& vec)
{
    mTuple[0] = vec[0];
    mTuple[1] = vec[1];
    mTuple[2] = vec[2];
	mTuple[3] = vec[3];
}
Color4b::Color4b( const unsigned char grey, unsigned char a )
	: Color4b(grey,grey,grey,a)
{
}
//----------------------------------------------------------------------------
Color4b::Color4b (const unsigned char & r, const unsigned char & g, const unsigned char & b, const unsigned char &a)
{
	mTuple[0] = r;
	mTuple[1] = g;
	mTuple[2] = b;
	mTuple[3] = a;
}
//----------------------------------------------------------------------------
Color4b::Color4b (const unsigned char * pValues)
{
	mTuple[0] = pValues[0];
	mTuple[1] = pValues[1];
	mTuple[2] = pValues[2];
	mTuple[3] = pValues[3];
}//----------------------------------------------------------------------------
Color4b::Color4b (const int& r, const int& g, const int& b, const int & a)
{
	mTuple[0] = clamp(r);
	mTuple[1] = clamp(g);
	mTuple[2] = clamp(b);
	mTuple[3] = clamp(a);
}
//----------------------------------------------------------------------------
Color4b::Color4b (const int * pValues)
{
	mTuple[0] = clamp(pValues[0]);
	mTuple[1] = clamp(pValues[1]);
	mTuple[2] = clamp(pValues[2]);
	mTuple[3] = clamp(pValues[3]);
}
//----------------------------------------------------------------------------
Color4b::Color4b (const float * pValues)
{
	mTuple[0] = clamp( (int)(pValues[0] * 255.0f) );
	mTuple[1] = clamp( (int)(pValues[1] * 255.0f) );
	mTuple[2] = clamp( (int)(pValues[2] * 255.0f) );
	mTuple[3] = clamp( (int)(pValues[3] * 255.0f) );
}
//----------------------------------------------------------------------------
Color4b::Color4b (const double * pValues)
{
	mTuple[0] = clamp( (int)(pValues[0] * 255.0) );
	mTuple[1] = clamp( (int)(pValues[1] * 255.0) );
	mTuple[2] = clamp( (int)(pValues[2] * 255.0) );
	mTuple[3] = clamp( (int)(pValues[3] * 255.0) );
}
//----------------------------------------------------------------------------
Color4b::Color4b( float r, float g, float b, float a )
{
	mTuple[0] = clamp( (int)(r * 255.0f) );
	mTuple[1] = clamp( (int)(g * 255.0f) );
	mTuple[2] = clamp( (int)(b * 255.0f) );
	mTuple[3] = clamp( (int)(a * 255.0f) );
}
//----------------------------------------------------------------------------
Color4b::Color4b( double r, double g, double b, double a )
{
	mTuple[0] = clamp( (int)(r * 255.0) );
	mTuple[1] = clamp( (int)(g * 255.0) );
	mTuple[2] = clamp( (int)(b * 255.0) );
	mTuple[3] = clamp( (int)(a * 255.0) );
}
//----------------------------------------------------------------------------
Color4b::Color4b( float g, float a )
	: Color4b(g,g,g,a)
{
}
//----------------------------------------------------------------------------
Color4b::Color4b (double g, double a)
	: Color4b(g,g,g,a)
{
}
//----------------------------------------------------------------------------
Color4b& Color4b::operator= (const Color4b& vec)
{
    mTuple[0] = vec.mTuple[0];
    mTuple[1] = vec.mTuple[1];
    mTuple[2] = vec.mTuple[2];
	mTuple[3] = vec.mTuple[3];
	return *this;
}
//----------------------------------------------------------------------------
Color4b& Color4b::operator= (const BVector<4>& vec)
{
    mTuple[0] = vec[0];
    mTuple[1] = vec[1];
    mTuple[2] = vec[2];
	mTuple[3] = vec[3];
	return *this;
}
//----------------------------------------------------------------------------
Color4b& Color4b::operator= (const Vector4f& vec)
{
	mTuple[0] = clamp( (int)(vec[0] * 255.0f) );
	mTuple[1] = clamp( (int)(vec[1] * 255.0f) );
	mTuple[2] = clamp( (int)(vec[2] * 255.0f) );
	mTuple[3] = clamp( (int)(vec[3] * 255.0f) );
	return *this;
}
//----------------------------------------------------------------------------
Color4b& Color4b::operator= (const Vector4d& vec)
{
	mTuple[0] = clamp( (int)(vec[0] * 255.0) );
	mTuple[1] = clamp( (int)(vec[1] * 255.0) );
	mTuple[2] = clamp( (int)(vec[2] * 255.0) );
	mTuple[3] = clamp( (int)(vec[3] * 255.0) );
	return *this;
}
//----------------------------------------------------------------------------
int Color4b::Dot (const Color4b& vec) const
{
    return (int)mTuple[0]*(int)vec.mTuple[0] + (int)mTuple[1]*(int)vec.mTuple[1] +
		(int)mTuple[2]*(int)vec.mTuple[2] + (int)mTuple[3]*(int)vec.mTuple[3];
}
//----------------------------------------------------------------------------
Vector4f Color4b::ToFloat() const
{
	return Vector4f( (float)mTuple[0]/255.0f, (float)mTuple[1]/255.0f, (float)mTuple[2]/255.0f, (float)mTuple[3]/255.0f );
}
//----------------------------------------------------------------------------
Vector4d Color4b::ToDouble() const
{
	return Vector4d( (float)mTuple[0]/255.0f, (float)mTuple[1]/255.0f, (float)mTuple[2]/255.0f, (float)mTuple[3]/255.0f );
}
//----------------------------------------------------------------------------
Vector3f Color4b::ToFloat3() const
{
	return Vector3f( (float)mTuple[0]/255.0f, (float)mTuple[1]/255.0f, (float)mTuple[2]/255.0f );
}
//----------------------------------------------------------------------------
Vector3d Color4b::ToDouble3() const
{
	return Vector3d( (float)mTuple[0]/255.0f, (float)mTuple[1]/255.0f, (float)mTuple[2]/255.0f );
}



static const char * Wm5Color4b_HexTable[256] = 
{
	"00", "01", "02", "03", "04", "05", "06", "07", "08", "09", "0a", "0b", "0c", "0d", "0e", "0f", "10", "11",
	"12", "13", "14", "15", "16", "17", "18", "19", "1a", "1b", "1c", "1d", "1e", "1f", "20", "21", "22", "23",
	"24", "25", "26", "27", "28", "29", "2a", "2b", "2c", "2d", "2e", "2f", "30", "31", "32", "33", "34", "35",
	"36", "37", "38", "39", "3a", "3b", "3c", "3d", "3e", "3f", "40", "41", "42", "43", "44", "45", "46", "47",
	"48", "49", "4a", "4b", "4c", "4d", "4e", "4f", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59",
	"5a", "5b", "5c", "5d", "5e", "5f", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "6a", "6b",
	"6c", "6d", "6e", "6f", "70", "71", "72", "73", "74", "75", "76", "77", "78", "79", "7a", "7b", "7c", "7d",
	"7e", "7f", "80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "8a", "8b", "8c", "8d", "8e", "8f",
	"90", "91", "92", "93", "94", "95", "96", "97", "98", "99", "9a", "9b", "9c", "9d", "9e", "9f", "a0", "a1",
	"a2", "a3", "a4", "a5", "a6", "a7", "a8", "a9", "aa", "ab", "ac", "ad", "ae", "af", "b0", "b1", "b2", "b3",
	"b4", "b5", "b6", "b7", "b8", "b9", "ba", "bb", "bc", "bd", "be", "bf", "c0", "c1", "c2", "c3", "c4", "c5",
	"c6", "c7", "c8", "c9", "ca", "cb", "cc", "cd", "ce", "cf", "d0", "d1", "d2", "d3", "d4", "d5", "d6", "d7",
	"d8", "d9", "da", "db", "dc", "dd", "de", "df", "e0", "e1", "e2", "e3", "e4", "e5", "e6", "e7", "e8", "e9",
	"ea", "eb", "ec", "ed", "ee", "ef", "f0", "f1", "f2", "f3", "f4", "f5", "f6", "f7", "f8", "f9", "fa", "fb",
	"fc", "fd", "fe", "ff"
};


HexColor Color4b::ToHexadecimal() const
{
	HexColor s;
	for (int j = 0; j < 4; ++j) {
		s.s[2*j] = Wm5Color4b_HexTable[mTuple[j]][0];
		s.s[2*j+1] = Wm5Color4b_HexTable[mTuple[j]][1];
	}
	s.s[8] = '\0';
	return s;
}



const Color4b::ColorMap Color4b::RGBA = { {0,1,2,3} };
const Color4b::ColorMap Color4b::BGRA = { {2,1,0,3} };
const Color4b::ColorMap Color4b::ARGB = { {3,0,1,2} };
Color4b Color4b::Remap( const ColorMap & from, const ColorMap & to )
{
	unsigned char c[4];
	c[from.m[0]] = mTuple[0];
	c[from.m[1]] = mTuple[1];
	c[from.m[2]] = mTuple[2];
	c[from.m[3]] = mTuple[3];

	//Color4b out;
	//out.mTuple[0] = c.mTuple[to.m[0]];
	//out.mTuple[1] = c.mTuple[to.m[1]];
	//out.mTuple[2] = c.mTuple[to.m[2]];
	//out.mTuple[3] = c.mTuple[to.m[3]];
	return Color4b( c[to.m[0]], c[to.m[1]], c[to.m[2]], c[to.m[3]] );
}



Color4b Color4b::Lerp( const Color4b & c0, const Color4b & c1, float f )
{
	Vector4f vc = (1.0f-f)*c0.ToFloat() + (f)*c1.ToFloat();
	return Color4b(vc);
}