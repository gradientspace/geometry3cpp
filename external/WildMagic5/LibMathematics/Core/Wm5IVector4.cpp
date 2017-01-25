// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.0 (2010/01/01)

#include "Wm5MathematicsPCH.h"
#include "Wm5IVector4.h"
using namespace Wm5;

//----------------------------------------------------------------------------
IVector4::IVector4 ()
{
    // uninitialized
}
//----------------------------------------------------------------------------
IVector4::IVector4 (const IVector4& vec)
{
    mTuple[0] = vec.mTuple[0];
    mTuple[1] = vec.mTuple[1];
    mTuple[2] = vec.mTuple[2];
	mTuple[3] = vec.mTuple[3];
}
//----------------------------------------------------------------------------
IVector4::IVector4 (const IVector<4>& vec)
{
    mTuple[0] = vec[0];
    mTuple[1] = vec[1];
    mTuple[2] = vec[2];
	mTuple[3] = vec[3];
}
//----------------------------------------------------------------------------
IVector4::IVector4 (const int& x, const int& y, const int& z, const int & w)
{
    mTuple[0] = x;
    mTuple[1] = y;
    mTuple[2] = z;
	mTuple[3] = w;
}
//----------------------------------------------------------------------------
IVector4::IVector4 (const int * pValues)
{
	mTuple[0] = pValues[0];
	mTuple[1] = pValues[1];
	mTuple[2] = pValues[2];
	mTuple[3] = pValues[3];
}
//----------------------------------------------------------------------------
IVector4& IVector4::operator= (const IVector4& vec)
{
    mTuple[0] = vec.mTuple[0];
    mTuple[1] = vec.mTuple[1];
    mTuple[2] = vec.mTuple[2];
	mTuple[3] = vec.mTuple[3];
	return *this;
}
//----------------------------------------------------------------------------
IVector4& IVector4::operator= (const IVector<4>& vec)
{
    mTuple[0] = vec[0];
    mTuple[1] = vec[1];
    mTuple[2] = vec[2];
	mTuple[3] = vec[3];
	return *this;
}
//----------------------------------------------------------------------------
int IVector4::Dot (const IVector4& vec) const
{
    return mTuple[0]*vec.mTuple[0] + mTuple[1]*vec.mTuple[1] +
        mTuple[2]*vec.mTuple[2] + mTuple[3]*vec.mTuple[3];
}
