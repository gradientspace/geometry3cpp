// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.1 (2011/03/27)

//----------------------------------------------------------------------------
template <int VSIZE>
BVector<VSIZE>::BVector ()
{
    // For efficiency in construction of large arrays of vectors, the
    // default constructor does not initialize the vector.
}
//----------------------------------------------------------------------------
template <int VSIZE>
BVector<VSIZE>::BVector (const BVector& vec)
{
    for (int i = 0; i < VSIZE; ++i)
    {
        mTuple[i] = vec.mTuple[i];
    }
}
//----------------------------------------------------------------------------
template <int VSIZE>
inline BVector<VSIZE>::operator const unsigned char * () const
{
    return mTuple;
}
//----------------------------------------------------------------------------
template <int VSIZE>
inline BVector<VSIZE>::operator unsigned char * ()
{
    return mTuple;
}
//----------------------------------------------------------------------------
template <int VSIZE>
inline const unsigned char& BVector<VSIZE>::operator[] (int i) const
{
    assertion(0 <= i && i < VSIZE, "Invalid index\n");
    return mTuple[i];
}
//----------------------------------------------------------------------------
template <int VSIZE>
inline unsigned char& BVector<VSIZE>::operator[] (int i)
{
    assertion(0 <= i && i < VSIZE, "Invalid index\n");
    return mTuple[i];
}
//----------------------------------------------------------------------------
template <int VSIZE>
BVector<VSIZE>& BVector<VSIZE>::operator= (const BVector& vec)
{
    for (int i = 0; i < VSIZE; ++i)
    {
        mTuple[i] = vec.mTuple[i];
    }
    return *this;
}
//----------------------------------------------------------------------------
template <int VSIZE>
bool BVector<VSIZE>::operator== (const BVector& vec) const
{
    for (int i = 0; i < VSIZE; ++i)
    {
        if (mTuple[i] != vec.mTuple[i])
        {
            return false;
        }
    }
    return true;
}
//----------------------------------------------------------------------------
template <int VSIZE>
bool BVector<VSIZE>::operator!= (const BVector& vec) const
{
    return !operator==(vec);
}
//----------------------------------------------------------------------------
template <int VSIZE>
int BVector<VSIZE>::CompareArrays (const BVector& vec) const
{
    for (int i = 0; i < VSIZE; ++i)
    {
        if (mTuple[i] < vec.mTuple[i])
        {
            return -1;
        }
        if (mTuple[i] > vec.mTuple[i])
        {
            return +1;
        }
    }
    return 0;
}
//----------------------------------------------------------------------------
template <int VSIZE>
bool BVector<VSIZE>::operator< (const BVector& vec) const
{
    return CompareArrays(vec) < 0;
}
//----------------------------------------------------------------------------
template <int VSIZE>
bool BVector<VSIZE>::operator<= (const BVector& vec) const
{
    return CompareArrays(vec) <= 0;
}
//----------------------------------------------------------------------------
template <int VSIZE>
bool BVector<VSIZE>::operator> (const BVector& vec) const
{
    return CompareArrays(vec) > 0;
}
//----------------------------------------------------------------------------
template <int VSIZE>
bool BVector<VSIZE>::operator>= (const BVector& vec) const
{
    return CompareArrays(vec) >= 0;
}
//----------------------------------------------------------------------------
template <int VSIZE>
unsigned char BVector<VSIZE>::clamp( int n )
{
	return n < 0 ? 0 : (n > 255 ? 255 : n);
}
//----------------------------------------------------------------------------
template <int VSIZE>
BVector<VSIZE> BVector<VSIZE>::operator+ (const BVector& vec) const
{
    BVector<VSIZE> iSum;
    for (int i = 0; i < VSIZE; ++i)
    {
        iSum.mTuple[i] = clamp((int)mTuple[i] + (int)vec.mTuple[i]);
    }
    return iSum;
}
//----------------------------------------------------------------------------
template <int VSIZE>
BVector<VSIZE> BVector<VSIZE>::operator- (const BVector& vec) const
{
    BVector<VSIZE> diff;
    for (int i = 0; i < VSIZE; ++i)
    {
        diff.mTuple[i] = clamp((int)mTuple[i] - (int)vec.mTuple[i]);
    }
    return diff;
}
//----------------------------------------------------------------------------
template <int VSIZE>
BVector<VSIZE> BVector<VSIZE>::operator* (const int& scalar) const
{
    BVector<VSIZE> prod;
    for (int i = 0; i < VSIZE; ++i)
    {
        prod.mTuple[i] = clamp(scalar*(int)mTuple[i]);
    }
    return prod;
}
//----------------------------------------------------------------------------
template <int VSIZE>
BVector<VSIZE> BVector<VSIZE>::operator/ (const int& scalar) const
{
    assertion(scalar != 0, "Division by zero\n");

    BVector<VSIZE> div;
    for (int i = 0; i < VSIZE; ++i)
    {
        div.mTuple[i] = clamp((int)mTuple[i]/scalar);
    }

    return div;
}
//----------------------------------------------------------------------------
template <int VSIZE>
BVector<VSIZE> operator* (const int& scalar, const BVector<VSIZE>& vec)
{
    BVector<VSIZE> prod;
    for (int i = 0; i < VSIZE; ++i)
    {
        prod[i] = BVector<VSIZE>::clamp(scalar*(int)vec[i]);
    }
    return prod;
}
//----------------------------------------------------------------------------
template <int VSIZE>
BVector<VSIZE>& BVector<VSIZE>::operator+= (const BVector& vec)
{
    for (int i = 0; i < VSIZE; ++i)
    {
		mTuple[i] = clamp((int)mTuple[i] + (int)vec.mTuple[i]);
    }
    return *this;
}
//----------------------------------------------------------------------------
template <int VSIZE>
BVector<VSIZE>& BVector<VSIZE>::operator-= (const BVector& vec)
{
    for (int i = 0; i < VSIZE; ++i)
    {
		mTuple[i] = clamp((int)mTuple[i] - (int)vec.mTuple[i]);
    }
    return *this;
}
//----------------------------------------------------------------------------
template <int VSIZE>
BVector<VSIZE>& BVector<VSIZE>::operator*= (const int& scalar)
{
    for (int i = 0; i < VSIZE; ++i)
    {
		mTuple[i] = clamp((int)mTuple[i] * scalar);
    }
    return *this;
}
//----------------------------------------------------------------------------
template <int VSIZE>
BVector<VSIZE>& BVector<VSIZE>::operator/= (const int& scalar)
{
    assertion(scalar != 0, "Division by zero\n");

    for (int i = 0; i < VSIZE; ++i)
    {
		mTuple[i] = clamp((int)mTuple[i] / scalar);
    }
    return *this;
}
//----------------------------------------------------------------------------
template <int VSIZE>
int BVector<VSIZE>::SquaredLength () const
{
    int sqrLen = 0;
    for (int i = 0; i < VSIZE; ++i)
    {
        sqrLen += (int)mTuple[i] * (int)mTuple[i];
    }
    return sqrLen;
}
//----------------------------------------------------------------------------
template <int VSIZE>
int BVector<VSIZE>::Dot (const BVector& vec) const
{
    int dot = 0;
    for (int i = 0; i < VSIZE; ++i)
    {
        dot += (int)mTuple[i] * (int)vec.mTuple[i];
    }
    return dot;
}
//----------------------------------------------------------------------------
