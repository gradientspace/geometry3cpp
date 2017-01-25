// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <GTEngineDEF.h>

namespace gte
{

// A template class to provide 2D array access that conforms to row-major
// order (RowMajor = true) or column-major order (RowMajor = false).  The
template <bool RowMajor, typename Real, int... Dimensions>
class Array2 {};

// The array dimensions are known only at run time.
template <typename Real>
class Array2<true, Real>
{
public:
    inline Array2(int numRows, int numCols, Real* matrix);

    inline int GetNumRows() const;
    inline int GetNumCols() const;
    inline Real& operator()(int r, int c);
    inline Real const& operator()(int r, int c) const;

private:
    int mNumRows, mNumCols;
    Real* mMatrix;
};

template <typename Real>
class Array2<false, Real>
{
public:
    inline Array2(int numRows, int numCols, Real* matrix);

    inline int GetNumRows() const;
    inline int GetNumCols() const;
    inline Real& operator()(int r, int c);
    inline Real const& operator()(int r, int c) const;

private:
    int mNumRows, mNumCols;
    Real* mMatrix;
};

// The array dimensions are known at compile time.
template <typename Real, int NumRows, int NumCols>
class Array2<true, Real, NumRows, NumCols>
{
public:
    inline Array2(Real* matrix);

    inline int GetNumRows() const;
    inline int GetNumCols() const;
    inline Real& operator()(int r, int c);
    inline Real const& operator()(int r, int c) const;

private:
    Real* mMatrix;
};

template <typename Real, int NumRows, int NumCols>
class Array2<false, Real, NumRows, NumCols>
{
public:
    inline Array2(Real* matrix);

    inline int GetNumRows() const;
    inline int GetNumCols() const;
    inline Real& operator()(int r, int c);
    inline Real const& operator()(int r, int c) const;

private:
    Real* mMatrix;
};


template <typename Real> inline
Array2<true, Real>::Array2(int numRows, int numCols, Real* matrix)
:
mNumRows(numRows),
mNumCols(numCols),
mMatrix(matrix)
{
}

template <typename Real> inline
int Array2<true, Real>::GetNumRows() const
{
    return mNumRows;
}

template <typename Real> inline
int Array2<true, Real>::GetNumCols() const
{
    return mNumCols;
}

template <typename Real> inline
Real& Array2<true, Real>::operator()(int r, int c)
{
    return mMatrix[c + mNumCols*r];
}

template <typename Real> inline
Real const& Array2<true, Real>::operator()(int r, int c) const
{
    return mMatrix[c + mNumCols*r];
}



template <typename Real> inline
Array2<false, Real>::Array2(int numRows, int numCols, Real* matrix)
:
mNumRows(numRows),
mNumCols(numCols),
mMatrix(matrix)
{
}

template <typename Real> inline
int Array2<false, Real>::GetNumRows() const
{
    return mNumRows;
}

template <typename Real> inline
int Array2<false, Real>::GetNumCols() const
{
    return mNumCols;
}

template <typename Real> inline
Real& Array2<false, Real>::operator()(int r, int c)
{
    return mMatrix[r + mNumRows*c];
}

template <typename Real> inline
Real const& Array2<false, Real>::operator()(int r, int c) const
{
    return mMatrix[r + mNumRows*c];
}



template <typename Real, int NumRows, int NumCols> inline
Array2<true, Real, NumRows, NumCols>::Array2(Real* matrix)
:
mMatrix(matrix)
{
}

template <typename Real, int NumRows, int NumCols> inline
int Array2<true, Real, NumRows, NumCols>::GetNumRows() const
{
    return NumRows;
}

template <typename Real, int NumRows, int NumCols> inline
int Array2<true, Real, NumRows, NumCols>::GetNumCols() const
{
    return NumCols;
}

template <typename Real, int NumRows, int NumCols> inline
Real& Array2<true, Real, NumRows, NumCols>::operator()(int r, int c)
{
    return mMatrix[c + NumCols*r];
}

template <typename Real, int NumRows, int NumCols> inline
Real const& Array2<true, Real, NumRows, NumCols>::operator()(int r, int c)
const
{
    return mMatrix[c + NumCols*r];
}



template <typename Real, int NumRows, int NumCols> inline
Array2<false, Real, NumRows, NumCols>::Array2(Real* matrix)
:
mMatrix(matrix)
{
}

template <typename Real, int NumRows, int NumCols> inline
int Array2<false, Real, NumRows, NumCols>::GetNumRows() const
{
    return NumRows;
}

template <typename Real, int NumRows, int NumCols> inline
int Array2<false, Real, NumRows, NumCols>::GetNumCols() const
{
    return NumCols;
}

template <typename Real, int NumRows, int NumCols> inline
Real& Array2<false, Real, NumRows, NumCols>::operator()(int r, int c)
{
    return mMatrix[r + NumRows*c];
}

template <typename Real, int NumRows, int NumCols> inline
Real const& Array2<false, Real, NumRows, NumCols>::operator()(int r, int c)
const
{
    return mMatrix[r + NumRows*c];
}


}
