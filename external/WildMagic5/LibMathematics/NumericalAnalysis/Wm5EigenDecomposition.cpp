// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.4 (2014/07/09)

#include "Wm5MathematicsPCH.h"
#include "Wm5EigenDecomposition.h"

namespace Wm5
{
//----------------------------------------------------------------------------
template <typename Real>
EigenDecomposition<Real>::EigenDecomposition (int size)
    :
    mMatrix(size, size),
    mSolver(size, MAX_ITERATIONS),
    mEigenvalues(size)
{
    assertion(size >= 2, "Invalid size in Eigendecomposition constructor\n");
}
//----------------------------------------------------------------------------
template <typename Real>
EigenDecomposition<Real>::EigenDecomposition (const Matrix2<Real>& mat)
    :
    mMatrix(2, 2, &mat[0][0]),
    mSolver(2, MAX_ITERATIONS),
    mEigenvalues(2)
{
}
//----------------------------------------------------------------------------
template <typename Real>
EigenDecomposition<Real>::EigenDecomposition (const Matrix3<Real>& mat)
    :
    mMatrix(3, 3, &mat[0][0]),
    mSolver(3, MAX_ITERATIONS),
    mEigenvalues(3)
{
}
//----------------------------------------------------------------------------
template <typename Real>
EigenDecomposition<Real>::EigenDecomposition (const GMatrix<Real>& mat)
    :
    mMatrix(mat),
    mSolver(mat.GetNumRows(), MAX_ITERATIONS),
    mEigenvalues(mat.GetNumRows())
{
    assertion(mat.GetNumRows() >= 2
        && mat.GetNumRows() == mat.GetNumColumns(),
        "Square matrix required in EigenDecomposition constructor\n");
}
//----------------------------------------------------------------------------
template <typename Real>
EigenDecomposition<Real>::~EigenDecomposition ()
{
}
//----------------------------------------------------------------------------
template <typename Real>
Real& EigenDecomposition<Real>::operator() (int row, int column)
{
    return mMatrix[row][column];
}
//----------------------------------------------------------------------------
template <typename Real>
EigenDecomposition<Real>& EigenDecomposition<Real>::operator= (
    const Matrix2<Real>& mat)
{
    mMatrix.SetMatrix(2, 2, &mat[0][0]);
    mSolver = SymmetricEigensolverGTE<Real>(2, MAX_ITERATIONS);
    mEigenvalues.resize(2);
    return *this;
}
//----------------------------------------------------------------------------
template <typename Real>
EigenDecomposition<Real>& EigenDecomposition<Real>::operator= (
    const Matrix3<Real>& mat)
{
    mMatrix.SetMatrix(3, 3, &mat[0][0]);
    mSolver = SymmetricEigensolverGTE<Real>(3, MAX_ITERATIONS);
    mEigenvalues.resize(3);
    return *this;
}
//----------------------------------------------------------------------------
template <typename Real>
EigenDecomposition<Real>& EigenDecomposition<Real>::operator= (
    const GMatrix<Real>& mat)
{
    mMatrix = mat;
    mSolver = SymmetricEigensolverGTE<Real>(mat.GetNumRows(), MAX_ITERATIONS);
    mEigenvalues.resize(mat.GetNumRows());
    return *this;
}
//----------------------------------------------------------------------------
template <typename Real>
void EigenDecomposition<Real>::Solve (bool increasingSort)
{
    mSolver.Solve(mMatrix.GetElements(), (increasingSort ? +1 : -1));
    mSolver.GetEigenvalues(&mEigenvalues[0]);
    mSolver.GetEigenvectors(mMatrix.GetElements());
    if (!mSolver.IsRotation())
    {
        // The return matrix of eigenvectors is a reflection.  Negate the
        // first column to make it a rotation.
        for (int r = 0; r < mMatrix.GetNumRows(); ++r)
        {
            mMatrix[r][0] = -mMatrix[r][0];
        }
    }
}
//----------------------------------------------------------------------------
template <typename Real>
Real EigenDecomposition<Real>::GetEigenvalue (int i) const
{
    assertion(0 <= i && i < mMatrix.GetNumRows(),
        "Invalid index in GetEigenvalue\n");
    return mEigenvalues[i];
}
//----------------------------------------------------------------------------
template <typename Real>
const Real* EigenDecomposition<Real>::GetEigenvalues () const
{
    return &mEigenvalues[0];
}
//----------------------------------------------------------------------------
template <typename Real>
Vector2<Real> EigenDecomposition<Real>::GetEigenvector2 (int i) const
{
    int size = mMatrix.GetNumRows();
    if (size == 2)
    {
        Vector2<Real> eigenvector;
        for (int row = 0; row < size; ++row)
        {
            eigenvector[row] = mMatrix[row][i];
        }
        return eigenvector;
    }

    assertion(false, "Mismatched dimension in GetEigenvector2\n");
    return Vector2<Real>::ZERO;
}
//----------------------------------------------------------------------------
template <typename Real>
Matrix2<Real> EigenDecomposition<Real>::GetEigenvectors2 () const
{
    assertion(mMatrix.GetNumRows() == 2,
        "Mismatched dimension in GetEigenvectors2\n");

    Matrix2<Real> eigenvectors;
    for (int row = 0; row < 2; ++row)
    {
        for (int column = 0; column < 2; ++column)
        {
            eigenvectors[row][column] = mMatrix[row][column];
        }
    }
    return eigenvectors;
}
//----------------------------------------------------------------------------
template <typename Real>
Vector3<Real> EigenDecomposition<Real>::GetEigenvector3 (int i) const
{
    int size = mMatrix.GetNumRows();
    if (size == 3)
    {
        Vector3<Real> eigenvector;
        for (int row = 0; row < size; ++row)
        {
            eigenvector[row] = mMatrix[row][i];
        }
        return eigenvector;
    }

    assertion(false, "Mismatched dimension in GetEigenvector3\n");
    return Vector3<Real>::ZERO;
}
//----------------------------------------------------------------------------
template <typename Real>
Matrix3<Real> EigenDecomposition<Real>::GetEigenvectors3 () const
{
    assertion(mMatrix.GetNumRows() == 3,
        "Mismatched dimension in GetEigenvectors2\n");

    Matrix3<Real> eigenvectors;
    for (int row = 0; row < 3; ++row)
    {
        for (int column = 0; column < 3; ++column)
        {
            eigenvectors[row][column] = mMatrix[row][column];
        }
    }
    return eigenvectors;
}
//----------------------------------------------------------------------------
template <typename Real>
GVector<Real> EigenDecomposition<Real>::GetEigenvector (int i) const
{
    return mMatrix.GetColumn(i);
}
//----------------------------------------------------------------------------
template <typename Real>
const GMatrix<Real>& EigenDecomposition<Real>::GetEigenvectors () const
{
    return mMatrix;
}
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Explicit instantiation.
//----------------------------------------------------------------------------
template WM5_MATHEMATICS_ITEM
class EigenDecomposition<float>;

template WM5_MATHEMATICS_ITEM
class EigenDecomposition<double>;
}
//----------------------------------------------------------------------------
