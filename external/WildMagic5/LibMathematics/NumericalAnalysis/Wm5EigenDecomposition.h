// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.2 (2014/07/02)

#ifndef WM5EIGENDECOMPOSITION_H
#define WM5EIGENDECOMPOSITION_H

#include "Wm5MathematicsLIB.h"
#include "Wm5Matrix2.h"
#include "Wm5Matrix3.h"
#include "Wm5Matrix4.h"
#include "Wm5GMatrix.h"
#include "Wm5SymmetricEigensolverGTE.h"

namespace Wm5
{

template <typename Real>
class WM5_MATHEMATICS_ITEM EigenDecomposition
{
public:
    // Construction and destruction.  The matrix of an eigensystem must be
    // symmetric.
    EigenDecomposition (int size);
    EigenDecomposition (const Matrix2<Real>& mat);
    EigenDecomposition (const Matrix3<Real>& mat);
    EigenDecomposition (const GMatrix<Real>& mat);
    ~EigenDecomposition ();

    // Set the matrix for the eigensystem.
    Real& operator() (int row, int column);
    EigenDecomposition& operator= (const Matrix2<Real>& mat);
    EigenDecomposition& operator= (const Matrix3<Real>& mat);
    EigenDecomposition& operator= (const GMatrix<Real>& mat);

    // Solve the eigensystem.  Set 'increasingSort' to 'true' when you want
    // the eigenvalues to be sorted in increasing order; otherwise, the
    // eigenvalues are sorted in decreasing order.
    void Solve (bool increasingSort);

    // Get the results.  The calls to GetEigenvector2, GetEigenvectors2,
    // GetEigenvector3, and GetEigenvector3 should be made only if you know
    // that the eigensystem is of the corresponding size.
    Real GetEigenvalue (int i) const;
    const Real* GetEigenvalues () const;
    Vector2<Real> GetEigenvector2 (int i) const;
    Matrix2<Real> GetEigenvectors2 () const;
    Vector3<Real> GetEigenvector3 (int i) const;
    Matrix3<Real> GetEigenvectors3 () const;
    GVector<Real> GetEigenvector (int i) const;
    const GMatrix<Real>& GetEigenvectors () const;

private:
    // This is a simple wrapper to preserve the Wild Magic 5 interface for
    // eigendecomposition of symmetric matrices.  The actual implementation
    // lives in Wm5SymmetricEigensolverGTE.{h,cpp}.
    GMatrix<Real> mMatrix;

    // The eigenvectors are stored in mMatrix after the decomposition.
    enum { MAX_ITERATIONS = 4096 };
    SymmetricEigensolverGTE<Real> mSolver;
    std::vector<Real> mEigenvalues;
};

typedef EigenDecomposition<float> EigenDecompositionf;
typedef EigenDecomposition<double> EigenDecompositiond;

}

#endif
