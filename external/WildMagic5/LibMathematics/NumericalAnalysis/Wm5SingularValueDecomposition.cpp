// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.3 (2014/07/02)

#include "Wm5MathematicsPCH.h"
#include "Wm5SingularValueDecomposition.h"
#include "Wm5SingularValueDecompositionGTE.h"

namespace Wm5
{
//----------------------------------------------------------------------------
template <typename Real>
SingularValueDecomposition<Real>::SingularValueDecomposition (
    const GMatrix<Real>& M, GMatrix<Real>& L, GMatrix<Real>& D,
    GMatrix<Real>& RTranspose)
{
    int numRows = M.GetNumRows();
    int numCols = M.GetNumColumns();
    assertion(numRows >= numCols, "Invalid inputs\n");

    unsigned int maxIterations = 4096;
    SingularValueDecompositionGTE<Real> svd(numRows, numCols, maxIterations);
    svd.Solve(M.GetElements(), -1);  // Sort from largest to smallest.

    // SetSize zeros out the elements.
    L.SetSize(numRows, numRows);
    D.SetSize(numRows, numCols);
    GMatrix<Real> R(numCols, numCols);

    std::vector<Real> singularValues(numCols);
    svd.GetSingularValues(&singularValues[0]);
    for (int i = 0; i < numCols; ++i)
    {
        D(i, i) = singularValues[i];
    }

    svd.GetU(L.GetElements());
    svd.GetV(R.GetElements());
    RTranspose = R.Transpose();
}
//----------------------------------------------------------------------------
template <typename Real>
SingularValueDecomposition<Real>::~SingularValueDecomposition ()
{
}
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Explicit instantiation.
//----------------------------------------------------------------------------
template WM5_MATHEMATICS_ITEM
class SingularValueDecomposition<float>;

template WM5_MATHEMATICS_ITEM
class SingularValueDecomposition<double>;
//----------------------------------------------------------------------------
}
