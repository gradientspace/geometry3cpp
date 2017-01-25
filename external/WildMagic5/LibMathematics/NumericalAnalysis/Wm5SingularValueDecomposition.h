// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.2 (2014/07/02)

#ifndef WM5SINGULARVALUEDECOMPOSITION
#define WM5SINGULARVALUEDECOMPOSITION

#include "Wm5MathematicsLIB.h"
#include "Wm5GMatrix.h"

namespace Wm5
{
template <typename Real>
class WM5_MATHEMATICS_ITEM SingularValueDecomposition
{
public:
    // Singular value decomposition, M = L*D*Transpose(R), where L and R are
    // orthogonal and D is a diagonal matrix whose diagonal entries are
    // nonnegative.  Observe that M is m-by-n with m >= n, L is m-by-m, R is
    // n-by-n, and D is m-by-n; that is, M and D are the same size and not
    // necessarily square.
    SingularValueDecomposition (const GMatrix<Real>& M, GMatrix<Real>& L,
        GMatrix<Real>& D, GMatrix<Real>& RTranspose);

    ~SingularValueDecomposition ();

private:
    // This is a simple wrapper to preserve the Wild Magic 5 interface for
    // singular value decomposition.  The actual implementation lives in
    // Wm5SingularValueDecompositionGTE.{h,cpp}.
};

typedef SingularValueDecomposition<float> SingularValueDecompositionf;
typedef SingularValueDecomposition<double> SingularValueDecompositiond;

}

#endif
