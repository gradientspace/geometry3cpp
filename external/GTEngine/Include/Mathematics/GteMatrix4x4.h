// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <Mathematics/GteMatrix.h>
#include <Mathematics/GteVector4.h>

namespace gte
{

// Template alias for convenience.
template <typename Real>
using Matrix4x4 = Matrix<4, 4, Real>;

// Geometric operations.
template <typename Real>
Matrix4x4<Real> Inverse(Matrix4x4<Real> const& M,
    bool* reportInvertibility = nullptr);

template <typename Real>
Matrix4x4<Real> Adjoint(Matrix4x4<Real> const& M);

template <typename Real>
Real Determinant(Matrix4x4<Real> const& M);

template <typename Real>
Real Trace(Matrix4x4<Real> const& M);


template <typename Real>
Matrix4x4<Real> Inverse(Matrix4x4<Real> const& M, bool* reportInvertibility)
{
    Matrix4x4<Real> inverse;
    bool invertible;
    Real a0 = M(0, 0)*M(1, 1) - M(0, 1)*M(1, 0);
    Real a1 = M(0, 0)*M(1, 2) - M(0, 2)*M(1, 0);
    Real a2 = M(0, 0)*M(1, 3) - M(0, 3)*M(1, 0);
    Real a3 = M(0, 1)*M(1, 2) - M(0, 2)*M(1, 1);
    Real a4 = M(0, 1)*M(1, 3) - M(0, 3)*M(1, 1);
    Real a5 = M(0, 2)*M(1, 3) - M(0, 3)*M(1, 2);
    Real b0 = M(2, 0)*M(3, 1) - M(2, 1)*M(3, 0);
    Real b1 = M(2, 0)*M(3, 2) - M(2, 2)*M(3, 0);
    Real b2 = M(2, 0)*M(3, 3) - M(2, 3)*M(3, 0);
    Real b3 = M(2, 1)*M(3, 2) - M(2, 2)*M(3, 1);
    Real b4 = M(2, 1)*M(3, 3) - M(2, 3)*M(3, 1);
    Real b5 = M(2, 2)*M(3, 3) - M(2, 3)*M(3, 2);
    Real det = a0*b5 - a1*b4 + a2*b3 + a3*b2 - a4*b1 + a5*b0;
    if (det != (Real)0)
    {
        Real invDet = ((Real)1) / det;
        inverse = Matrix4x4<Real>
        {
            (+M(1, 1)*b5 - M(1, 2)*b4 + M(1, 3)*b3)*invDet,
                (-M(0, 1)*b5 + M(0, 2)*b4 - M(0, 3)*b3)*invDet,
                (+M(3, 1)*a5 - M(3, 2)*a4 + M(3, 3)*a3)*invDet,
                (-M(2, 1)*a5 + M(2, 2)*a4 - M(2, 3)*a3)*invDet,
                (-M(1, 0)*b5 + M(1, 2)*b2 - M(1, 3)*b1)*invDet,
                (+M(0, 0)*b5 - M(0, 2)*b2 + M(0, 3)*b1)*invDet,
                (-M(3, 0)*a5 + M(3, 2)*a2 - M(3, 3)*a1)*invDet,
                (+M(2, 0)*a5 - M(2, 2)*a2 + M(2, 3)*a1)*invDet,
                (+M(1, 0)*b4 - M(1, 1)*b2 + M(1, 3)*b0)*invDet,
                (-M(0, 0)*b4 + M(0, 1)*b2 - M(0, 3)*b0)*invDet,
                (+M(3, 0)*a4 - M(3, 1)*a2 + M(3, 3)*a0)*invDet,
                (-M(2, 0)*a4 + M(2, 1)*a2 - M(2, 3)*a0)*invDet,
                (-M(1, 0)*b3 + M(1, 1)*b1 - M(1, 2)*b0)*invDet,
                (+M(0, 0)*b3 - M(0, 1)*b1 + M(0, 2)*b0)*invDet,
                (-M(3, 0)*a3 + M(3, 1)*a1 - M(3, 2)*a0)*invDet,
                (+M(2, 0)*a3 - M(2, 1)*a1 + M(2, 2)*a0)*invDet
        };
        invertible = true;
    }
    else
    {
        inverse.MakeZero();
        invertible = false;
    }

    if (reportInvertibility)
    {
        *reportInvertibility = invertible;
    }
    return inverse;
}

template <typename Real>
Matrix4x4<Real> Adjoint(Matrix4x4<Real> const& M)
{
    Real a0 = M(0, 0)*M(1, 1) - M(0, 1)*M(1, 0);
    Real a1 = M(0, 0)*M(1, 2) - M(0, 2)*M(1, 0);
    Real a2 = M(0, 0)*M(1, 3) - M(0, 3)*M(1, 0);
    Real a3 = M(0, 1)*M(1, 2) - M(0, 2)*M(1, 1);
    Real a4 = M(0, 1)*M(1, 3) - M(0, 3)*M(1, 1);
    Real a5 = M(0, 2)*M(1, 3) - M(0, 3)*M(1, 2);
    Real b0 = M(2, 0)*M(3, 1) - M(2, 1)*M(3, 0);
    Real b1 = M(2, 0)*M(3, 2) - M(2, 2)*M(3, 0);
    Real b2 = M(2, 0)*M(3, 3) - M(2, 3)*M(3, 0);
    Real b3 = M(2, 1)*M(3, 2) - M(2, 2)*M(3, 1);
    Real b4 = M(2, 1)*M(3, 3) - M(2, 3)*M(3, 1);
    Real b5 = M(2, 2)*M(3, 3) - M(2, 3)*M(3, 2);

    return Matrix4x4<Real>
    {
        +M(1, 1)*b5 - M(1, 2)*b4 + M(1, 3)*b3,
            -M(0, 1)*b5 + M(0, 2)*b4 - M(0, 3)*b3,
            +M(3, 1)*a5 - M(3, 2)*a4 + M(3, 3)*a3,
            -M(2, 1)*a5 + M(2, 2)*a4 - M(2, 3)*a3,
            -M(1, 0)*b5 + M(1, 2)*b2 - M(1, 3)*b1,
            +M(0, 0)*b5 - M(0, 2)*b2 + M(0, 3)*b1,
            -M(3, 0)*a5 + M(3, 2)*a2 - M(3, 3)*a1,
            +M(2, 0)*a5 - M(2, 2)*a2 + M(2, 3)*a1,
            +M(1, 0)*b4 - M(1, 1)*b2 + M(1, 3)*b0,
            -M(0, 0)*b4 + M(0, 1)*b2 - M(0, 3)*b0,
            +M(3, 0)*a4 - M(3, 1)*a2 + M(3, 3)*a0,
            -M(2, 0)*a4 + M(2, 1)*a2 - M(2, 3)*a0,
            -M(1, 0)*b3 + M(1, 1)*b1 - M(1, 2)*b0,
            +M(0, 0)*b3 - M(0, 1)*b1 + M(0, 2)*b0,
            -M(3, 0)*a3 + M(3, 1)*a1 - M(3, 2)*a0,
            +M(2, 0)*a3 - M(2, 1)*a1 + M(2, 2)*a0
    };
}

template <typename Real>
Real Determinant(Matrix4x4<Real> const& M)
{
    Real a0 = M(0, 0)*M(1, 1) - M(0, 1)*M(1, 0);
    Real a1 = M(0, 0)*M(1, 2) - M(0, 2)*M(1, 0);
    Real a2 = M(0, 0)*M(1, 3) - M(0, 3)*M(1, 0);
    Real a3 = M(0, 1)*M(1, 2) - M(0, 2)*M(1, 1);
    Real a4 = M(0, 1)*M(1, 3) - M(0, 3)*M(1, 1);
    Real a5 = M(0, 2)*M(1, 3) - M(0, 3)*M(1, 2);
    Real b0 = M(2, 0)*M(3, 1) - M(2, 1)*M(3, 0);
    Real b1 = M(2, 0)*M(3, 2) - M(2, 2)*M(3, 0);
    Real b2 = M(2, 0)*M(3, 3) - M(2, 3)*M(3, 0);
    Real b3 = M(2, 1)*M(3, 2) - M(2, 2)*M(3, 1);
    Real b4 = M(2, 1)*M(3, 3) - M(2, 3)*M(3, 1);
    Real b5 = M(2, 2)*M(3, 3) - M(2, 3)*M(3, 2);
    Real det = a0*b5 - a1*b4 + a2*b3 + a3*b2 - a4*b1 + a5*b0;
    return det;
}

template <typename Real>
Real Trace(Matrix4x4<Real> const& M)
{
    Real trace = M(0, 0) + M(1, 1) + M(2, 2) + M(3, 3);
    return trace;
}


}
