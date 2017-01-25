// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <GTEngineDEF.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

// This is an implementation of the iterative algorithm described in
// http://www.geometrictools.com/Documentation/RobustEigenSymmetric3x3.pdf
// The code does not use GTEngine objects.

namespace gte
{

template <typename Real>
class SymmetricEigensolver3x3
{
public:
    // The input matrix must be symmetric, so only the unique elements must
    // be specified: a00, a01, a02, a11, a12, and a22.
    //
    // If 'aggressive' is 'true', the iterations occur until a superdiagonal
    // entry is exactly zero.  If 'aggressive' is 'false', the iterations
    // occur until a superdiagonal entry is effectively zero compared to the
    // sum of magnitudes of its diagonal neighbors.  Generally, the
    // nonaggressive convergence is acceptable.
    //
    // The order of the eigenvalues is specified by sortType: -1 (decreasing),
    // 0 (no sorting), or +1 (increasing).  When sorted, the eigenvectors are
    // ordered accordingly, and {evec[0], evec[1], evec[2]} is guaranteed to
    // be a right-handed orthonormal set.  The return value is the number of
    // iterations used by the algorithm.

    int operator()(Real a00, Real a01, Real a02, Real a11, Real a12, Real a22,
        bool aggressive, int sortType, std::array<Real, 3>& eval,
        std::array<std::array<Real, 3>, 3>& evec) const;

private:
    // Update Q = Q*G in-place using G = {{c,0,-s},{s,0,c},{0,0,1}}.
    void Update0(Real Q[3][3], Real c, Real s) const;

    // Update Q = Q*G in-place using G = {{0,1,0},{c,0,s},{-s,0,c}}.
    void Update1(Real Q[3][3], Real c, Real s) const;

    // Update Q = Q*H in-place using H = {{c,s,0},{s,-c,0},{0,0,1}}.
    void Update2(Real Q[3][3], Real c, Real s) const;

    // Update Q = Q*H in-place using H = {{1,0,0},{0,c,s},{0,s,-c}}.
    void Update3(Real Q[3][3], Real c, Real s) const;

    // Normalize (u,v) robustly, avoiding floating-point overflow in the sqrt
    // call.  The normalized pair is (cs,sn) with cs <= 0.  If (u,v) = (0,0),
    // the function returns (cs,sn) = (-1,0).  When used to generate a
    // Householder reflection, it does not matter whether (cs,sn) or (-cs,-sn)
    // is used.  When generating a Givens reflection, cs = cos(2*theta) and
    // sn = sin(2*theta).  Having a negative cosine for the double-angle
    // term ensures that the single-angle terms c = cos(theta) and
    // s = sin(theta) satisfy |c| <= |s|.
    void GetCosSin(Real u, Real v, Real& cs, Real& sn) const;

    // The convergence test.  When 'aggressive' is 'true', the superdiagonal
    // test is "bSuper == 0".  When 'aggressive' is 'false', the superdiagonal
    // test is "|bDiag0| + |bDiag1| + |bSuper| == |bDiag0| + |bDiag1|, which
    // means bSuper is effectively zero compared to the sizes of the diagonal
    // entries.
    bool Converged(bool aggressive, Real bDiag0, Real bDiag1,
        Real bSuper) const;

    // Support for sorting the eigenvalues and eigenvectors.  The output
    // (i0,i1,i2) is a permutation of (0,1,2) so that d[i0] <= d[i1] <= d[i2].
    // The 'bool' return indicates whether the permutation is odd.  If it is
    // not, the handedness of the Q matrix must be adjusted.
    bool Sort(std::array<Real, 3> const& d, int& i0, int& i1, int& i2) const;
};


template <typename Real>
int SymmetricEigensolver3x3<Real>::operator()(Real a00, Real a01, Real a02,
    Real a11, Real a12, Real a22, bool aggressive, int sortType,
    std::array<Real, 3>& eval, std::array<std::array<Real, 3>, 3>& evec) const
{
    // Compute the Householder reflection H and B = H*A*H, where b02 = 0.
    Real const zero = (Real)0, one = (Real)1, half = (Real)0.5;
    bool isRotation = false;
    Real c, s;
    GetCosSin(a12, -a02, c, s);
    Real Q[3][3] = { { c, s, zero }, { s, -c, zero }, { zero, zero, one } };
    Real term0 = c * a00 + s * a01;
    Real term1 = c * a01 + s * a11;
    Real b00 = c * term0 + s * term1;
    Real b01 = s * term0 - c * term1;
    term0 = s * a00 - c * a01;
    term1 = s * a01 - c * a11;
    Real b11 = s * term0 - c * term1;
    Real b12 = s * a02 - c * a12;
    Real b22 = a22;

    // Givens reflections, B' = G^T*B*G, preserve tridiagonal matrices.
    int const maxIteration = 2 * (1 + std::numeric_limits<Real>::digits -
        std::numeric_limits<Real>::min_exponent);
    int iteration;
    Real c2, s2;

    if (std::abs(b12) <= std::abs(b01))
    {
        Real saveB00, saveB01, saveB11;
        for (iteration = 0; iteration < maxIteration; ++iteration)
        {
            // Compute the Givens reflection.
            GetCosSin(half * (b00 - b11), b01, c2, s2);
            s = sqrt(half * (one - c2));  // >= 1/sqrt(2)
            c = half * s2 / s;

            // Update Q by the Givens reflection.
            Update0(Q, c, s);
            isRotation = !isRotation;

            // Update B <- Q^T*B*Q, ensuring that b02 is zero and |b12| has
            // strictly decreased.
            saveB00 = b00;
            saveB01 = b01;
            saveB11 = b11;
            term0 = c * saveB00 + s * saveB01;
            term1 = c * saveB01 + s * saveB11;
            b00 = c * term0 + s * term1;
            b11 = b22;
            term0 = c * saveB01 - s * saveB00;
            term1 = c * saveB11 - s * saveB01;
            b22 = c * term1 - s * term0;
            b01 = s * b12;
            b12 = c * b12;

            if (Converged(aggressive, b00, b11, b01))
            {
                // Compute the Householder reflection.
                GetCosSin(half * (b00 - b11), b01, c2, s2);
                s = sqrt(half * (one - c2));
                c = half * s2 / s;  // >= 1/sqrt(2)

                // Update Q by the Householder reflection.
                Update2(Q, c, s);
                isRotation = !isRotation;

                // Update D = Q^T*B*Q.
                saveB00 = b00;
                saveB01 = b01;
                saveB11 = b11;
                term0 = c * saveB00 + s * saveB01;
                term1 = c * saveB01 + s * saveB11;
                b00 = c * term0 + s * term1;
                term0 = s * saveB00 - c * saveB01;
                term1 = s * saveB01 - c * saveB11;
                b11 = s * term0 - c * term1;
                break;
            }
        }
    }
    else
    {
        Real saveB11, saveB12, saveB22;
        for (iteration = 0; iteration < maxIteration; ++iteration)
        {
            // Compute the Givens reflection.
            GetCosSin(half * (b22 - b11), b12, c2, s2);
            s = sqrt(half * (one - c2));  // >= 1/sqrt(2)
            c = half * s2 / s;

            // Update Q by the Givens reflection.
            Update1(Q, c, s);
            isRotation = !isRotation;

            // Update B <- Q^T*B*Q, ensuring that b02 is zero and |b12| has
            // strictly decreased.  MODIFY...
            saveB11 = b11;
            saveB12 = b12;
            saveB22 = b22;
            term0 = c * saveB22 + s * saveB12;
            term1 = c * saveB12 + s * saveB11;
            b22 = c * term0 + s * term1;
            b11 = b00;
            term0 = c * saveB12 - s * saveB22;
            term1 = c * saveB11 - s * saveB12;
            b00 = c * term1 - s * term0;
            b12 = s * b01;
            b01 = c * b01;

            if (Converged(aggressive, b11, b22, b12))
            {
                // Compute the Householder reflection.
                GetCosSin(half * (b11 - b22), b12, c2, s2);
                s = sqrt(half * (one - c2));
                c = half * s2 / s;  // >= 1/sqrt(2)

                // Update Q by the Householder reflection.
                Update3(Q, c, s);
                isRotation = !isRotation;

                // Update D = Q^T*B*Q.
                saveB11 = b11;
                saveB12 = b12;
                saveB22 = b22;
                term0 = c * saveB11 + s * saveB12;
                term1 = c * saveB12 + s * saveB22;
                b11 = c * term0 + s * term1;
                term0 = s * saveB11 - c * saveB12;
                term1 = s * saveB12 - c * saveB22;
                b22 = s * term0 - c * term1;
                break;
            }
        }
    }

    std::array<Real, 3> diagonal = { b00, b11, b22 };
    int i0, i1, i2;
    if (sortType >= 1)
    {
        // diagonal[i0] <= diagonal[i1] <= diagonal[i2]
        bool isOdd = Sort(diagonal, i0, i1, i2);
        if (!isOdd)
        {
            isRotation = !isRotation;
        }
    }
    else if (sortType <= -1)
    {
        // diagonal[i0] >= diagonal[i1] >= diagonal[i2]
        bool isOdd = Sort(diagonal, i0, i1, i2);
        std::swap(i0, i2);  // (i0,i1,i2)->(i2,i1,i0) is odd
        if (isOdd)
        {
            isRotation = !isRotation;
        }
    }
    else
    {
        i0 = 0;
        i1 = 1;
        i2 = 2;
    }

    eval[0] = diagonal[i0];
    eval[1] = diagonal[i1];
    eval[2] = diagonal[i2];
    evec[0][0] = Q[0][i0];
    evec[0][1] = Q[1][i0];
    evec[0][2] = Q[2][i0];
    evec[1][0] = Q[0][i1];
    evec[1][1] = Q[1][i1];
    evec[1][2] = Q[2][i1];
    evec[2][0] = Q[0][i2];
    evec[2][1] = Q[1][i2];
    evec[2][2] = Q[2][i2];

    // Ensure the columns of Q form a right-handed set.
    if (!isRotation)
    {
        for (int j = 0; j < 3; ++j)
        {
            evec[2][j] = -evec[2][j];
        }
    }

    return iteration;
}

template <typename Real>
void SymmetricEigensolver3x3<Real>::Update0(Real Q[3][3], Real c, Real s)
const
{
    for (int r = 0; r < 3; ++r)
    {
        Real tmp0 = c * Q[r][0] + s * Q[r][1];
        Real tmp1 = Q[r][2];
        Real tmp2 = c * Q[r][1] - s * Q[r][0];
        Q[r][0] = tmp0;
        Q[r][1] = tmp1;
        Q[r][2] = tmp2;
    }
}

template <typename Real>
void SymmetricEigensolver3x3<Real>::Update1(Real Q[3][3], Real c, Real s)
const
{
    for (int r = 0; r < 3; ++r)
    {
        Real tmp0 = c * Q[r][1] - s * Q[r][2];
        Real tmp1 = Q[r][0];
        Real tmp2 = c * Q[r][2] + s * Q[r][1];
        Q[r][0] = tmp0;
        Q[r][1] = tmp1;
        Q[r][2] = tmp2;
    }
}

template <typename Real>
void SymmetricEigensolver3x3<Real>::Update2(Real Q[3][3], Real c, Real s)
const
{
    for (int r = 0; r < 3; ++r)
    {
        Real tmp0 = c * Q[r][0] + s * Q[r][1];
        Real tmp1 = s * Q[r][0] - c * Q[r][1];
        Q[r][0] = tmp0;
        Q[r][1] = tmp1;
    }
}

template <typename Real>
void SymmetricEigensolver3x3<Real>::Update3(Real Q[3][3], Real c, Real s)
const
{
    for (int r = 0; r < 3; ++r)
    {
        Real tmp0 = c * Q[r][1] + s * Q[r][2];
        Real tmp1 = s * Q[r][1] - c * Q[r][2];
        Q[r][1] = tmp0;
        Q[r][2] = tmp1;
    }
}

template <typename Real>
void SymmetricEigensolver3x3<Real>::GetCosSin(Real u, Real v, Real& cs,
    Real& sn) const
{
    Real maxAbsComp = std::max(std::abs(u), std::abs(v));
    if (maxAbsComp > (Real)0)
    {
        u /= maxAbsComp;  // in [-1,1]
        v /= maxAbsComp;  // in [-1,1]
        Real length = sqrt(u * u + v * v);
        cs = u / length;
        sn = v / length;
        if (cs > (Real)0)
        {
            cs = -cs;
            sn = -sn;
        }
    }
    else
    {
        cs = (Real)-1;
        sn = (Real)0;
    }
}

template <typename Real>
bool SymmetricEigensolver3x3<Real>::Converged(bool aggressive, Real bDiag0,
    Real bDiag1, Real bSuper) const
{
    if (aggressive)
    {
        return bSuper == (Real)0;
    }
    else
    {
        Real sum = std::abs(bDiag0) + std::abs(bDiag1);
        return sum + std::abs(bSuper) == sum;
    }
}

template <typename Real>
bool SymmetricEigensolver3x3<Real>::Sort(std::array<Real, 3> const& d,
    int& i0, int& i1, int& i2) const
{
    bool odd;
    if (d[0] < d[1])
    {
        if (d[2] < d[0])
        {
            i0 = 2; i1 = 0; i2 = 1; odd = true;
        }
        else if (d[2] < d[1])
        {
            i0 = 0; i1 = 2; i2 = 1; odd = false;
        }
        else
        {
            i0 = 0; i1 = 1; i2 = 2; odd = true;
        }
    }
    else
    {
        if (d[2] < d[1])
        {
            i0 = 2; i1 = 1; i2 = 0; odd = false;
        }
        else if (d[2] < d[0])
        {
            i0 = 1; i1 = 2; i2 = 0; odd = true;
        }
        else
        {
            i0 = 1; i1 = 0; i2 = 2; odd = false;
        }
    }
    return odd;
}


}
