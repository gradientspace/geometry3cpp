// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.12.1 (2015/11/21)

// NOTE: This code was written for the upcoming Geometric Tools Engine but
// has been back-ported to Wild Magic 5 to replace its badly implemented
// version.

#include "Wm5MathematicsPCH.h"
#include "Wm5SingularValueDecompositionGTE.h"

namespace Wm5
{
//----------------------------------------------------------------------------
template <typename Real>
SingularValueDecompositionGTE<Real>::SingularValueDecompositionGTE(
    int numRows, int numCols, unsigned int maxIterations)
    :
    mNumRows(0),
    mNumCols(0),
    mMaxIterations(0)
{
    if (numCols > 1 && numRows >= numCols && maxIterations > 0)
    {
        mNumRows = numRows;
        mNumCols = numCols;
        mMaxIterations = maxIterations;
        mMatrix.resize(numRows * numCols);
        mDiagonal.resize(numCols);
        mSuperdiagonal.resize(numCols - 1);
        mRGivens.reserve(maxIterations*(numCols - 1));
        mLGivens.reserve(maxIterations*(numCols - 1));
        mFixupDiagonal.resize(numCols);
        mPermutation.resize(numCols);
        mVisited.resize(numCols);
        mTwoInvUTU.resize(numCols);
        mTwoInvVTV.resize(numCols - 2);
        mUVector.resize(numRows);
        mVVector.resize(numCols);
        mWVector.resize(numRows);
    }
}
//----------------------------------------------------------------------------
template <typename Real>
unsigned int SingularValueDecompositionGTE<Real>::Solve(Real const* input,
    int sortType)
{
    if (mNumRows > 0)
    {
        int numElements = mNumRows * mNumCols;
        std::copy(input, input + numElements, mMatrix.begin());
        Bidiagonalize();

        // Compute 'threshold = multiplier*epsilon*|B|' as the threshold for
        // diagonal entries effectively zero; that is, |d| <= |threshold|
        // implies that d is (effectively) zero.  TODO: Allow the caller to
        // pass 'multiplier' to the constructor.
        //
        // We will use the L2-norm |B|, which is the length of the elements
        // of B treated as an NM-tuple.  The following code avoids overflow
        // when accumulating the squares of the elements when those elements
        // are large.
        Real maxAbsComp = std::abs(input[0]);
        for (int i = 1; i < numElements; ++i)
        {
            Real absComp = std::abs(input[i]);
            if (absComp > maxAbsComp)
            {
                maxAbsComp = absComp;
            }
        }

        Real norm = (Real)0;
        if (maxAbsComp > (Real)0)
        {
            Real invMaxAbsComp = ((Real)1) / maxAbsComp;
            for (int i = 0; i < numElements; ++i)
            {
                Real ratio = input[i] * invMaxAbsComp;
                norm += ratio * ratio;
            }
            norm = maxAbsComp*sqrt(norm);
        }

        Real const multiplier = (Real)8;  // TODO: Expose to caller.
        Real const epsilon = std::numeric_limits<Real>::epsilon();
        Real const threshold = multiplier * epsilon * norm;

        mRGivens.clear();
        mLGivens.clear();
        for (unsigned int j = 0; j < mMaxIterations; ++j)
        {
            int imin = -1, imax = -1;
            for (int i = mNumCols - 2; i >= 0; --i)
            {
                // When a01 is much smaller than its diagonal neighbors, it is
                // effectively zero.
                Real a00 = mDiagonal[i];
                Real a01 = mSuperdiagonal[i];
                Real a11 = mDiagonal[i + 1];
                Real sum = std::abs(a00) + std::abs(a11);
                if (sum + std::abs(a01) != sum)
                {
                    if (imax == -1)
                    {
                        imax = i;
                    }
                    imin = i;
                }
                else
                {
                    // The superdiagonal term is effectively zero compared to
                    // the neighboring diagonal terms.
                    if (imin >= 0)
                    {
                        break;
                    }
                }
            }

            if (imax == -1)
            {
                // The algorithm has converged.
                EnsureNonnegativeDiagonal();
                ComputePermutation(sortType);
                return j;
            }

            // We need to test diagonal entries of B for zero.  For each zero
            // diagonal entry, zero the superdiagonal.
            if (DiagonalEntriesNonzero(imin, imax, threshold))
            {
                // Process the lower-right-most unreduced bidiagonal block.
                DoGolubKahanStep(imin, imax);
            }
        }
        return 0xFFFFFFFF;
    }
    else
    {
        return 0;
    }
}
//----------------------------------------------------------------------------
template <typename Real>
void SingularValueDecompositionGTE<Real>::GetSingularValues(
    Real* singularValues) const
{
    if (singularValues && mNumCols > 0)
    {
        if (mPermutation[0] >= 0)
        {
            // Sorting was requested.
            for (int i = 0; i < mNumCols; ++i)
            {
                int p = mPermutation[i];
                singularValues[i] = mDiagonal[p];
            }
        }
        else
        {
            // Sorting was not requested.
            memcpy(singularValues, &mDiagonal[0], mNumCols*sizeof(Real));
        }
    }
}
//----------------------------------------------------------------------------
template <typename Real>
void SingularValueDecompositionGTE<Real>::GetU(Real* uMatrix) const
{
    if (!uMatrix || mNumCols == 0)
    {
        // Invalid input or the constructor failed.
        return;
    }

    // Start with the identity matrix.
    std::fill(uMatrix, uMatrix + mNumRows*mNumRows, (Real)0);
    for (int d = 0; d < mNumRows; ++d)
    {
        uMatrix[d + mNumRows*d] = (Real)1;
    }

    // Multiply the Householder reflections using backward accumulation.
    int r, c;
    for (int i0 = mNumCols - 1, i1 = i0 + 1; i0 >= 0; --i0, --i1)
    {
        // Copy the u vector and 2/Dot(u,u) from the matrix.
        Real twoinvudu = mTwoInvUTU[i0];
        Real const* column = &mMatrix[i0];
        mUVector[i0] = (Real)1;
        for (r = i1; r < mNumRows; ++r)
        {
            mUVector[r] = column[mNumCols*r];
        }

        // Compute the w vector.
        mWVector[i0] = twoinvudu;
        for (r = i1; r < mNumRows; ++r)
        {
            mWVector[r] = (Real)0;
            for (c = i1; c < mNumRows; ++c)
            {
                mWVector[r] += mUVector[c] * uMatrix[r + mNumRows*c];
            }
            mWVector[r] *= twoinvudu;
        }

        // Update the matrix, U <- U - u*w^T.
        for (r = i0; r < mNumRows; ++r)
        {
            for (c = i0; c < mNumRows; ++c)
            {
                uMatrix[c + mNumRows*r] -= mUVector[r] * mWVector[c];
            }
        }
    }

    // Multiply the Givens rotations.
    typename std::vector<GivensRotation>::const_iterator givens
        = mLGivens.begin();
    for (/**/; givens != mLGivens.end(); ++givens)
    {
        int j0 = givens->index0;
        int j1 = givens->index1;
        for (r = 0; r < mNumRows; ++r, j0 += mNumRows, j1 += mNumRows)
        {
            Real& q0 = uMatrix[j0];
            Real& q1 = uMatrix[j1];
            Real prd0 = givens->cs * q0 - givens->sn * q1;
            Real prd1 = givens->sn * q0 + givens->cs * q1;
            q0 = prd0;
            q1 = prd1;
        }
    }

    if (mPermutation[0] >= 0)
    {
        // Sorting was requested.
        std::fill(mVisited.begin(), mVisited.end(), 0);
        for (c = 0; c < mNumCols; ++c)
        {
            if (mVisited[c] == 0 && mPermutation[c] != c)
            {
                // The item starts a cycle with 2 or more elements.
                int start = c, current = c, next;
                for (r = 0; r < mNumRows; ++r)
                {
                    mWVector[r] = uMatrix[c + mNumRows*r];
                }
                while ((next = mPermutation[current]) != start)
                {
                    mVisited[current] = 1;
                    for (r = 0; r < mNumRows; ++r)
                    {
                        uMatrix[current + mNumRows*r] =
                            uMatrix[next + mNumRows*r];
                    }
                    current = next;
                }
                mVisited[current] = 1;
                for (r = 0; r < mNumRows; ++r)
                {
                    uMatrix[current + mNumRows*r] = mWVector[r];
                }
            }
        }
    }
}
//----------------------------------------------------------------------------
template <typename Real>
void SingularValueDecompositionGTE<Real>::GetV(Real* vMatrix) const
{
    if (!vMatrix || mNumCols == 0)
    {
        // Invalid input or the constructor failed.
        return;
    }

    // Start with the identity matrix.
    std::fill(vMatrix, vMatrix + mNumCols*mNumCols, (Real)0);
    for (int d = 0; d < mNumCols; ++d)
    {
        vMatrix[d + mNumCols*d] = (Real)1;
    }

    // Multiply the Householder reflections using backward accumulation.
    int i0 = mNumCols - 3;
    int i1 = i0 + 1;
    int i2 = i0 + 2;
    int r, c;
    for (/**/; i0 >= 0; --i0, --i1, --i2)
    {
        // Copy the v vector and 2/Dot(v,v) from the matrix.
        Real twoinvvdv = mTwoInvVTV[i0];
        Real const* row = &mMatrix[mNumCols*i0];
        mVVector[i1] = (Real)1;
        for (r = i2; r < mNumCols; ++r)
        {
            mVVector[r] = row[r];
        }

        // Compute the w vector.
        mWVector[i1] = twoinvvdv;
        for (r = i2; r < mNumCols; ++r)
        {
            mWVector[r] = (Real)0;
            for (c = i2; c < mNumCols; ++c)
            {
                mWVector[r] += mVVector[c] * vMatrix[r + mNumCols*c];
            }
            mWVector[r] *= twoinvvdv;
        }

        // Update the matrix, V <- V - v*w^T.
        for (r = i1; r < mNumCols; ++r)
        {
            for (c = i1; c < mNumCols; ++c)
            {
                vMatrix[c + mNumCols*r] -= mVVector[r] * mWVector[c];
            }
        }
    }

    // Multiply the Givens rotations.
    typename std::vector<GivensRotation>::const_iterator givens
        = mRGivens.begin();
    for (/**/; givens != mRGivens.end(); ++givens)
    {
        int j0 = givens->index0;
        int j1 = givens->index1;
        for (c = 0; c < mNumCols; ++c, j0 += mNumCols, j1 += mNumCols)
        {
            Real& q0 = vMatrix[j0];
            Real& q1 = vMatrix[j1];
            Real prd0 = givens->cs * q0 - givens->sn * q1;
            Real prd1 = givens->sn * q0 + givens->cs * q1;
            q0 = prd0;
            q1 = prd1;
        }
    }

    // Fix-up the diagonal.
    for (r = 0; r < mNumCols; ++r)
    {
        for (c = 0; c < mNumCols; ++c)
        {
            vMatrix[c + mNumCols*r] *= mFixupDiagonal[c];
        }
    }

    if (mPermutation[0] >= 0)
    {
        // Sorting was requested.
        std::fill(mVisited.begin(), mVisited.end(), 0);
        for (c = 0; c < mNumCols; ++c)
        {
            if (mVisited[c] == 0 && mPermutation[c] != c)
            {
                // The item starts a cycle with 2 or more elements.
                int start = c, current = c, next;
                for (r = 0; r < mNumCols; ++r)
                {
                    mWVector[r] = vMatrix[c + mNumCols*r];
                }
                while ((next = mPermutation[current]) != start)
                {
                    mVisited[current] = 1;
                    for (r = 0; r < mNumCols; ++r)
                    {
                        vMatrix[current + mNumCols*r] =
                            vMatrix[next + mNumCols*r];
                    }
                    current = next;
                }
                mVisited[current] = 1;
                for (r = 0; r < mNumCols; ++r)
                {
                    vMatrix[current + mNumCols*r] = mWVector[r];
                }
            }
        }
    }
}
//----------------------------------------------------------------------------
template <typename Real>
void SingularValueDecompositionGTE<Real>::Bidiagonalize()
{
    int r, c;
    for (int i = 0, ip1 = 1; i < mNumCols; ++i, ++ip1)
    {
        // Compute the U-Householder vector.
        Real length = (Real)0;
        for (r = i; r < mNumRows; ++r)
        {
            Real ur = mMatrix[i + mNumCols*r];
            mUVector[r] = ur;
            length += ur * ur;
        }
        Real udu = (Real)1;
        length = sqrt(length);
        if (length >(Real)0)
        {
            Real& u1 = mUVector[i];
            Real sgn = (u1 >= (Real)0 ? (Real)1 : (Real)-1);
            Real invDenom = ((Real)1) / (u1 + sgn * length);
            u1 = (Real)1;
            for (r = ip1; r < mNumRows; ++r)
            {
                Real& ur = mUVector[r];
                ur *= invDenom;
                udu += ur * ur;
            }
        }

        // Compute the rank-1 offset u*w^T.
        Real invudu = (Real)1 / udu;
        Real twoinvudu = invudu * (Real)2;
        for (c = i; c < mNumCols; ++c)
        {
            mWVector[c] = (Real)0;
            for (r = i; r < mNumRows; ++r)
            {
                mWVector[c] += mMatrix[c + mNumCols*r] * mUVector[r];
            }
            mWVector[c] *= twoinvudu;
        }

        // Update the input matrix.
        for (r = i; r < mNumRows; ++r)
        {
            for (c = i; c < mNumCols; ++c)
            {
                mMatrix[c + mNumCols*r] -= mUVector[r] * mWVector[c];
            }
        }

        if (i < mNumCols - 2)
        {
            // Compute the V-Householder vectors.
            length = (Real)0;
            for (c = ip1; c < mNumCols; ++c)
            {
                Real vc = mMatrix[c + mNumCols*i];
                mVVector[c] = vc;
                length += vc * vc;
            }
            Real vdv = (Real)1;
            length = sqrt(length);
            if (length >(Real)0)
            {
                Real& v1 = mVVector[ip1];
                Real sgn = (v1 >= (Real)0 ? (Real)1 : (Real)-1);
                Real invDenom = ((Real)1) / (v1 + sgn * length);
                v1 = (Real)1;
                for (c = ip1 + 1; c < mNumCols; ++c)
                {
                    Real& vc = mVVector[c];
                    vc *= invDenom;
                    vdv += vc * vc;
                }
            }

            // Compute the rank-1 offset w*v^T.
            Real invvdv = (Real)1 / vdv;
            Real twoinvvdv = invvdv * (Real)2;
            for (r = i; r < mNumRows; ++r)
            {
                mWVector[r] = (Real)0;
                for (c = ip1; c < mNumCols; ++c)
                {
                    mWVector[r] += mMatrix[c + mNumCols*r] * mVVector[c];
                }
                mWVector[r] *= twoinvvdv;
            }

            // Update the input matrix.
            for (r = i; r < mNumRows; ++r)
            {
                for (c = ip1; c < mNumCols; ++c)
                {
                    mMatrix[c + mNumCols*r] -= mWVector[r] * mVVector[c];
                }
            }

            mTwoInvVTV[i] = twoinvvdv;
            for (c = i + 2; c < mNumCols; ++c)
            {
                mMatrix[c + mNumCols*i] = mVVector[c];
            }
        }

        mTwoInvUTU[i] = twoinvudu;
        for (r = ip1; r < mNumRows; ++r)
        {
            mMatrix[i + mNumCols*r] = mUVector[r];
        }
    }

    // Copy the diagonal and subdiagonal for cache coherence in the
    // Golub-Kahan iterations.
    int k, ksup = mNumCols - 1, index = 0, delta = mNumCols + 1;
    for (k = 0; k < ksup; ++k, index += delta)
    {
        mDiagonal[k] = mMatrix[index];
        mSuperdiagonal[k] = mMatrix[index + 1];
    }
    mDiagonal[k] = mMatrix[index];
}
//----------------------------------------------------------------------------
template <typename Real>
void SingularValueDecompositionGTE<Real>::GetSinCos(Real x, Real y, Real& cs,
    Real& sn)
{
    // Solves sn*x + cs*y = 0 robustly.
    Real tau;
    if (y != (Real)0)
    {
        if (std::abs(y) > std::abs(x))
        {
            tau = -x / y;
            sn = ((Real)1) / sqrt(((Real)1) + tau*tau);
            cs = sn * tau;
        }
        else
        {
            tau = -y / x;
            cs = ((Real)1) / sqrt(((Real)1) + tau*tau);
            sn = cs * tau;
        }
    }
    else
    {
        cs = (Real)1;
        sn = (Real)0;
    }
}
//----------------------------------------------------------------------------
template <typename Real>
bool SingularValueDecompositionGTE<Real>::DiagonalEntriesNonzero(int imin,
    int imax, Real threshold)
{
    for (int i = imin; i <= imax; ++i)
    {
        if (std::abs(mDiagonal[i]) <= threshold)
        {
            // Use planar rotations to case the superdiagonal entry out of
            // the matrix, thus producing a row of zeros.
            Real x, z, cs, sn;
            Real y = mSuperdiagonal[i];
            mSuperdiagonal[i] = (Real)0;
            for (int j = i + 1; j <= imax + 1; ++j)
            {
                x = mDiagonal[j];
                GetSinCos(x, y, cs, sn);
                mLGivens.push_back(GivensRotation(i, j, cs, sn));
                mDiagonal[j] = cs*x - sn*y;
                if (j <= imax)
                {
                    z = mSuperdiagonal[j];
                    mSuperdiagonal[j] = cs*z;
                    y = sn*z;
                }
            }
            return false;
        }
    }
    return true;
}
//----------------------------------------------------------------------------
template <typename Real>
void SingularValueDecompositionGTE<Real>::DoGolubKahanStep(int imin, int imax)
{
    // The implicit shift.  Compute the eigenvalue u of the lower-right 2x2
    // block of A = B^T*B that is closer to b11.
    Real f0 = (imax >= (Real)1 ? mSuperdiagonal[imax - 1] : (Real)0);
    Real d1 = mDiagonal[imax];
    Real f1 = mSuperdiagonal[imax];
    Real d2 = mDiagonal[imax + 1];
    Real a00 = d1*d1 + f0*f0;
    Real a01 = d1*f1;
    Real a11 = d2*d2 + f1*f1;
    Real dif = (a00 - a11) * (Real)0.5;
    Real sgn = (dif >= (Real)0 ? (Real)1 : (Real)-1);
    Real a01sqr = a01 * a01;
    Real u = a11 - a01sqr / (dif + sgn*sqrt(dif*dif + a01sqr));
    Real x = mDiagonal[imin] * mDiagonal[imin] - u;
    Real y = mDiagonal[imin] * mSuperdiagonal[imin];

    Real a12, a21, a22, a23, cs, sn;
    Real a02 = (Real)0;
    int i0 = imin - 1, i1 = imin, i2 = imin + 1;
    for (/**/; i1 <= imax; ++i0, ++i1, ++i2)
    {
        // Compute the Givens rotation G and save it for use in computing
        // V in U^T*A*V = S.
        GetSinCos(x, y, cs, sn);
        mRGivens.push_back(GivensRotation(i1, i2, cs, sn));

        // Update B0 = B*G.
        if (i1 > imin)
        {
            mSuperdiagonal[i0] = cs*mSuperdiagonal[i0] - sn*a02;
        }

        a11 = mDiagonal[i1];
        a12 = mSuperdiagonal[i1];
        a22 = mDiagonal[i2];
        mDiagonal[i1] = cs*a11 - sn*a12;
        mSuperdiagonal[i1] = sn*a11 + cs*a12;
        mDiagonal[i2] = cs*a22;
        a21 = -sn*a22;

        // Update the parameters for the next Givens rotations.
        x = mDiagonal[i1];
        y = a21;

        // Compute the Givens rotation G and save it for use in computing
        // U in U^T*A*V = S.
        GetSinCos(x, y, cs, sn);
        mLGivens.push_back(GivensRotation(i1, i2, cs, sn));

        // Update B1 = G^T*B0.
        a11 = mDiagonal[i1];
        a12 = mSuperdiagonal[i1];
        a22 = mDiagonal[i2];
        mDiagonal[i1] = cs*a11 - sn*a21;
        mSuperdiagonal[i1] = cs*a12 - sn*a22;
        mDiagonal[i2] = sn*a12 + cs*a22;

        if (i1 < imax)
        {
            a23 = mSuperdiagonal[i2];
            a02 = -sn*a23;
            mSuperdiagonal[i2] = cs*a23;

            // Update the parameters for the next Givens rotations.
            x = mSuperdiagonal[i1];
            y = a02;
        }
    }
}
//----------------------------------------------------------------------------
template <typename Real>
void SingularValueDecompositionGTE<Real>::EnsureNonnegativeDiagonal()
{
    for (int i = 0; i < mNumCols; ++i)
    {
        if (mDiagonal[i] >= (Real)0)
        {
            mFixupDiagonal[i] = (Real)1;
        }
        else
        {
            mDiagonal[i] = -mDiagonal[i];
            mFixupDiagonal[i] = (Real)-1;
        }
    }
}
//----------------------------------------------------------------------------
template <typename Real>
void SingularValueDecompositionGTE<Real>::ComputePermutation(int sortType)
{
    if (sortType == 0)
    {
        // Set a flag for GetSingularValues() and GetOrthogonalMatrices() to
        // know that sorted output was not requested.
        mPermutation[0] = -1;
        return;
    }

    // Compute the permutation induced by sorting.  Initially, we start with
    // the identity permutation I = (0,1,...,N-1).
    std::vector<SortItem> items(mNumCols);
    int i;
    for (i = 0; i < mNumCols; ++i)
    {
        items[i].singularValue = mDiagonal[i];
        items[i].index = i;
    }

    if (sortType > 0)
    {
        std::sort(items.begin(), items.end(), std::less<SortItem>());
    }
    else
    {
        std::sort(items.begin(), items.end(), std::greater<SortItem>());
    }

    typename std::vector<SortItem>::const_iterator item = items.begin();
    for (i = 0; item != items.end(); ++item, ++i)
    {
        mPermutation[i] = item->index;
    }

    // GetOrthogonalMatrices() has nontrivial code for computing the
    // orthogonal U and V from the reflections and rotations.  To avoid
    // complicating the code further when sorting is requested, U and V are
    // computed as in the unsorted case.  We then need to swap columns of
    // U and V to be consistent with the sorting of the singular values.  To
    // minimize copying due to column swaps, we use permutation P.  The
    // minimum number of transpositions to obtain P from I is N minus the
    // number of cycles of P.  Each cycle is reordered with a minimum number
    // of transpositions; that is, the singular items are cyclically swapped,
    // leading to a minimum amount of copying.  For example, if there is a
    // cycle i0 -> i1 -> i2 -> i3, then the copying is
    //   save = singularitem[i0];
    //   singularitem[i1] = singularitem[i2];
    //   singularitem[i2] = singularitem[i3];
    //   singularitem[i3] = save;
}
//----------------------------------------------------------------------------
template <typename Real>
SingularValueDecompositionGTE<Real>::GivensRotation::GivensRotation()
{
    // No default initialization for fast creation of std::vector of objects
    // of this type.
}
//----------------------------------------------------------------------------
template <typename Real>
SingularValueDecompositionGTE<Real>::GivensRotation::GivensRotation(
    int inIndex0, int inIndex1, Real inCs, Real inSn)
    :
    index0(inIndex0),
    index1(inIndex1),
    cs(inCs),
    sn(inSn)
{
}
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Explicit instantiation.
//----------------------------------------------------------------------------
template WM5_MATHEMATICS_ITEM
class SingularValueDecompositionGTE<float>;

template WM5_MATHEMATICS_ITEM
class SingularValueDecompositionGTE<double>;
//----------------------------------------------------------------------------
}
