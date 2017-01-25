// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.12.1 (2014/07/09)

// NOTE: This code was written for the upcoming Geometric Tools Engine but
// has been back-ported to Wild Magic 5 because it has better quality than
// its previous version.

#include "Wm5MathematicsPCH.h"
#include "Wm5SymmetricEigensolverGTE.h"

namespace Wm5
{
//----------------------------------------------------------------------------
template <typename Real>
SymmetricEigensolverGTE<Real>::SymmetricEigensolverGTE(int size,
    unsigned int maxIterations)
    :
    mSize(0),
    mMaxIterations(0),
    mIsRotation(-1)
{
    if (size > 1 && maxIterations > 0)
    {
        mSize = size;
        mMaxIterations = maxIterations;
        mMatrix.resize(size*size);
        mDiagonal.resize(size);
        mSuperdiagonal.resize(size - 1);
        mGivens.reserve(maxIterations * (size - 1));
        mPermutation.resize(size);
        mVisited.resize(size);
        mPVector.resize(size);
        mVVector.resize(size);
        mWVector.resize(size);
    }
}
//----------------------------------------------------------------------------
template <typename Real>
unsigned int SymmetricEigensolverGTE<Real>::Solve(Real const* input,
    int sortType)
{
    if (mSize > 0)
    {
        std::copy(input, input + mSize*mSize, mMatrix.begin());
        Tridiagonalize();

        mGivens.clear();
        for (unsigned int j = 0; j < mMaxIterations; ++j)
        {
            int imin = -1, imax = -1;
            for (int i = mSize - 2; i >= 0; --i)
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
                ComputePermutation(sortType);
                return j;
            }

            // Process the lower-right-most unreduced tridiagonal block.
            DoQRImplicitShift(imin, imax);
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
void SymmetricEigensolverGTE<Real>::GetEigenvalues(Real* eigenvalues) const
{
    if (eigenvalues && mSize > 0)
    {
        if (mPermutation[0] >= 0)
        {
            // Sorting was requested.
            for (int i = 0; i < mSize; ++i)
            {
                int p = mPermutation[i];
                eigenvalues[i] = mDiagonal[p];
            }
        }
        else
        {
            // Sorting was not requested.
            memcpy(eigenvalues, &mDiagonal[0], mSize*sizeof(Real));
        }
    }
}
//----------------------------------------------------------------------------
template <typename Real>
void SymmetricEigensolverGTE<Real>::GetEigenvectors(Real* eigenvectors) const
{
    if (eigenvectors && mSize > 0)
    {
        // Start with the identity matrix.
        std::fill(eigenvectors, eigenvectors + mSize*mSize, (Real)0);
        for (int d = 0; d < mSize; ++d)
        {
            eigenvectors[d + mSize*d] = (Real)1;
        }

        // Multiply the Householder reflections using backward accumulation.
        int r, c;
        for (int i = mSize - 3, rmin = i + 1; i >= 0; --i, --rmin)
        {
            // Copy the v vector and 2/Dot(v,v) from the matrix.
            Real const* column = &mMatrix[i];
            Real twoinvvdv = column[mSize*(i + 1)];
            for (r = 0; r < i + 1; ++r)
            {
                mVVector[r] = (Real)0;
            }
            mVVector[r] = (Real)1;
            for (++r; r < mSize; ++r)
            {
                mVVector[r] = column[mSize*r];
            }

            // Compute the w vector.
            for (r = 0; r < mSize; ++r)
            {
                mWVector[r] = (Real)0;
                for (c = rmin; c < mSize; ++c)
                {
                    mWVector[r] += mVVector[c] * eigenvectors[r + mSize*c];
                }
                mWVector[r] *= twoinvvdv;
            }

            // Update the matrix, Q <- Q - v*w^T.
            for (r = rmin; r < mSize; ++r)
            {
                for (c = 0; c < mSize; ++c)
                {
                    eigenvectors[c + mSize*r] -= mVVector[r] * mWVector[c];
                }
            }
        }

        // Multiply the Givens rotations.
		typename std::vector<GivensRotation>::const_iterator givens
            = mGivens.begin();
        for (/**/; givens != mGivens.end(); ++givens)
        {
            for (r = 0; r < mSize; ++r)
            {
                int j = givens->index + mSize*r;
                Real& q0 = eigenvectors[j];
                Real& q1 = eigenvectors[j + 1];
                Real prd0 = givens->cs * q0 - givens->sn * q1;
                Real prd1 = givens->sn * q0 + givens->cs * q1;
                q0 = prd0;
                q1 = prd1;
            }
        }

        mIsRotation = 1 - (mSize & 1);
        if (mPermutation[0] >= 0)
        {
            // Sorting was requested.
            std::fill(mVisited.begin(), mVisited.end(), 0);
            for (int i = 0; i < mSize; ++i)
            {
                if (mVisited[i] == 0 && mPermutation[i] != i)
                {
                    // The item starts a cycle with 2 or more elements.
                    mIsRotation = 1 - mIsRotation;
                    int start = i, current = i, j, next;
                    for (j = 0; j < mSize; ++j)
                    {
                        mPVector[j] = eigenvectors[i + mSize*j];
                    }
                    while ((next = mPermutation[current]) != start)
                    {
                        mVisited[current] = 1;
                        for (j = 0; j < mSize; ++j)
                        {
                            eigenvectors[current + mSize*j] =
                                eigenvectors[next + mSize*j];
                        }
                        current = next;
                    }
                    mVisited[current] = 1;
                    for (j = 0; j < mSize; ++j)
                    {
                        eigenvectors[current + mSize*j] = mPVector[j];
                    }
                }
            }
        }
    }
}
//----------------------------------------------------------------------------
template <typename Real>
bool SymmetricEigensolverGTE<Real>::IsRotation() const
{
    if (mSize > 0)
    {
        if (mIsRotation == -1)
        {
            // Without sorting, the matrix is a rotation when size is even.
            mIsRotation = 1 - (mSize & 1);
            if (mPermutation[0] >= 0)
            {
                // With sorting, the matrix is a rotation when the number of
                // cycles in the permutation is even.
                std::fill(mVisited.begin(), mVisited.end(), 0);
                for (int i = 0; i < mSize; ++i)
                {
                    if (mVisited[i] == 0 && mPermutation[i] != i)
                    {
                        // The item starts a cycle with 2 or more elements.
                        int start = i, current = i, next;
                        while ((next = mPermutation[current]) != start)
                        {
                            mVisited[current] = 1;
                            current = next;
                        }
                        mVisited[current] = 1;
                    }
                }
            }
        }
        return mIsRotation == 1;
    }
    else
    {
        return false;
    }
}
//----------------------------------------------------------------------------
template <typename Real>
void SymmetricEigensolverGTE<Real>::GetEigenvector(int c, Real* eigenvector)
    const
{
    if (0 <= c && c < mSize)
    {
        // y = H*x, then x and y are swapped for the next H
        Real* x = eigenvector;
        Real* y = &mPVector[0];

        // Start with the Euclidean basis vector.
        memset(x, 0, mSize*sizeof(Real));
        if (mPermutation[c] >= 0)
        {
            x[mPermutation[c]] = (Real)1;
        }
        else
        {
            x[c] = (Real)1;
        }

        // Apply the Givens rotations.
        typename std::vector<GivensRotation>::const_reverse_iterator givens =
            mGivens.rbegin();
        for (/**/; givens != mGivens.rend(); ++givens)
        {
            Real& xr = x[givens->index];
            Real& xrp1 = x[givens->index + 1];
            Real tmp0 = givens->cs * xr + givens->sn * xrp1;
            Real tmp1 = -givens->sn * xr + givens->cs * xrp1;
            xr = tmp0;
            xrp1 = tmp1;
        }

        // Apply the Householder reflections.
        for (int i = mSize - 3; i >= 0; --i)
        {
            // Get the Householder vector v.
            Real const* column = &mMatrix[i];
            Real twoinvvdv = column[mSize*(i + 1)];
            int r;
            for (r = 0; r < i + 1; ++r)
            {
                y[r] = x[r];
            }

            // Compute s = Dot(x,v) * 2/v^T*v.
            Real s = x[r];  // r = i+1, v[i+1] = 1
            for (int j = r + 1; j < mSize; ++j)
            {
                s += x[j] * column[mSize*j];
            }
            s *= twoinvvdv;

            y[r] = x[r] - s;  // v[i+1] = 1

            // Compute the remaining components of y.
            for (++r; r < mSize; ++r)
            {
                y[r] = x[r] - s * column[mSize*r];
            }

            std::swap(x, y);
        }
        // The final product is stored in x.

        if (x != eigenvector)
        {
            memcpy(eigenvector, x, mSize * sizeof(Real));
        }
    }
}
//----------------------------------------------------------------------------
template <typename Real>
void SymmetricEigensolverGTE<Real>::Tridiagonalize()
{
    int r, c;
    for (int i = 0, ip1 = 1; i < mSize - 2; ++i, ++ip1)
    {
        // Compute the Householder vector.  Read the initial vector from the
        // row of the matrix.
        Real length = (Real)0;
        for (r = 0; r < ip1; ++r)
        {
            mVVector[r] = (Real)0;
        }
        for (r = ip1; r < mSize; ++r)
        {
            Real vr = mMatrix[r + mSize*i];
            mVVector[r] = vr;
            length += vr * vr;
        }
        Real vdv = (Real)1;
        length = sqrt(length);
        if (length >(Real)0)
        {
            Real& v1 = mVVector[ip1];
            Real sgn = (v1 >= (Real)0 ? (Real)1 : (Real)-1);
            Real invDenom = ((Real)1) / (v1 + sgn * length);
            v1 = (Real)1;
            for (r = ip1 + 1; r < mSize; ++r)
            {
                Real& vr = mVVector[r];
                vr *= invDenom;
                vdv += vr * vr;
            }
        }

        // Compute the rank-1 offsets v*w^T and w*v^T.
        Real invvdv = (Real)1 / vdv;
        Real twoinvvdv = invvdv * (Real)2;
        Real pdvtvdv = (Real)0;
        for (r = i; r < mSize; ++r)
        {
            mPVector[r] = (Real)0;
            for (c = i; c < r; ++c)
            {
                mPVector[r] += mMatrix[r + mSize*c] * mVVector[c];
            }
            for (/**/; c < mSize; ++c)
            {
                mPVector[r] += mMatrix[c + mSize*r] * mVVector[c];
            }
            mPVector[r] *= twoinvvdv;
            pdvtvdv += mPVector[r] * mVVector[r];
        }

        pdvtvdv *= invvdv;
        for (r = i; r < mSize; ++r)
        {
            mWVector[r] = mPVector[r] - pdvtvdv * mVVector[r];
        }

        // Update the input matrix.
        for (r = i; r < mSize; ++r)
        {
            Real vr = mVVector[r];
            Real wr = mWVector[r];
            Real offset = vr * wr * (Real)2;
            mMatrix[r + mSize*r] -= offset;
            for (c = r + 1; c < mSize; ++c)
            {
                offset = vr * mWVector[c] + wr * mVVector[c];
                mMatrix[c + mSize*r] -= offset;
            }
        }

        // Copy the vector to column i of the matrix.  The 0-valued components
        // at indices 0 through i are not stored.  The 1-valued component at
        // index i+1 is also not stored; instead, the quantity 2/Dot(v,v) is
        // stored for use in eigenvector construction. That construction must
        // take into account the implied components that are not stored.
        mMatrix[i + mSize*ip1] = twoinvvdv;
        for (r = ip1 + 1; r < mSize; ++r)
        {
            mMatrix[i + mSize*r] = mVVector[r];
        }
    }

    // Copy the diagonal and subdiagonal entries for cache coherence in
    // the QR iterations.
    int k, ksup = mSize - 1, index = 0, delta = mSize + 1;
    for (k = 0; k < ksup; ++k, index += delta)
    {
        mDiagonal[k] = mMatrix[index];
        mSuperdiagonal[k] = mMatrix[index + 1];
    }
    mDiagonal[k] = mMatrix[index];
}
//----------------------------------------------------------------------------
template <typename Real>
void SymmetricEigensolverGTE<Real>::GetSinCos(Real x, Real y, Real& cs,
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
void SymmetricEigensolverGTE<Real>::DoQRImplicitShift(int imin, int imax)
{
    // The implicit shift.  Compute the eigenvalue u of the lower-right 2x2
    // block that is closer to a11.
    Real a00 = mDiagonal[imax];
    Real a01 = mSuperdiagonal[imax];
    Real a11 = mDiagonal[imax + 1];
    Real dif = (a00 - a11) * (Real)0.5;
    Real sgn = (dif >= (Real)0 ? (Real)1 : (Real)-1);
    Real a01sqr = a01 * a01;
    Real u = a11 - a01sqr / (dif + sgn*sqrt(dif*dif + a01sqr));
    Real x = mDiagonal[imin] - u;
    Real y = mSuperdiagonal[imin];

    Real a12, a22, a23, tmp11, tmp12, tmp21, tmp22, cs, sn;
    Real a02 = (Real)0;
    int i0 = imin - 1, i1 = imin, i2 = imin + 1;
    for (/**/; i1 <= imax; ++i0, ++i1, ++i2)
    {
        // Compute the Givens rotation and save it for use in computing the
        // eigenvectors.
        GetSinCos(x, y, cs, sn);
        mGivens.push_back(GivensRotation(i1, cs, sn));

        // Update the tridiagonal matrix.  This amounts to updating a 4x4
        // subblock,
        //   b00 b01 b02 b03
        //   b01 b11 b12 b13
        //   b02 b12 b22 b23
        //   b03 b13 b23 b33
        // The four corners (b00, b03, b33) do not change values.  The
        // The interior block {{b11,b12},{b12,b22}} is updated on each pass.
        // For the first pass, the b0c values are out of range, so only
        // the values (b13, b23) change.  For the last pass, the br3 values
        // are out of range, so only the values (b01, b02) change.  For
        // passes between first and last, the values (b01, b02, b13, b23)
        // change.
        if (i1 > imin)
        {
            mSuperdiagonal[i0] = cs*mSuperdiagonal[i0] - sn*a02;
        }

        a11 = mDiagonal[i1];
        a12 = mSuperdiagonal[i1];
        a22 = mDiagonal[i2];
        tmp11 = cs*a11 - sn*a12;
        tmp12 = cs*a12 - sn*a22;
        tmp21 = sn*a11 + cs*a12;
        tmp22 = sn*a12 + cs*a22;
        mDiagonal[i1] = cs*tmp11 - sn*tmp12;
        mSuperdiagonal[i1] = sn*tmp11 + cs*tmp12;
        mDiagonal[i2] = sn*tmp21 + cs*tmp22;

        if (i1 < imax)
        {
            a23 = mSuperdiagonal[i2];
            a02 = -sn*a23;
            mSuperdiagonal[i2] = cs*a23;

            // Update the parameters for the next Givens rotation.
            x = mSuperdiagonal[i1];
            y = a02;
        }
    }
}
//----------------------------------------------------------------------------
template <typename Real>
void SymmetricEigensolverGTE<Real>::ComputePermutation(int sortType)
{
    mIsRotation = -1;

    if (sortType == 0)
    {
        // Set a flag for GetEigenvalues() and GetEigenvectors() to know
        // that sorted output was not requested.
        mPermutation[0] = -1;
        return;
    }

    // Compute the permutation induced by sorting.  Initially, we start with
    // the identity permutation I = (0,1,...,N-1).
    std::vector<SortItem> items(mSize);
    int i;
    for (i = 0; i < mSize; ++i)
    {
        items[i].eigenvalue = mDiagonal[i];
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

    // GetEigenvectors() has nontrivial code for computing the orthogonal Q
    // from the reflections and rotations.  To avoid complicating the code
    // further when sorting is requested, Q is computed as in the unsorted
    // case.  We then need to swap columns of Q to be consistent with the
    // sorting of the eigenvalues.  To minimize copying due to column swaps,
    // we use permutation P.  The minimum number of transpositions to obtain
    // P from I is N minus the number of cycles of P.  Each cycle is reordered
    // with a minimum number of transpositions; that is, the eigenitems are
    // cyclically swapped, leading to a minimum amount of copying.  For
    // example, if there is a cycle i0 -> i1 -> i2 -> i3, then the copying is
    //   save = eigenitem[i0];
    //   eigenitem[i1] = eigenitem[i2];
    //   eigenitem[i2] = eigenitem[i3];
    //   eigenitem[i3] = save;
}
//----------------------------------------------------------------------------
template <typename Real>
SymmetricEigensolverGTE<Real>::GivensRotation::GivensRotation()
{
    // No default initialization for fast creation of std::vector of objects
    // of this type.
}
//----------------------------------------------------------------------------
template <typename Real>
SymmetricEigensolverGTE<Real>::GivensRotation::GivensRotation(int inIndex,
    Real inCs, Real inSn)
    :
    index(inIndex),
    cs(inCs),
    sn(inSn)
{
}
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Explicit instantiation.
//----------------------------------------------------------------------------
template WM5_MATHEMATICS_ITEM
class SymmetricEigensolverGTE<float>;

template WM5_MATHEMATICS_ITEM
class SymmetricEigensolverGTE<double>;
//----------------------------------------------------------------------------
}
