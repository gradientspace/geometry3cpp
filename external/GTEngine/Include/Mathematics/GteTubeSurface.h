// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <Mathematics/GteFrenetFrame.h>
#include <Mathematics/GteConstants.h>
#include <functional>
#include <memory>

namespace gte
{

template <typename Real>
class TubeSurface
{
public:
    // Construction.
    TubeSurface(std::shared_ptr<ParametricCurve<3, Real>> const& medial,
        std::function<Real(Real)> const& radial);

    // Support for generating triangle meshes that approximate the tube
    // surface.  Use Set(...) first.  The information is persistent for calls
    // to the Get* functions.  Use the GetNumVertices() and GetNumTriangles()
    // functions second so you can allocate the correct number of elements in
    // the arrays passed to the Get* functions.
    void Set(int numMedialSamples, int numRadialSamples, bool closed,
        bool sampleByArcLength, bool outsideView);

    // The return values of these functions should be used to allocate the
    // array inputs to the other Get* functions.
    inline int GetNumVertices() const;
    inline int GetNumTriangles() const;

    // The positions must be 3-tuples separated by 'stride' bytes.  The
    // normals must be 3-tuples separated by 'stride' bytes.  The tcoords
    // must be 2-tuples separated by 'stride' bytes.
    void GetPositions(void* positions, size_t stride) const;
    void GetNormals(void* normals, size_t stride) const;
    void GetTCoords(void* tcoords, size_t stride) const;
    void GetIndices(int* indices) const;

private:
    inline int Index(int medialIndex, int radialIndex) const;
    inline Vector3<Real>& Position(int i) const;
    inline Vector3<Real>& Normal(int i) const;
    inline Vector2<Real>& TCoord(int i) const;

    std::shared_ptr<ParametricCurve<3, Real>> mMedial;
    std::function<Real(Real)> mRadial;
    int mNumMedialSamples;
    int mNumRadialSamples;
    bool mClosed, mSampleByArcLength, mOutsideView;
    int mNumVertices;
    int mNumTriangles;
    std::vector<Real> mCosAngle, mSinAngle;
    mutable char* mPositions;
    mutable char* mNormals;
    mutable char* mTCoords;
    mutable size_t mStride;
};


template <typename Real>
TubeSurface<Real>::TubeSurface(
    std::shared_ptr<ParametricCurve<3, Real>> const& medial,
    std::function<Real(Real)> const& radial)
    :
    mMedial(medial),
    mRadial(radial),
    mNumMedialSamples(0),
    mNumRadialSamples(0),
    mClosed(false),
    mSampleByArcLength(false),
    mOutsideView(false),
    mNumVertices(0),
    mNumTriangles(0),
    mPositions(nullptr),
    mNormals(nullptr),
    mTCoords(nullptr),
    mStride(0)
{
}

template <typename Real>
void TubeSurface<Real>::Set(int numMedialSamples, int numRadialSamples,
    bool closed, bool sampleByArcLength, bool outsideView)
{
    mNumMedialSamples = numMedialSamples;
    mNumRadialSamples = numRadialSamples;
    mClosed = closed;
    mSampleByArcLength = sampleByArcLength;
    mOutsideView = outsideView;

    int offset = (closed ? 0 : 1);
    mNumVertices = (numRadialSamples + 1) * (numMedialSamples + offset);
    mNumTriangles = 2 * numRadialSamples * (numMedialSamples - offset);

    mCosAngle.resize(numRadialSamples + 1);
    mSinAngle.resize(numRadialSamples + 1);
    Real invRadialSamples = ((Real)1) / (Real)numRadialSamples;
    for (int i = 0; i < numRadialSamples; ++i)
    {
        float angle = i * invRadialSamples * (Real)GTE_C_TWO_PI;
        mCosAngle[i] = cos(angle);
        mSinAngle[i] = sin(angle);
    }
    mCosAngle[numRadialSamples] = mCosAngle[0];
    mSinAngle[numRadialSamples] = mSinAngle[0];
}

template <typename Real> inline
int TubeSurface<Real>::GetNumVertices() const
{
    return mNumVertices;
}

template <typename Real> inline
int TubeSurface<Real>::GetNumTriangles() const
{
    return mNumTriangles;
}

template <typename Real>
void TubeSurface<Real>::GetPositions(void* positions, size_t stride) const
{
    mPositions = static_cast<char*>(positions);
    mStride = stride;

    int offset = (mClosed ? 0 : 1);
    Real invDenom = ((Real)1) / (Real)(mNumMedialSamples - offset);
    FrenetFrame3<Real> frenet(*mMedial.get());
    Vector3<Real> position, tangent, normal, binormal;
    int save, v, row, col;
    Real factor, t, radius;

    if (mSampleByArcLength)
    {
        factor = mMedial->GetTotalLength() * invDenom;
        for (row = 0, v = 0; row < mNumMedialSamples; ++row, ++v)
        {
            t = mMedial->GetTime(row * factor);
            radius = mRadial(t);
            frenet(t, position, tangent, normal, binormal);
            for (col = 0, save = v; col < mNumRadialSamples; ++col, ++v)
            {
                Position(v) = position + radius * (mCosAngle[col] * normal +
                    mSinAngle[col] * binormal);
            }
            Position(v) = Position(save);
        }
    }
    else
    {
        Real tmin = mMedial->GetTMin();
        factor = (mMedial->GetTMax() - tmin) * invDenom;
        for (row = 0, v = 0; row < mNumMedialSamples; ++row, ++v)
        {
            t = tmin + row * factor;
            radius = mRadial(t);
            frenet(t, position, tangent, normal, binormal);
            for (col = 0, save = v; col < mNumRadialSamples; ++col, ++v)
            {
                Position(v) = position + radius * (mCosAngle[col] * normal +
                    mSinAngle[col] * binormal);
            }
            Position(v) = Position(save);
        }
    }

    if (mClosed)
    {
        for (col = 0; col <= mNumRadialSamples; ++col)
        {
            int i0 = Index(0, col);
            int i1 = Index(mNumMedialSamples, col);
            Position(i1) = Position(i0);
        }
    }

    mPositions = nullptr;
    mStride = 0;
}

template <typename Real>
void TubeSurface<Real>::GetNormals(void* normals, size_t stride) const
{
    mNormals = static_cast<char*>(normals);
    mStride = stride;

    int row, rowM1, rowP1, col, colM1, colP1;
    Vector3<Real> dir0, dir1;

    // Compute the interior normals (central differences).
    for (row = 1; row <= mNumMedialSamples - 2; ++row)
    {
        for (col = 0; col < mNumRadialSamples; ++col)
        {
            colM1 = (col > 0 ? col - 1 : mNumRadialSamples - 1);
            colP1 = col + 1;
            rowM1 = row - 1;
            rowP1 = row + 1;
            dir0 = Position(Index(row, colM1)) - Position(Index(row, colP1));
            dir1 = Position(Index(rowM1, col)) - Position(Index(rowP1, col));
            Normal(Index(row, col)) = UnitCross(dir0, dir1);
        }
        Normal(Index(row, mNumRadialSamples)) = Normal(Index(row, 0));
    }

    // Compute the boundary normals.
    if (mClosed)
    {
        // Compute with central differences.
        for (col = 0; col < mNumRadialSamples; ++col)
        {
            colM1 = (col > 0 ? col - 1 : mNumRadialSamples - 1);
            colP1 = col + 1;

            // row = 0
            dir0 = Position(Index(0, colM1)) - Position(Index(0, colP1));
            dir1 = Position(Index(mNumMedialSamples - 1, col)) -
                Position(Index(1, col));
            Normal(col) = UnitCross(dir0, dir1);

            // row = max
            Normal(Index(mNumMedialSamples, col)) = Normal(Index(0, col));
        }
        Normal(Index(0, mNumRadialSamples)) = Normal(Index(0, 0));
        Normal(Index(mNumMedialSamples, mNumRadialSamples)) =
            Normal(Index(mNumMedialSamples, 0));
    }
    else
    {
        // Compute with one-sided differences.

        // row = 0
        for (col = 0; col < mNumRadialSamples; ++col)
        {
            colM1 = (col > 0 ? col - 1 : mNumRadialSamples - 1);
            colP1 = col + 1;
            dir0 = Position(Index(0, colM1)) - Position(Index(0, colP1));
            dir1 = Position(Index(0, col)) - Position(Index(1, col));
            Normal(Index(0, col)) = UnitCross(dir0, dir1);
        }
        Normal(Index(0, mNumRadialSamples)) = Normal(Index(0, 0));

        // row = max-1
        row = mNumMedialSamples - 1;
        for (col = 0; col < mNumRadialSamples; ++col)
        {
            colM1 = (col > 0 ? col - 1 : mNumRadialSamples - 1);
            colP1 = col + 1;
            dir0 = Position(Index(row, colM1)) - Position(Index(row, colP1));
            dir1 = Position(Index(row - 1, col)) - Position(Index(row, col));
            Normal(col) = UnitCross(dir0, dir1);
        }
        Normal(Index(row, mNumRadialSamples)) = Normal(Index(row, 0));
    }

    mNormals = nullptr;
    mStride = 0;
}

template <typename Real>
void TubeSurface<Real>::GetTCoords(void* tcoords, size_t stride) const
{
    mTCoords = static_cast<char*>(tcoords);
    mStride = stride;

    int maxRow = mNumMedialSamples - (mClosed ? 0 : 1);
    int maxCol = mNumRadialSamples;
    Real invMaxRow = ((Real)1) / (Real)maxRow;
    Real invMaxCol = ((Real)1) / (Real)maxCol;
    Vector2<Real> tcoord;
    for (int row = 0, v = 0; row <= maxRow; ++row)
    {
        tcoord[1] = row * invMaxRow;
        for (int col = 0; col <= maxCol; ++col, ++v)
        {
            tcoord[0] = col * invMaxCol;
            TCoord(Index(row, col)) = tcoord;
        }
    }

    mTCoords = nullptr;
    mStride = 0;
}

template <typename Real>
void TubeSurface<Real>::GetIndices(int* indices) const
{
    int rowStart = 0;
    for (int row = 0; row < mNumMedialSamples - 1; ++row)
    {
        int i0 = rowStart;
        int i1 = i0 + 1;
        rowStart += mNumRadialSamples + 1;
        int i2 = rowStart;
        int i3 = i2 + 1;
        if (mOutsideView)
        {
            for (int i = 0; i < mNumRadialSamples; ++i, indices += 6)
            {
                indices[0] = i0++;
                indices[1] = i1;
                indices[2] = i2;
                indices[3] = i1++;
                indices[4] = i3++;
                indices[5] = i2++;
            }
        }
        else
        {
            for (int i = 0; i < mNumRadialSamples; ++i, indices += 6)
            {
                indices[0] = i0++;
                indices[1] = i2;
                indices[2] = i1;
                indices[3] = i1++;
                indices[4] = i2++;
                indices[5] = i3++;
            }
        }
    }

    if (mClosed)
    {
        int i0 = rowStart;
        int i1 = i0 + 1;
        int i2 = 0;
        int i3 = i2 + 1;
        if (mOutsideView)
        {
            for (int i = 0; i < mNumRadialSamples; ++i, indices += 6)
            {
                indices[0] = i0++;
                indices[1] = i1;
                indices[2] = i2;
                indices[3] = i1++;
                indices[4] = i3++;
                indices[5] = i2++;
            }
        }
        else
        {
            for (int i = 0; i < mNumRadialSamples; ++i, indices += 6)
            {
                indices[0] = i0++;
                indices[1] = i2;
                indices[2] = i1;
                indices[3] = i1++;
                indices[4] = i2++;
                indices[5] = i3++;
            }
        }
    }
}

template <typename Real> inline
int TubeSurface<Real>::Index(int medialIndex, int radialIndex) const
{
    return radialIndex + (mNumRadialSamples + 1) * medialIndex;
}

template <typename Real> inline
Vector3<Real>& TubeSurface<Real>::Position(int i) const
{
    return *reinterpret_cast<Vector3<Real>*>(mPositions + i * mStride);
}

template <typename Real> inline
Vector3<Real>& TubeSurface<Real>::Normal(int i) const
{
    return *reinterpret_cast<Vector3<Real>*>(mNormals + i * mStride);
}

template <typename Real> inline
Vector2<Real>& TubeSurface<Real>::TCoord(int i) const
{
    return *reinterpret_cast<Vector2<Real>*>(mTCoords + i * mStride);
}


}
