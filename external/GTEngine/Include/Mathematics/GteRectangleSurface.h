// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <Mathematics/GteDarbouxFrame.h>
#include <memory>

namespace gte
{

template <typename Real>
class RectangleSurface
{
public:
    // Construction.
    RectangleSurface(
        std::shared_ptr<ParametricSurface<3, Real>> const& surface);

    // Support for generating triangle meshes that approximate the rectangle
    // surface.  The input surface must be rectangular (IsRectangular() must
    // return 'true'.)  Use Set(...) first.  The information is persistent for
    // calls to the Get* functions.  Use the GetNumVertices() and
    // GetNumTriangles() functions second so you can allocate the correct
    // number of elements in the arrays passed to the Get* functions.
    void Set(int numUSamples, int numVSamples);

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
    inline Vector3<Real>& Position(int i) const;
    inline Vector3<Real>& Normal(int i) const;
    inline Vector2<Real>& TCoord(int i) const;

    std::shared_ptr<ParametricSurface<3, Real>> mSurface;
    int mNumUSamples, mNumVSamples;
    int mNumVertices, mNumTriangles;
    mutable char* mPositions;
    mutable char* mNormals;
    mutable char* mTCoords;
    mutable size_t mStride;
};


template <typename Real>
RectangleSurface<Real>::RectangleSurface(
    std::shared_ptr<ParametricSurface<3, Real>> const& surface)
    :
    mSurface(surface),
    mNumUSamples(0),
    mNumVSamples(0),
    mNumVertices(0),
    mNumTriangles(0),
    mPositions(nullptr),
    mNormals(nullptr),
    mTCoords(nullptr),
    mStride(0)
{
}

template <typename Real>
void RectangleSurface<Real>::Set(int numUSamples, int numVSamples)
{
    mNumUSamples = numUSamples;
    mNumVSamples = numVSamples;

    mNumVertices = numUSamples * numVSamples;
    mNumTriangles = 2 * (numUSamples - 1) * (numVSamples - 1);
}

template <typename Real> inline
int RectangleSurface<Real>::GetNumVertices() const
{
    return mNumVertices;
}

template <typename Real> inline
int RectangleSurface<Real>::GetNumTriangles() const
{
    return mNumTriangles;
}

template <typename Real>
void RectangleSurface<Real>::GetPositions(void* positions, size_t stride)
const
{
    mPositions = static_cast<char*>(positions);
    mStride = stride;

    Real uMin = mSurface->GetUMin();
    Real uDelta = (mSurface->GetUMax() - uMin) / (Real)(mNumUSamples - 1);
    Real vMin = mSurface->GetVMin();
    Real vDelta = (mSurface->GetVMax() - vMin) / (Real)(mNumVSamples - 1);
    for (int vIndex = 0, i = 0; vIndex < mNumVSamples; ++vIndex)
    {
        Real v = vMin + vDelta * vIndex;
        for (int uIndex = 0; uIndex < mNumUSamples; ++uIndex, ++i)
        {
            Real u = uMin + uDelta * uIndex;
            Position(i) = mSurface->GetPosition(u, v);
        }
    }

    mPositions = nullptr;
    mStride = 0;
}

template <typename Real>
void RectangleSurface<Real>::GetNormals(void* normals, size_t stride) const
{
    mNormals = static_cast<char*>(normals);
    mStride = stride;

    Real uMin = mSurface->GetUMin();
    Real uDelta = (mSurface->GetUMax() - uMin) / (Real)(mNumUSamples - 1);
    Real vMin = mSurface->GetVMin();
    Real vDelta = (mSurface->GetVMax() - vMin) / (Real)(mNumVSamples - 1);
    DarbouxFrame3<Real> darboux(*mSurface.get());
    for (int vIndex = 0, i = 0; vIndex < mNumVSamples; ++vIndex)
    {
        Real v = vMin + vDelta * vIndex;
        for (int uIndex = 0; uIndex < mNumUSamples; ++uIndex, ++i)
        {
            Real u = uMin + uDelta * uIndex;
            Vector3<Real> position, tangent0, tangent1, normal;
            darboux(u, v, position, tangent0, tangent1, normal);
            Normal(i) = normal;
        }
    }

    mNormals = nullptr;
    mStride = 0;
}

template <typename Real>
void RectangleSurface<Real>::GetTCoords(void* tcoords, size_t stride) const
{
    mTCoords = static_cast<char*>(tcoords);
    mStride = stride;

    Real uMin = mSurface->GetUMin();
    Real uDelta = (mSurface->GetUMax() - uMin) / (Real)(mNumUSamples - 1);
    Real vMin = mSurface->GetVMin();
    Real vDelta = (mSurface->GetVMax() - vMin) / (Real)(mNumVSamples - 1);
    Vector2<Real> tcoord;
    for (int vIndex = 0, i = 0; vIndex < mNumVSamples; ++vIndex)
    {
        tcoord[1] = vMin + vDelta * vIndex;
        for (int uIndex = 0; uIndex < mNumUSamples; ++uIndex, ++i)
        {
            tcoord[0] = uMin + uDelta * uIndex;
            TCoord(i) = tcoord;
        }
    }

    mTCoords = nullptr;
    mStride = 0;
}

template <typename Real>
void RectangleSurface<Real>::GetIndices(int* indices) const
{
    for (int v = 0, i = 0; v < mNumVSamples - 1; ++v)
    {
        int i0 = i;
        int i1 = i0 + 1;
        i += mNumUSamples;
        int i2 = i;
        int i3 = i2 + 1;
        for (int u = 0; u < mNumUSamples - 1; ++u, indices += 6)
        {
            indices[0] = i0;
            indices[1] = i2;
            indices[2] = i1;
            indices[3] = i1;
            indices[4] = i2;
            indices[5] = i3;
            i0++;
            i1++;
            i2++;
            i3++;
        }
    }
}

template <typename Real> inline
Vector3<Real>& RectangleSurface<Real>::Position(int i) const
{
    return *reinterpret_cast<Vector3<Real>*>(mPositions + i * mStride);
}

template <typename Real> inline
Vector3<Real>& RectangleSurface<Real>::Normal(int i) const
{
    return *reinterpret_cast<Vector3<Real>*>(mNormals + i * mStride);
}

template <typename Real> inline
Vector2<Real>& RectangleSurface<Real>::TCoord(int i) const
{
    return *reinterpret_cast<Vector2<Real>*>(mTCoords + i * mStride);
}


}
