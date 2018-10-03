// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <LowLevel/GteLogger.h>
#include <vector>

// Expose this define to test for out-of-range indices in the classes
// Image2, Image3, and Image4.
//#define GTE_IMAGICS_ASSERT_ON_INVALID_INDEX

namespace gte
{

template <typename PixelType>
class Image
{
public:
    // Construction and destruction.  The variadic argument of the last
    // constructor must contain 'numDimensions' positive values of type
    // 'int', each defining a dimension of the image.
    virtual ~Image();
    Image();
    Image(std::vector<int> const& dimensions);

    // Support for copy semantics.
    Image(Image const& image);
    Image& operator=(Image const& image);

    // Support for move semantics.
    Image(Image&& image);
    Image& operator=(Image&& image);

    // Support for changing the image dimensions.  All pixel data is lost by
    // this operation.
    void Reconstruct(std::vector<int> const& dimensions);

    // Access to image data.
    inline std::vector<int> const& GetDimensions() const;
    inline int GetNumDimensions() const;
    inline int GetDimension(int d) const;
    inline std::vector<size_t> const& GetOffsets() const;
    inline size_t GetOffset(int d) const;
    inline std::vector<PixelType> const& GetPixels() const;
    inline std::vector<PixelType>& GetPixels();
    inline size_t GetNumPixels() const;

    // Conversions between n-dim and 1-dim structures.  The 'coord' arrays
    // must have GetNumDimensions() elements.
    size_t GetIndex(int const* coord) const;
    void GetCoordinates(size_t index, int* coord) const;

    // Access the data as a 1-dimensional array.  The operator[] functions
    // test for valid i when iterator checking is enabled and assert on
    // invalid i.  The Get() functions test for valid i and clamp when
    // invalid; these functions cannot fail.
    inline PixelType& operator[] (size_t i);
    inline PixelType const& operator[] (size_t i) const;
    PixelType& Get(size_t i);
    PixelType const& Get(size_t i) const;

protected:
    std::vector<int> mDimensions;
    std::vector<size_t> mOffsets;
    std::vector<PixelType> mPixels;
};


template <typename PixelType>
Image<PixelType>::~Image()
{
}

template <typename PixelType>
Image<PixelType>::Image()
{
}

template <typename PixelType>
Image<PixelType>::Image(std::vector<int> const& dimensions)
{
    Reconstruct(dimensions);
}

template <typename PixelType>
Image<PixelType>::Image(Image const& image)
{
    *this = image;
}

template <typename PixelType>
Image<PixelType>& Image<PixelType>::operator=(Image const& image)
{
    mDimensions = image.mDimensions;
    mOffsets = image.mOffsets;
    mPixels = image.mPixels;
    return *this;
}

template <typename PixelType>
Image<PixelType>::Image(Image&& image)
{
    *this = std::move(image);
}

template <typename PixelType>
Image<PixelType>& Image<PixelType>::operator=(Image&& image)
{
    mDimensions = std::move(image.mDimensions);
    mOffsets = std::move(image.mOffsets);
    mPixels = std::move(image.mPixels);
    return *this;
}

template <typename PixelType>
void Image<PixelType>::Reconstruct(std::vector<int> const& dimensions)
{
    mDimensions.clear();
    mOffsets.clear();
    mPixels.clear();

    if (dimensions.size() > 0)
    {
        for (auto dim : dimensions)
        {
            if (dim <= 0)
            {
                return;
            }
        }

        mDimensions = dimensions;
        mOffsets.resize(dimensions.size());

        size_t numPixels = 1;
        for (size_t d = 0; d < dimensions.size(); ++d)
        {
            numPixels *= static_cast<size_t>(mDimensions[d]);
        }

        mOffsets[0] = 1;
        for (size_t d = 1; d < dimensions.size(); ++d)
        {
            mOffsets[d] = static_cast<size_t>(mDimensions[d - 1]) * mOffsets[d - 1];
        }

        mPixels.resize(numPixels);
    }
}

template <typename PixelType> inline
std::vector<int> const& Image<PixelType>::GetDimensions() const
{
    return mDimensions;
}

template <typename PixelType> inline
int Image<PixelType>::GetNumDimensions() const
{
    return static_cast<int>(mDimensions.size());
}

template <typename PixelType> inline
int Image<PixelType>::GetDimension(int d) const
{
    return mDimensions[d];
}

template <typename PixelType> inline
std::vector<size_t> const& Image<PixelType>::GetOffsets() const
{
    return mOffsets;
}

template <typename PixelType> inline
size_t Image<PixelType>::GetOffset(int d) const
{
    return mOffsets[d];
}

template <typename PixelType> inline
std::vector<PixelType> const& Image<PixelType>::GetPixels() const
{
    return mPixels;
}

template <typename PixelType> inline
size_t Image<PixelType>::GetNumPixels() const
{
    return mPixels.size();
}

template <typename PixelType> inline
std::vector<PixelType>& Image<PixelType>::GetPixels()
{
    return mPixels;
}

template <typename PixelType>
size_t Image<PixelType>::GetIndex(int const* coord) const
{
    // assert:  coord is array of mNumDimensions elements
    int const numDimensions = static_cast<int>(mDimensions.size());
    size_t index = coord[0];
    for (int d = 1; d < numDimensions; ++d)
    {
        index += mOffsets[d] * coord[d];
    }
    return index;
}

template <typename PixelType>
void Image<PixelType>::GetCoordinates(size_t index, int* coord) const
{
    // assert:  coord is array of numDimensions elements
    int const numDimensions = static_cast<int>(mDimensions.size());
    for (int d = 0; d < numDimensions; ++d)
    {
        coord[d] = index % mDimensions[d];
        index /= mDimensions[d];
    }
}

template <typename PixelType> inline
PixelType& Image<PixelType>::operator[] (size_t i)
{
    return mPixels[i];
}

template <typename PixelType> inline
PixelType const& Image<PixelType>::operator[] (size_t i) const
{
    return mPixels[i];
}

template <typename PixelType>
PixelType& Image<PixelType>::Get(size_t i)
{
    return (i < mPixels.size() ? mPixels[i] : mPixels.front());
}

template <typename PixelType>
PixelType const& Image<PixelType>::Get(size_t i) const
{
    return (i < mPixels.size() ? mPixels[i] : mPixels.front());
}

}
