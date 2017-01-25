// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <Imagics/GteImage.h>
#include <array>
#include <type_traits>

namespace gte
{

template <typename PixelType>
class Image1 : public Image
{
public:
    // The default constructor creates a null image.  You can resize the image
    // later with an explicit Resize call, an assignment, or by loading the
    // image from disk.
    virtual ~Image1();
    Image1();

    // Copy the input image using the assignment operator.
    Image1(Image1 const& image);

    // The input dimension must be positive; otherwise, a null image is
    // created.
    Image1(int dimension);

    // If the input image is compatible with 'this', a copy of the input
    // image data occurs.  If the input image is not compatible, 'this' is
    // recreated to be a copy of 'image'.
    Image1& operator= (Image1 const& image);

    // Access the data.  The operator[] function tests for valid i in debug
    // configurations and assert on invalid i.  The Get() functions test for
    // valid i and clamp when invalid (debug and release); these functions
    // cannot fail.
    inline PixelType* GetPixels1D();
    inline PixelType& operator[] (size_t i);
    inline PixelType const& operator[] (size_t i) const;
    PixelType& Get(size_t i);
    PixelType const& Get(size_t i) const;

    // Access the data consistent with the functions that occur for images of
    // higher dimension.  The operator() function tests for valid (x) in debug
    // configurations and assert on invalid (x).  The Get() functions test for
    // valid (x) and clamp when invalid (debug and release); these functions
    // cannot fail.
    inline PixelType& operator() (int x);
    inline PixelType const& operator() (int x) const;
    inline PixelType& operator() (std::array<int, 1> const& coord);
    inline PixelType const& operator() (std::array<int, 1> const& coord) const;
    inline PixelType& Get(std::array<int, 1> coord);
    inline PixelType const& Get(std::array<int, 1> coord) const;

    // Resize an image.  All data is lost from the original image.  The
    // function is convenient for taking a default-constructed image and
    // setting its dimension once it is known.  This avoids an irrelevant
    // memory copy that occurs if instead you were to use the statement
    // image = Image1<PixelType>(dimension).  The return value is 'true'
    // whenever the image is resized (reallocations occurred).
    bool Resize(int dimension);

    // Set all pixels to the specified value.
    void SetAllPixels(PixelType const& value);

    // The required dimensions and pixel type are that of the current image
    // object.  The pixel type for the image class is identified by run-time
    // type information.  If you save an image in one application, it is
    // possible that an equivalent pixel structure is used in another
    // application but the RTTI does not match.  To ignore this mismatch,
    // set 'ignorePixelType' to 'true'.
    bool Load(std::string const& filename, bool ignorePixelType = false);

private:
    // A typed pointer to Image::mRawPixels.
    PixelType* mPixels;

    // Uninitialized, used in the Get(int) calls.
    PixelType mInvalidPixel;
};


template <typename PixelType>
Image1<PixelType>::~Image1()
{
}

template <typename PixelType>
Image1<PixelType>::Image1()
    :
    Image(typeid(PixelType).name(), sizeof(PixelType), 1, 0),
    mPixels(nullptr)
{
}

template <typename PixelType>
Image1<PixelType>::Image1(Image1 const& image)
    :
    Image(image),
    mPixels((PixelType*)mRawPixels)
{
}

template <typename PixelType>
Image1<PixelType>::Image1(int dimension)
    :
    Image(typeid(PixelType).name(), sizeof(PixelType), 1, dimension),
    mPixels((PixelType*)mRawPixels)
{
}

template <typename PixelType>
Image1<PixelType>& Image1<PixelType>::operator= (Image1 const& image)
{
    bool compatible = Copy(image);
    if (!compatible)
    {
        mPixels = (PixelType*)mRawPixels;
    }
    return *this;
}

template <typename PixelType>
inline PixelType* Image1<PixelType>::GetPixels1D()
{
    return mPixels;
}

template <typename PixelType>
inline PixelType& Image1<PixelType>::operator[] (size_t i)
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (mPixels && i < mNumPixels)
    {
        return mPixels[i];
    }
    LogError("No pixels or invalid index " + std::to_string(i) + ".");
    return mInvalidPixel;
#else
    return mPixels[i];
#endif
}

template <typename PixelType>
inline PixelType const& Image1<PixelType>::operator[] (size_t i) const
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (mPixels && i < mNumPixels)
    {
        return mPixels[i];
    }
    LogError("No pixels or invalid index " + std::to_string(i) + ".");
    return mInvalidPixel;
#else
    return mPixels[i];
#endif
}

template <typename PixelType>
PixelType& Image1<PixelType>::Get(size_t i)
{
    if (mPixels)
    {
        if (i >= mNumPixels)
        {
            i = mNumPixels - 1;
        }
        return mPixels[i];
    }
    return mInvalidPixel;
}

template <typename PixelType>
PixelType const& Image1<PixelType>::Get(size_t i) const
{
    if (mPixels)
    {
        if (i >= mNumPixels)
        {
            i = mNumPixels - 1;
        }
        return mPixels[i];
    }
    return mInvalidPixel;
}

template <typename PixelType>
inline PixelType& Image1<PixelType>::operator() (int x)
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (mPixels != 0 && 0 <= x && x < mDimensions[0])
    {
        return mPixels[x];
    }
    LogError("No pixels or invalid index " + std::to_string(x) + ".");
    return mInvalidPixel;
#else
    return mPixels[x];
#endif
}

template <typename PixelType>
inline PixelType const& Image1<PixelType>::operator() (int x) const
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (mPixels != 0 && 0 <= x && x < mDimensions[0])
    {
        return mPixels[x];
    }
    LogError("No pixels or invalid index " + std::to_string(x) + ".");
    return mInvalidPixel;
#else
    return mPixels[x];
#endif
}

template <typename PixelType>
inline PixelType& Image1<PixelType>::operator() (
    std::array<int, 1> const& coord)
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (mPixels != 0 && 0 <= coord[0] && coord[0] < mDimensions[0])
    {
        return mPixels[coord[0]];
    }
    LogError("No pixels or invalid index " + std::to_string(coord[0]) + ".");
    return mInvalidPixel;
#else
    return mPixels[coord[0]];
#endif
}

template <typename PixelType>
inline PixelType const& Image1<PixelType>::operator() (
    std::array<int, 1> const& coord) const
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (mPixels != 0 && 0 <= coord[0] && coord[0] < mDimensions[0])
    {
        return mPixels[coord[0]];
    }
    LogError("No pixels or invalid index " + std::to_string(coord[0]) + ".");
    return mInvalidPixel;
#else
    return mPixels[coord[0]];
#endif
}

template <typename PixelType>
PixelType& Image1<PixelType>::Get(std::array<int, 1> coord)
{
    if (mPixels)
    {
        // Clamp to valid (x).
        if (coord[0] < 0)
        {
            coord[0] = 0;
        }
        else if (coord[0] >= mDimensions[0])
        {
            coord[0] = mDimensions[0] - 1;
        }

        return mPixels[coord[0]];
    }
    return mInvalidPixel;
}

template <typename PixelType>
PixelType const& Image1<PixelType>::Get(std::array<int, 1> coord) const
{
    if (mPixels)
    {
        // Clamp to valid (x).
        if (coord[0] < 0)
        {
            coord[0] = 0;
        }
        else if (coord[0] >= mDimensions[0])
        {
            coord[0] = mDimensions[0] - 1;
        }

        return mPixels[coord[0]];
    }
    return mInvalidPixel;
}

template <typename PixelType>
bool Image1<PixelType>::Resize(int dimension)
{
    if (Image::Resize(1, dimension))
    {
        mPixels = (PixelType*)mRawPixels;
        return true;
    }
    return false;
}

template <typename PixelType>
void Image1<PixelType>::SetAllPixels(PixelType const& value)
{
    if (mPixels)
    {
        for (size_t i = 0; i < mNumPixels; ++i)
        {
            mPixels[i] = value;
        }
    }
}

template <typename PixelType>
bool Image1<PixelType>::Load(std::string const& filename,
    bool ignorePixelType)
{
    std::vector<int> numDimensions(1);
    numDimensions[0] = mNumDimensions;

    std::vector<std::string> pixelTypes(1);
    pixelTypes[0] = mPixelType;
    std::vector<std::string> const* ptr =
        (ignorePixelType ? nullptr : &pixelTypes);

    if (Image::Load(filename, &numDimensions, ptr))
    {
        mPixels = (PixelType*)mRawPixels;
        return true;
    }

    return false;
}


}
