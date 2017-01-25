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
class Image2 : public Image
{
public:
    // The default constructor creates a null image.  You can resize the image
    // later with an explicit Resize call, an assignment, or by loading the
    // image from disk.
    virtual ~Image2();
    Image2();

    // Copy the input image using the assignment operator.
    Image2(Image2 const& image);

    // The input dimensions must be positive; otherwise, a null image is
    // created.
    Image2(int dimension0, int dimension1);

    // If the input image is compatible with 'this', a copy of the input
    // image data occurs.  If the input image is not compatible, 'this' is
    // recreated to be a copy of 'image'.
    Image2& operator= (Image2 const& image);

    // The input array must have the correct number of pixels as determined by
    // the image parameters.  Use at your own risk, because we cannot verify
    // the compatibility.
    virtual void SetRawPixels(char* rawPixels);

    // Conversion between 1-dimensional indices and 2-dimensional coordinates.
    inline size_t GetIndex(int x, int y) const;
    inline size_t GetIndex(std::array<int, 2> const& coord) const;
    inline void GetCoordinates(size_t index, int& x, int& y) const;
    inline std::array<int, 2> GetCoordinates(size_t index) const;

    // Access the data as a 1-dimensional array.  The operator[] functions
    // test for valid i in debug configurations and assert on invalid i.  The
    // Get() functions test for valid i and clamp when invalid (debug and
    // release); these functions cannot fail.
    inline PixelType* GetPixels1D() const;
    inline PixelType& operator[] (size_t i);
    inline PixelType const& operator[] (size_t i) const;
    PixelType& Get(size_t i);
    PixelType const& Get(size_t i) const;

    // Access the data as a 2-dimensional array.  Pixel (x,y) is accessed
    // as "pixels2D[y][x]".  The operator() functions test for valid (x,y) in
    // debug configurations and assert on invalid (x,y).  The Get() functions
    // test for valid (x,y) and clamp when invalid (debug and release); these
    // functions cannot fail.
    inline PixelType** GetPixels2D() const;
    inline PixelType& operator() (int x, int y);
    inline PixelType const& operator() (int x, int y) const;
    inline PixelType& operator() (std::array<int, 2> const& coord);
    inline PixelType const& operator() (std::array<int,2> const& coord) const;
    inline PixelType& Get(int x, int y);
    inline PixelType const& Get(int x, int y) const;
    inline PixelType& Get(std::array<int, 2> coord);
    inline PixelType const& Get(std::array<int, 2> coord) const;

    // In the following discussion, u and v are in {-1,1}.  Given a pixel
    // (x,y), the 4-connected neighbors have relative offsets (u,0) and
    // (0,v).  The 8-connected neighbors include the 4-connected neighbors
    // and have additional relative offsets (u,v).  The corner neighbors
    // have relative offsets (0,0), (1,0), (0,1), and (1,1) in that order.
    // The full neighborhood is the set of 3x3 pixels centered at (x,y).

    // The neighborhoods can be accessed as 1-dimensional indices using these
    // functions.  The first four functions provide 1-dimensional indices
    // relative to any pixel location; these depend only on the image
    // dimensions.  The last four functions provide 1-dimensional indices
    // for the actual pixels in the neighborhood; no clamping is used when
    // (x,y) is on the boundary.
    void GetNeighborhood(std::array<int, 4>& nbr) const;
    void GetNeighborhood(std::array<int, 8>& nbr) const;
    void GetCorners(std::array<int, 4>& nbr) const;
    void GetFull(std::array<int, 9>& nbr) const;
    void GetNeighborhood(int x, int y, std::array<size_t, 4>& nbr) const;
    void GetNeighborhood(int x, int y, std::array<size_t, 8>& nbr) const;
    void GetCorners(int x, int y, std::array<size_t, 4>& nbr) const;
    void GetFull(int x, int y, std::array<size_t, 9>& nbr) const;

    // The neighborhoods can be accessed as 2-tuples using these functions.
    // The first four functions provide 2-tuples relative to any pixel
    // location; these depend only on the image dimensions.  The last four
    // functions provide 2-tuples for the actual pixels in the neighborhood;
    // no clamping is used when (x,y) is on the boundary.
    void GetNeighborhood(std::array<std::array<int, 2>, 4>& nbr) const;
    void GetNeighborhood(std::array<std::array<int, 2>, 8>& nbr) const;
    void GetCorners(std::array<std::array<int, 2>, 4>& nbr) const;
    void GetFull(std::array<std::array<int, 2>, 9>& nbr) const;
    void GetNeighborhood(int x, int y,
        std::array<std::array<size_t, 2>, 4>& nbr) const;
    void GetNeighborhood(int x, int y,
        std::array<std::array<size_t, 2>, 8>& nbr) const;
    void GetCorners(int x, int y,
        std::array<std::array<size_t, 2>, 4>& nbr) const;
    void GetFull(int x, int y,
        std::array<std::array<size_t, 2>, 9>& nbr) const;

    // Resize an image.  All data is lost from the original image.  The
    // function is convenient for taking a default-constructed image and
    // setting its dimension once it is known.  This avoids an irrelevant
    // memory copy that occurs if instead you were to use the statement
    // image = Image1<PixelType>(dimension0, dimension1).  The return value
    // is 'true' whenever the image is resized (reallocations occurred).
    bool Resize(int dimension0, int dimension1);

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
    void AllocatePointers();
    void DeallocatePointers();

    // Typed pointers to Image::mRawPixels.
    PixelType** mPixels;

    // Uninitialized, used in the Get(int) calls.
    PixelType mInvalidPixel;
};


template <typename PixelType>
Image2<PixelType>::~Image2()
{
    DeallocatePointers();
}

template <typename PixelType>
Image2<PixelType>::Image2()
    :
    Image(typeid(PixelType).name(), sizeof(PixelType), 2, 0, 0),
    mPixels(nullptr)
{
}

template <typename PixelType>
Image2<PixelType>::Image2(Image2 const& image)
    :
    Image(image),
    mPixels(nullptr)
{
    AllocatePointers();
}

template <typename PixelType>
Image2<PixelType>::Image2(int dimension0, int dimension1)
    :
    Image(typeid(PixelType).name(), sizeof(PixelType), 2, dimension0,
    dimension1),
    mPixels(nullptr)
{
    AllocatePointers();
}

template <typename PixelType>
Image2<PixelType>& Image2<PixelType>::operator= (Image2 const& image)
{
    bool compatible = Copy(image);
    if (!compatible)
    {
        AllocatePointers();
    }
    return *this;
}

template <typename PixelType>
void Image2<PixelType>::SetRawPixels(char* rawPixels)
{
    Image::SetRawPixels(rawPixels);
    AllocatePointers();
}

template <typename PixelType> inline
size_t Image2<PixelType>::GetIndex(int x, int y) const
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (0 <= x && x < mDimensions[0] && 0 <= y && y < mDimensions[1])
    {
        return (size_t)x + (size_t)mDimensions[0] * (size_t)y;
    }
    LogError("Invalid coordinates (" + std::to_string(x) + "," +
        std::to_string(y) + ").");
    return 0;
#else
    return (size_t)x + (size_t)mDimensions[0] * (size_t)y;
#endif
}

template <typename PixelType> inline
size_t Image2<PixelType>::GetIndex(std::array<int, 2> const& coord) const
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (0 <= coord[0] && coord[0] < mDimensions[0]
        && 0 <= coord[1] && coord[1] < mDimensions[1])
    {
        return (size_t)coord[0] + (size_t)mDimensions[0] * (size_t)coord[1];
    }
    LogError("Invalid coordinates (" + std::to_string(coord[0]) + "," +
        std::to_string(coord[1]) + ").");
    return 0;
#else
    return (size_t)coord[0] + (size_t)mDimensions[0] * (size_t)coord[1];
#endif
}

template <typename PixelType> inline
void Image2<PixelType>::GetCoordinates(size_t index, int& x, int& y) const
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (index < mNumPixels)
    {
        x = (int)(index % mDimensions[0]);
        y = (int)(index / mDimensions[0]);
    }
    else
    {
        LogError("Invalid index " + std::to_string(index) + ".");
        x = 0;
        y = 0;
    }
#else
    x = (int)(index % mDimensions[0]);
    y = (int)(index / mDimensions[0]);
#endif
}

template <typename PixelType> inline
std::array<int, 2> Image2<PixelType>::GetCoordinates(size_t index) const
{
    std::array<int, 2> coord;
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (index < mNumPixels)
    {
        coord[0] = (int)(index % mDimensions[0]);
        coord[1] = (int)(index / mDimensions[0]);
    }
    else
    {
        LogError("Invalid index " + std::to_string(index) + ".");
        coord[0] = 0;
        coord[1] = 0;
    }
#else
    coord[0] = (int)(index % mDimensions[0]);
    coord[1] = (int)(index / mDimensions[0]);
#endif
    return coord;
}

template <typename PixelType> inline
PixelType* Image2<PixelType>::GetPixels1D() const
{
    if (mPixels)
    {
        return mPixels[0];
    }
    LogError("Pixels do not exist.");
    return 0;
}

template <typename PixelType> inline
PixelType& Image2<PixelType>::operator[] (size_t i)
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (mPixels != 0 && i < mNumPixels)
    {
        return mPixels[0][i];
    }
    LogError("No pixels or invalid index " + std::to_string(i) + ".");
    return mInvalidPixel;
#else
    return mPixels[0][i];
#endif
}

template <typename PixelType> inline
PixelType const& Image2<PixelType>::operator[] (size_t i) const
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (mPixels != 0 && i < mNumPixels)
    {
        return mPixels[0][i];
    }
    LogError("No pixels or invalid index " + std::to_string(i) + ".");
    return mInvalidPixel;
#else
    return mPixels[0][i];
#endif
}

template <typename PixelType>
PixelType& Image2<PixelType>::Get(size_t i)
{
    if (mPixels)
    {
        if (i >= mNumPixels)
        {
            i = mNumPixels - 1;
        }
        return mPixels[0][i];
    }
    return mInvalidPixel;
}

template <typename PixelType>
PixelType const& Image2<PixelType>::Get(size_t i) const
{
    if (mPixels)
    {
        if (i >= mNumPixels)
        {
            i = mNumPixels - 1;
        }
        return mPixels[0][i];
    }
    return mInvalidPixel;
}

template <typename PixelType> inline
PixelType** Image2<PixelType>::GetPixels2D() const
{
    LogAssert(mPixels, "Pixels do not exist.");
    return mPixels;
}

template <typename PixelType> inline
PixelType& Image2<PixelType>::operator() (int x, int y)
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (mPixels != 0
        && 0 <= x && x < mDimensions[0]
        && 0 <= y && y < mDimensions[1])
    {
        return mPixels[y][x];
    }
    LogError("Invalid coordinates (" + std::to_string(x) + "," +
        std::to_string(y) + ").");
    return mInvalidPixel;
#else
    return mPixels[y][x];
#endif
}

template <typename PixelType> inline
PixelType const& Image2<PixelType>::operator() (int x, int y) const
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (mPixels != 0
        && 0 <= x && x < mDimensions[0]
        && 0 <= y && y < mDimensions[1])
    {
        return mPixels[y][x];
    }
    LogError("Invalid coordinates (" + std::to_string(x) + "," +
        std::to_string(y) + ").");
    return mInvalidPixel;
#else
    return mPixels[y][x];
#endif
}

template <typename PixelType> inline
PixelType& Image2<PixelType>::operator() (std::array<int, 2> const& coord)
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (mPixels != 0
        && 0 <= coord[0] && coord[0] < mDimensions[0]
        && 0 <= coord[1] && coord[1] < mDimensions[1])
    {
        return mPixels[coord[1]][coord[0]];
    }
    LogError("Invalid coordinates (" + std::to_string(coord[0]) + "," +
        std::to_string(coord[1]) + ").");
    return mInvalidPixel;
#else
    return mPixels[coord[1]][coord[0]];
#endif
}

template <typename PixelType> inline
PixelType const& Image2<PixelType>::operator() (
std::array<int, 2> const& coord) const
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (mPixels != 0
        && 0 <= coord[0] && coord[0] < mDimensions[0]
        && 0 <= coord[1] && coord[1] < mDimensions[1])
    {
        return mPixels[coord[1]][coord[0]];
    }
    LogError("Invalid coordinates (" + std::to_string(coord[0]) + "," +
        std::to_string(coord[1]) + ").");
    return mInvalidPixel;
#else
    return mPixels[coord[1]][coord[0]];
#endif
}

template <typename PixelType>
PixelType& Image2<PixelType>::Get(int x, int y)
{
    if (mPixels)
    {
        // Clamp to valid (x,y).
        if (x < 0)
        {
            x = 0;
        }
        else if (x >= mDimensions[0])
        {
            x = mDimensions[0] - 1;
        }

        if (y < 0)
        {
            y = 0;
        }
        else if (y >= mDimensions[1])
        {
            y = mDimensions[1] - 1;
        }

        return mPixels[y][x];
    }
    return mInvalidPixel;
}

template <typename PixelType>
PixelType const& Image2<PixelType>::Get(int x, int y) const
{
    if (mPixels)
    {
        // Clamp to valid (x,y).
        if (x < 0)
        {
            x = 0;
        }
        else if (x >= mDimensions[0])
        {
            x = mDimensions[0] - 1;
        }

        if (y < 0)
        {
            y = 0;
        }
        else if (y >= mDimensions[1])
        {
            y = mDimensions[1] - 1;
        }

        return mPixels[y][x];
    }
    return mInvalidPixel;
}

template <typename PixelType>
PixelType& Image2<PixelType>::Get(std::array<int, 2> coord)
{
    if (mPixels)
    {
        // Clamp to valid (x,y).
        for (int i = 0; i < 2; ++i)
        {
            if (coord[i] < 0)
            {
                coord[i] = 0;
            }
            else if (coord[i] >= mDimensions[i])
            {
                coord[i] = mDimensions[i] - 1;
            }
        }

        return mPixels[coord[1]][coord[0]];
    }
    return mInvalidPixel;
}

template <typename PixelType>
PixelType const& Image2<PixelType>::Get(std::array<int, 2> coord) const
{
    if (mPixels)
    {
        // Clamp to valid (x,y).
        for (int i = 0; i < 2; ++i)
        {
            if (coord[i] < 0)
            {
                coord[i] = 0;
            }
            else if (coord[i] >= mDimensions[i])
            {
                coord[i] = mDimensions[i] - 1;
            }
        }

        return mPixels[coord[1]][coord[0]];
    }
    return mInvalidPixel;
}

template <typename PixelType>
void Image2<PixelType>::GetNeighborhood(std::array<int, 4>& nbr) const
{
    int dim0 = mDimensions[0];
    nbr[0] = -1;        // (x-1,y)
    nbr[1] = +1;        // (x+1,y)
    nbr[2] = -dim0;     // (x,y-1)
    nbr[3] = +dim0;     // (x,y+1)
}

template <typename PixelType>
void Image2<PixelType>::GetNeighborhood(std::array<int, 8>& nbr) const
{
    int dim0 = mDimensions[0];
    nbr[0] = -1;            // (x-1,y)
    nbr[1] = +1;            // (x+1,y)
    nbr[2] = -dim0;         // (x,y-1)
    nbr[3] = +dim0;         // (x,y+1)
    nbr[4] = -1 - dim0;     // (x-1,y-1)
    nbr[5] = +1 - dim0;     // (x+1,y-1)
    nbr[6] = -1 + dim0;     // (x-1,y+1)
    nbr[7] = +1 + dim0;     // (x+1,y+1)
}

template <typename PixelType>
void Image2<PixelType>::GetCorners(std::array<int, 4>& nbr) const
{
    int dim0 = mDimensions[0];
    nbr[0] = 0;         // (x,y)
    nbr[1] = 1;         // (x+1,y)
    nbr[2] = dim0;      // (x,y+1)
    nbr[3] = dim0 + 1;  // (x+1,y+1)
}

template <typename PixelType>
void Image2<PixelType>::GetFull(std::array<int, 9>& nbr) const
{
    int dim0 = mDimensions[0];
    nbr[0] = -1 - dim0;     // (x-1,y-1)
    nbr[1] = -dim0;         // (x,y-1)
    nbr[2] = +1 - dim0;     // (x+1,y-1)
    nbr[3] = -1;            // (x-1,y)
    nbr[4] = 0;             // (x,y)
    nbr[5] = +1;            // (x+1,y)
    nbr[6] = -1 + dim0;     // (x-1,y+1)
    nbr[7] = +dim0;         // (x,y+1)
    nbr[8] = +1 + dim0;     // (x+1,y+1)
}

template <typename PixelType>
void Image2<PixelType>::GetNeighborhood(int x, int y,
    std::array<size_t, 4>& nbr) const
{
    size_t index = GetIndex(x, y);
    std::array<int, 4> inbr;
    GetNeighborhood(inbr);
    for (int i = 0; i < 4; ++i)
    {
        nbr[i] = index + inbr[i];
    }
}

template <typename PixelType>
void Image2<PixelType>::GetNeighborhood(int x, int y,
    std::array<size_t, 8>& nbr) const
{
    size_t index = GetIndex(x, y);
    std::array<int, 8> inbr;
    GetNeighborhood(inbr);
    for (int i = 0; i < 8; ++i)
    {
        nbr[i] = index + inbr[i];
    }
}

template <typename PixelType>
void Image2<PixelType>::GetCorners(int x, int y,
    std::array<size_t, 4>& nbr) const
{
    size_t index = GetIndex(x, y);
    std::array<int, 4> inbr;
    GetCorners(inbr);
    for (int i = 0; i < 4; ++i)
    {
        nbr[i] = index + inbr[i];
    }
}

template <typename PixelType>
void Image2<PixelType>::GetFull(int x, int y,
    std::array<size_t, 9>& nbr) const
{
    size_t index = GetIndex(x, y);
    std::array<int, 9> inbr;
    GetFull(inbr);
    for (int i = 0; i < 9; ++i)
    {
        nbr[i] = index + inbr[i];
    }
}

template <typename PixelType>
void Image2<PixelType>::GetNeighborhood(
    std::array<std::array<int, 2>, 4>& nbr) const
{
    nbr[0] = { { -1, 0 } };
    nbr[1] = { { +1, 0 } };
    nbr[2] = { { 0, -1 } };
    nbr[3] = { { 0, +1 } };
}

template <typename PixelType>
void Image2<PixelType>::GetNeighborhood(
    std::array<std::array<int, 2>, 8>& nbr) const
{
    nbr[0] = { { -1, -1 } };
    nbr[1] = { { 0, -1 } };
    nbr[2] = { { +1, -1 } };
    nbr[3] = { { -1, 0 } };
    nbr[4] = { { +1, 0 } };
    nbr[5] = { { -1, +1 } };
    nbr[6] = { { 0, +1 } };
    nbr[7] = { { +1, +1 } };
}

template <typename PixelType>
void Image2<PixelType>::GetCorners(
    std::array<std::array<int, 2>, 4>& nbr) const
{
    nbr[0] = { { 0, 0 } };
    nbr[1] = { { 1, 0 } };
    nbr[2] = { { 0, 1 } };
    nbr[3] = { { 1, 1 } };
}

template <typename PixelType>
void Image2<PixelType>::GetFull(std::array<std::array<int, 2>, 9>& nbr) const
{
    nbr[0] = { { -1, -1 } };
    nbr[1] = { { 0, -1 } };
    nbr[2] = { { +1, -1 } };
    nbr[3] = { { -1, 0 } };
    nbr[4] = { { 0, 0 } };
    nbr[5] = { { +1, 0 } };
    nbr[6] = { { -1, +1 } };
    nbr[7] = { { 0, +1 } };
    nbr[8] = { { +1, +1 } };
}

template <typename PixelType>
void Image2<PixelType>::GetNeighborhood(int x, int y,
    std::array<std::array<size_t, 2>, 4>& nbr) const
{
    std::array<std::array<int, 2>, 4> inbr;
    GetNeighborhood(inbr);
    for (int i = 0; i < 4; ++i)
    {
        nbr[i][0] = (size_t)x + inbr[i][0];
        nbr[i][1] = (size_t)y + inbr[i][1];
    }
}

template <typename PixelType>
void Image2<PixelType>::GetNeighborhood(int x, int y,
    std::array<std::array<size_t, 2>, 8>& nbr) const
{
    std::array<std::array<int, 2>, 8> inbr;
    GetNeighborhood(inbr);
    for (int i = 0; i < 8; ++i)
    {
        nbr[i][0] = (size_t)x + inbr[i][0];
        nbr[i][1] = (size_t)y + inbr[i][1];
    }
}

template <typename PixelType>
void Image2<PixelType>::GetCorners(int x, int y,
    std::array<std::array<size_t, 2>, 4>& nbr) const
{
    std::array<std::array<int, 2>, 4> inbr;
    GetCorners(inbr);
    for (int i = 0; i < 4; ++i)
    {
        nbr[i][0] = (size_t)x + inbr[i][0];
        nbr[i][1] = (size_t)y + inbr[i][1];
    }
}

template <typename PixelType>
void Image2<PixelType>::GetFull(int x, int y,
    std::array<std::array<size_t, 2>, 9>& nbr) const
{
    std::array<std::array<int, 2>, 9> inbr;
    GetFull(inbr);
    for (int i = 0; i < 9; ++i)
    {
        nbr[i][0] = (size_t)x + inbr[i][0];
        nbr[i][1] = (size_t)y + inbr[i][1];
    }
}

template <typename PixelType>
bool Image2<PixelType>::Resize(int dimension0, int dimension1)
{
    if (Image::Resize(2, dimension0, dimension1))
    {
        AllocatePointers();
        return true;
    }
    return false;
}

template <typename PixelType>
void Image2<PixelType>::SetAllPixels(PixelType const& value)
{
    if (mPixels)
    {
        for (size_t i = 0; i < mNumPixels; ++i)
        {
            mPixels[0][i] = value;
        }
    }
}

template <typename PixelType>
bool Image2<PixelType>::Load(std::string const& filename,
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
        AllocatePointers();
        return true;
    }

    return false;
}

template <typename PixelType>
void Image2<PixelType>::AllocatePointers()
{
    if (mPixels)
    {
        delete[] mPixels;
    }

    mPixels = new PixelType*[mDimensions[1]];
    mPixels[0] = reinterpret_cast<PixelType*>(mRawPixels);
    for (int i1 = 1; i1 < mDimensions[1]; ++i1)
    {
        int j0 = mDimensions[0] * i1;  // = mDimensions[0]*(i1 + j1), j1 = 0
        mPixels[i1] = &mPixels[0][j0];
    }
}

template <typename PixelType>
void Image2<PixelType>::DeallocatePointers()
{
    delete[] mPixels;
    mPixels = nullptr;
}


}
