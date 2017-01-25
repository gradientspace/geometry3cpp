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
class Image3 : public Image
{
public:
    // The default constructor creates a null image.  You can resize the image
    // later with an explicit Resize call, an assignment, or by loading the
    // image from disk.
    virtual ~Image3();
    Image3();

    // Copy the input image using the assignment operator.
    Image3(Image3 const& image);

    // The input dimensions must be positive; otherwise, a null image is
    // created.
    Image3(int dimension0, int dimension1, int dimension2);

    // If the input image is compatible with 'this', a copy of the input
    // image data occurs.  If the input image is not compatible, 'this' is
    // recreated to be a copy of 'image'.
    Image3& operator= (Image3 const& image);

    // The input array must have the correct number of pixels as determined by
    // the image parameters.  Use at your own risk, because we cannot verify
    // the compatibility.
    virtual void SetRawPixels(char* rawPixels);

    // Conversion between 1-dimensional indices and 2-dimensional coordinates.
    inline size_t GetIndex(int x, int y, int z) const;
    inline size_t GetIndex(std::array<int, 3> const& coord) const;
    inline void GetCoordinates(size_t index, int& x, int& y, int& z) const;
    inline std::array<int, 3> GetCoordinates(size_t index) const;

    // Access the data as a 1-dimensional array.  The operator[] functions
    // test for valid i in debug configurations and assert on invalid i.  The
    // Get() functions test for valid i and clamp when invalid (debug and
    // release); these functions cannot fail.
    inline PixelType* GetPixels1D() const;
    inline PixelType& operator[] (size_t i);
    inline PixelType const& operator[] (size_t i) const;
    PixelType& Get(size_t i);
    PixelType const& Get(size_t i) const;

    // Access the data as a 3-dimensional array.  Pixel (x,y,z) is accessed
    // as "pixels3D[z][y][x]".  The operator() functions test for valid
    // (x,y,z) in debug configurations and assert on invalid (x,y,z).  The
    // Get() functions test for valid (x,y,z) and clamp when invalid (debug
    // and release); these functions cannot fail.
    inline PixelType*** GetPixels3D() const;
    inline PixelType& operator() (int x, int y, int z);
    inline PixelType const& operator() (int x, int y, int z) const;
    inline PixelType& operator() (std::array<int, 3> const& coord);
    inline PixelType const& operator() (std::array<int,3> const& coord) const;
    inline PixelType& Get(int x, int y, int z);
    inline PixelType const& Get(int x, int y, int z) const;
    inline PixelType& Get(std::array<int, 3> coord);
    inline PixelType const& Get(std::array<int, 3> coord) const;

    // In the following discussion, u, v and w are in {-1,1}.  Given a voxel
    // (x,y,z), the 6-connected neighbors have relative offsets (u,0,0),
    // (0,v,0), and (0,0,w).  The 18-connected neighbors include the
    // 6-connected neighbors and have additional relative offsets (u,v,0),
    // (u,0,w), and (0,v,w).  The 26-connected neighbors include the
    // 18-connected neighbors and have additional relative offsets (u,v,w).
    // The corner neighbors have offsets (0,0,0), (1,0,0), (0,1,0), (1,1,0),
    // (0,0,1), (1,0,1), (0,1,1), and (1,1,1) in that order.  The full
    // neighborhood is the set of 3x3x3 pixels centered at (x,y).

    // The neighborhoods can be accessed as 1-dimensional indices using these
    // functions.  The first five functions provide 1-dimensional indices
    // relative to any voxel location; these depend only on the image
    // dimensions.  The last five functions provide 1-dimensional indices
    // for the actual voxels in the neighborhood; no clamping is used when
    // (x,y,z) is on the boundary.
    void GetNeighborhood(std::array<int, 6>& nbr) const;
    void GetNeighborhood(std::array<int, 18>& nbr) const;
    void GetNeighborhood(std::array<int, 26>& nbr) const;
    void GetCorners(std::array<int, 8>& nbr) const;
    void GetFull(std::array<int, 27>& nbr) const;
    void GetNeighborhood(int x, int y, int z,
        std::array<size_t, 6>& nbr) const;
    void GetNeighborhood(int x, int y, int z,
        std::array<size_t, 18>& nbr) const;
    void GetNeighborhood(int x, int y, int z,
        std::array<size_t, 26>& nbr) const;
    void GetCorners(int x, int y, int z,
        std::array<size_t, 8>& nbr) const;
    void GetFull(int x, int y, int z,
        std::array<size_t, 27>& nbr) const;

    // The neighborhoods can be accessed as 3-tuples using these functions.
    // The first five functions provide 3-tuples relative to any voxel
    // location; these depend only on the image dimensions.  The last five
    // functions provide 3-tuples for the actual voxels in the neighborhood;
    // no clamping is used when (x,y,z) is on the boundary.
    void GetNeighborhood(std::array<std::array<int, 3>, 6>& nbr) const;
    void GetNeighborhood(std::array<std::array<int, 3>, 18>& nbr) const;
    void GetNeighborhood(std::array<std::array<int, 3>, 26>& nbr) const;
    void GetCorners(std::array<std::array<int, 3>, 8>& nbr) const;
    void GetFull(std::array<std::array<int, 3>, 27>& nbr) const;
    void GetNeighborhood(int x, int y, int z,
        std::array<std::array<size_t, 3>, 6>& nbr) const;
    void GetNeighborhood(int x, int y, int z,
        std::array<std::array<size_t, 3>, 18>& nbr) const;
    void GetNeighborhood(int x, int y, int z,
        std::array<std::array<size_t, 3>, 26>& nbr) const;
    void GetCorners(int x, int y, int z,
        std::array<std::array<size_t, 3>, 8>& nbr) const;
    void GetFull(int x, int y, int z,
        std::array<std::array<size_t, 3>, 27>& nbr) const;

    // Resize an image.  All data is lost from the original image.  The
    // function is convenient for taking a default-constructed image and
    // setting its dimension once it is known.  This avoids an irrelevant
    // memory copy that occurs if instead you were to use the statement
    // image = Image1<PixelType>(dimension0, dimension1, dimension2).  The
    // return value is 'true' whenever the image is resized (reallocations
    // occurred).
    bool Resize(int dimension0, int dimension1, int dimension2);

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
    PixelType*** mPixels;

    // Uninitialized, used in the Get(int) calls.
    PixelType mInvalidPixel;
};


template <typename PixelType>
Image3<PixelType>::~Image3()
{
    DeallocatePointers();
}

template <typename PixelType>
Image3<PixelType>::Image3()
    :
    Image(typeid(PixelType).name(), sizeof(PixelType), 3, 0, 0, 0),
    mPixels(nullptr)
{
}

template <typename PixelType>
Image3<PixelType>::Image3(Image3 const& image)
    :
    Image(image),
    mPixels(nullptr)
{
    AllocatePointers();
}

template <typename PixelType>
Image3<PixelType>::Image3(int dimension0, int dimension1, int dimension2)
    :
    Image(typeid(PixelType).name(), sizeof(PixelType), 3, dimension0,
    dimension1, dimension2),
    mPixels(nullptr)
{
    AllocatePointers();
}

template <typename PixelType>
Image3<PixelType>& Image3<PixelType>::operator= (Image3 const& image)
{
    bool compatible = Copy(image);
    if (!compatible)
    {
        AllocatePointers();
    }
    return *this;
}

template <typename PixelType>
void Image3<PixelType>::SetRawPixels(char* rawPixels)
{
    Image::SetRawPixels(rawPixels);
    AllocatePointers();
}

template <typename PixelType> inline
size_t Image3<PixelType>::GetIndex(int x, int y, int z) const
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (0 <= x && x < mDimensions[0]
        && 0 <= y && y < mDimensions[1]
        && 0 <= z && z < mDimensions[2])
    {
        return (size_t)x + (size_t)mDimensions[0] * ((size_t)y +
            (size_t)mDimensions[1] * (size_t)z);
    }
    LogError("Invalid coordinates (" + std::to_string(x) + "," +
        std::to_string(y) + "," + std::to_string(z) + ").");
    return 0;
#else
    return (size_t)x + (size_t)mDimensions[0] * ((size_t)y +
        (size_t)mDimensions[1] * (size_t)z);
#endif
}

template <typename PixelType> inline
size_t Image3<PixelType>::GetIndex(std::array<int, 3> const& coord) const
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (0 <= coord[0] && coord[0] < mDimensions[0]
        && 0 <= coord[1] && coord[1] < mDimensions[1]
        && 0 <= coord[2] && coord[2] < mDimensions[2])
    {
        return (size_t)coord[0] + (size_t)mDimensions[0] * ((size_t)coord[1] +
            (size_t)mDimensions[1] * (size_t)coord[2]);
    }
    LogError("Invalid coordinates (" + std::to_string(coord[0]) + "," +
        std::to_string(coord[1]) + "," + std::to_string(coord[2]) + ").");
    return 0;
#else
    return (size_t)coord[0] + (size_t)mDimensions[0] * ((size_t)coord[1] +
        (size_t)mDimensions[1] * (size_t)coord[2]);
#endif
}

template <typename PixelType> inline
void Image3<PixelType>::GetCoordinates(size_t index, int& x, int& y, int& z)
const
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (index < mNumPixels)
    {
        x = (int)(index % mDimensions[0]);
        index /= mDimensions[0];
        y = (int)(index % mDimensions[1]);
        z = (int)(index / mDimensions[1]);
    }
    else
    {
        LogError("Invalid index " + std::to_string(index) + ".");
        x = 0;
        y = 0;
        z = 0;
    }
#else
    x = (int)(index % mDimensions[0]);
    index /= mDimensions[0];
    y = (int)(index % mDimensions[1]);
    z = (int)(index / mDimensions[1]);
#endif
}

template <typename PixelType> inline
std::array<int, 3> Image3<PixelType>::GetCoordinates(size_t index) const
{
    std::array<int, 3> coord;
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (index < mNumPixels)
    {
        coord[0] = (int)(index % mDimensions[0]);
        index /= mDimensions[0];
        coord[1] = (int)(index % mDimensions[1]);
        coord[2] = (int)(index / mDimensions[1]);
    }
    else
    {
        LogError("Invalid index " + std::to_string(index) + ".");
        coord[0] = 0;
        coord[1] = 0;
        coord[2] = 0;
    }
#else
    coord[0] = (int)(index % mDimensions[0]);
    index /= mDimensions[0];
    coord[1] = (int)(index % mDimensions[1]);
    coord[2] = (int)(index / mDimensions[1]);
#endif
    return coord;
}

template <typename PixelType> inline
PixelType* Image3<PixelType>::GetPixels1D() const
{
    if (mPixels)
    {
        return mPixels[0][0];
    }
    LogError("Pixels do not exist.");
    return 0;
}

template <typename PixelType> inline
PixelType& Image3<PixelType>::operator[] (size_t i)
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (mPixels != 0 && i < mNumPixels)
    {
        return mPixels[0][0][i];
    }
    LogError("No pixels or invalid index " + std::to_string(i) + ".");
    return mInvalidPixel;
#else
    return mPixels[0][0][i];
#endif
}

template <typename PixelType> inline
PixelType const& Image3<PixelType>::operator[] (size_t i) const
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (mPixels != 0 && i < mNumPixels)
    {
        return mPixels[0][0][i];
    }
    LogError("No pixels or invalid index " + std::to_string(i) + ".");
    return mInvalidPixel;
#else
    return mPixels[0][0][i];
#endif
}

template <typename PixelType>
PixelType& Image3<PixelType>::Get(size_t i)
{
    if (mPixels)
    {
        if (i >= mNumPixels)
        {
            i = mNumPixels - 1;
        }
        return mPixels[0][0][i];
    }
    return mInvalidPixel;
}

template <typename PixelType>
PixelType const& Image3<PixelType>::Get(size_t i) const
{
    if (mPixels)
    {
        if (i >= mNumPixels)
        {
            i = mNumPixels - 1;
        }
        return mPixels[0][0][i];
    }
    return mInvalidPixel;
}

template <typename PixelType> inline
PixelType*** Image3<PixelType>::GetPixels3D() const
{
    LogAssert(mPixels, "Pixels do not exist.");
    return mPixels;
}

template <typename PixelType> inline
PixelType& Image3<PixelType>::operator() (int x, int y, int z)
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (mPixels != 0
        && 0 <= x && x < mDimensions[0]
        && 0 <= y && y < mDimensions[1]
        && 0 <= z && z < mDimensions[2])
    {
        return mPixels[z][y][x];
    }
    LogError("Invalid coordinates (" + std::to_string(x) + "," +
        std::to_string(y) + "," + std::to_string(z) + ").");
    return mInvalidPixel;
#else
    return mPixels[z][y][x];
#endif
}

template <typename PixelType> inline
PixelType const& Image3<PixelType>::operator() (int x, int y, int z) const
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (mPixels != 0
        && 0 <= x && x < mDimensions[0]
        && 0 <= y && y < mDimensions[1]
        && 0 <= z && z < mDimensions[2])
    {
        return mPixels[z][y][x];
    }
    LogError("Invalid coordinates (" + std::to_string(x) + "," +
        std::to_string(y) + "," + std::to_string(z) + ").");
    return mInvalidPixel;
#else
    return mPixels[z][y][x];
#endif
}

template <typename PixelType> inline
PixelType& Image3<PixelType>::operator() (std::array<int, 3> const& coord)
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (mPixels != 0
        && 0 <= coord[0] && coord[0] < mDimensions[0]
        && 0 <= coord[1] && coord[1] < mDimensions[1]
        && 0 <= coord[2] && coord[2] < mDimensions[2])
    {
        return mPixels[coord[2]][coord[1]][coord[0]];
    }
    LogError("Invalid coordinates (" + std::to_string(coord[0]) + "," +
        std::to_string(coord[1]) + "," + std::to_string(coord[2]) + ").");
    return mInvalidPixel;
#else
    return mPixels[coord[2]][coord[1]][coord[0]];
#endif
}

template <typename PixelType> inline
PixelType const& Image3<PixelType>::operator() (
std::array<int, 3> const& coord) const
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (mPixels != 0
        && 0 <= coord[0] && coord[0] < mDimensions[0]
        && 0 <= coord[1] && coord[1] < mDimensions[1]
        && 0 <= coord[2] && coord[2] < mDimensions[2])
    {
        return mPixels[coord[2]][coord[1]][coord[0]];
    }
    LogError("Invalid coordinates (" + std::to_string(coord[0]) + "," +
        std::to_string(coord[1]) + "," + std::to_string(coord[2]) + ").");
    return mInvalidPixel;
#else
    return mPixels[coord[2]][coord[1]][coord[0]];
#endif
}

template <typename PixelType>
PixelType& Image3<PixelType>::Get(int x, int y, int z)
{
    if (mPixels)
    {
        // Clamp to valid (x,y,z).
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

        if (z < 0)
        {
            z = 0;
        }
        else if (z >= mDimensions[2])
        {
            z = mDimensions[2] - 1;
        }

        return mPixels[z][y][x];
    }
    return mInvalidPixel;
}

template <typename PixelType>
PixelType const& Image3<PixelType>::Get(int x, int y, int z) const
{
    if (mPixels)
    {
        // Clamp to valid (x,y,z).
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

        if (z < 0)
        {
            z = 0;
        }
        else if (z >= mDimensions[2])
        {
            z = mDimensions[2] - 1;
        }

        return mPixels[z][y][x];
    }
    return mInvalidPixel;
}

template <typename PixelType>
PixelType& Image3<PixelType>::Get(std::array<int, 3> coord)
{
    if (mPixels)
    {
        // Clamp to valid (x,y,z).
        for (int i = 0; i < 3; ++i)
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

        return mPixels[coord[2]][coord[1]][coord[0]];
    }
    return mInvalidPixel;
}

template <typename PixelType>
PixelType const& Image3<PixelType>::Get(std::array<int, 3> coord) const
{
    if (mPixels)
    {
        // Clamp to valid (x,y,z).
        for (int i = 0; i < 3; ++i)
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

        return mPixels[coord[2]][coord[1]][coord[0]];
    }
    return mInvalidPixel;
}

template <typename PixelType>
void Image3<PixelType>::GetNeighborhood(std::array<int, 6>& nbr) const
{
    int dim0 = mDimensions[0];
    int dim01 = mDimensions[0] * mDimensions[1];
    nbr[0] = -1;        // (x-1,y,z)
    nbr[1] = +1;        // (x+1,y,z)
    nbr[2] = -dim0;     // (x,y-1,z)
    nbr[3] = +dim0;     // (x,y+1,z)
    nbr[4] = -dim01;    // (x,y,z-1)
    nbr[5] = +dim01;    // (x,y,z+1)
}

template <typename PixelType>
void Image3<PixelType>::GetNeighborhood(std::array<int, 18>& nbr) const
{
    int dim0 = mDimensions[0];
    int dim01 = mDimensions[0] * mDimensions[1];
    nbr[0] = -1;                // (x-1,y,z)
    nbr[1] = +1;                // (x+1,y,z)
    nbr[2] = -dim0;             // (x,y-1,z)
    nbr[3] = +dim0;             // (x,y+1,z)
    nbr[4] = -dim01;            // (x,y,z-1)
    nbr[5] = +dim01;            // (x,y,z+1)
    nbr[6] = -1 - dim0;         // (x-1,y-1,z)
    nbr[7] = +1 - dim0;         // (x+1,y-1,z)
    nbr[8] = -1 + dim0;         // (x-1,y+1,z)
    nbr[9] = +1 + dim0;         // (x+1,y+1,z)
    nbr[10] = -1 + dim01;       // (x-1,y,z+1)
    nbr[11] = +1 + dim01;       // (x+1,y,z+1)
    nbr[12] = -dim0 + dim01;    // (x,y-1,z+1)
    nbr[13] = +dim0 + dim01;    // (x,y+1,z+1)
    nbr[14] = -1 - dim01;       // (x-1,y,z-1)
    nbr[15] = +1 - dim01;       // (x+1,y,z-1)
    nbr[16] = -dim0 - dim01;    // (x,y-1,z-1)
    nbr[17] = +dim0 - dim01;    // (x,y+1,z-1)
}

template <typename PixelType>
void Image3<PixelType>::GetNeighborhood(std::array<int, 26>& nbr) const
{
    int dim0 = mDimensions[0];
    int dim01 = mDimensions[0] * mDimensions[1];
    nbr[0] = -1;                    // (x-1,y,z)
    nbr[1] = +1;                    // (x+1,y,z)
    nbr[2] = -dim0;                 // (x,y-1,z)
    nbr[3] = +dim0;                 // (x,y+1,z)
    nbr[4] = -dim01;                // (x,y,z-1)
    nbr[5] = +dim01;                // (x,y,z+1)
    nbr[6] = -1 - dim0;             // (x-1,y-1,z)
    nbr[7] = +1 - dim0;             // (x+1,y-1,z)
    nbr[8] = -1 + dim0;             // (x-1,y+1,z)
    nbr[9] = +1 + dim0;             // (x+1,y+1,z)
    nbr[10] = -1 + dim01;           // (x-1,y,z+1)
    nbr[11] = +1 + dim01;           // (x+1,y,z+1)
    nbr[12] = -dim0 + dim01;        // (x,y-1,z+1)
    nbr[13] = +dim0 + dim01;        // (x,y+1,z+1)
    nbr[14] = -1 - dim01;           // (x-1,y,z-1)
    nbr[15] = +1 - dim01;           // (x+1,y,z-1)
    nbr[16] = -dim0 - dim01;        // (x,y-1,z-1)
    nbr[17] = +dim0 - dim01;        // (x,y+1,z-1)
    nbr[18] = -1 - dim0 - dim01;    // (x-1,y-1,z-1)
    nbr[19] = +1 - dim0 - dim01;    // (x+1,y-1,z-1)
    nbr[20] = -1 + dim0 - dim01;    // (x-1,y+1,z-1)
    nbr[21] = +1 + dim0 - dim01;    // (x+1,y+1,z-1)
    nbr[22] = -1 - dim0 + dim01;    // (x-1,y-1,z+1)
    nbr[23] = +1 - dim0 + dim01;    // (x+1,y-1,z+1)
    nbr[24] = -1 + dim0 + dim01;    // (x-1,y+1,z+1)
    nbr[25] = +1 + dim0 + dim01;    // (x+1,y+1,z+1)
}

template <typename PixelType>
void Image3<PixelType>::GetCorners(std::array<int, 8>& nbr) const
{
    int dim0 = mDimensions[0];
    int dim01 = mDimensions[0] * mDimensions[1];
    nbr[0] = 0;                 // (x,y,z)
    nbr[1] = 1;                 // (x+1,y,z)
    nbr[2] = dim0;              // (x,y+1,z)
    nbr[3] = dim0 + 1;          // (x+1,y+1,z)
    nbr[4] = dim01;             // (x,y,z+1)
    nbr[5] = dim01 + 1;         // (x+1,y,z+1)
    nbr[6] = dim01 + dim0;      // (x,y+1,z+1)
    nbr[7] = dim01 + dim0 + 1;  // (x+1,y+1,z+1)
}

template <typename PixelType>
void Image3<PixelType>::GetFull(std::array<int, 27>& nbr) const
{
    int dim0 = mDimensions[0];
    int dim01 = mDimensions[0] * mDimensions[1];
    nbr[0] = -1 - dim0 - dim01;     // (x-1,y-1,z-1)
    nbr[1] = -dim0 - dim01;         // (x,  y-1,z-1)
    nbr[2] = +1 - dim0 - dim01;     // (x+1,y-1,z-1)
    nbr[3] = -1 - dim01;            // (x-1,y,  z-1)
    nbr[4] = -dim01;                // (x,  y,  z-1)
    nbr[5] = +1 - dim01;            // (x+1,y,  z-1)
    nbr[6] = -1 + dim0 - dim01;     // (x-1,y+1,z-1)
    nbr[7] = +dim0 - dim01;         // (x,  y+1,z-1)
    nbr[8] = +1 + dim0 - dim01;     // (x+1,y+1,z-1)
    nbr[9] = -1 - dim0;             // (x-1,y-1,z)
    nbr[10] = -dim0;                // (x,  y-1,z)
    nbr[11] = +1 - dim0;            // (x+1,y-1,z)
    nbr[12] = -1;                   // (x-1,y,  z)
    nbr[13] = 0;                    // (x,  y,  z)
    nbr[14] = +1;                   // (x+1,y,  z)
    nbr[15] = -1 + dim0;            // (x-1,y+1,z)
    nbr[16] = +dim0;                // (x,  y+1,z)
    nbr[17] = +1 + dim0;            // (x+1,y+1,z)
    nbr[18] = -1 - dim0 + dim01;    // (x-1,y-1,z+1)
    nbr[19] = -dim0 + dim01;        // (x,  y-1,z+1)
    nbr[20] = +1 - dim0 + dim01;    // (x+1,y-1,z+1)
    nbr[21] = -1 + dim01;           // (x-1,y,  z+1)
    nbr[22] = +dim01;               // (x,  y,  z+1)
    nbr[23] = +1 + dim01;           // (x+1,y,  z+1)
    nbr[24] = -1 + dim0 + dim01;    // (x-1,y+1,z+1)
    nbr[25] = +dim0 + dim01;        // (x,  y+1,z+1)
    nbr[26] = +1 + dim0 + dim01;    // (x+1,y+1,z+1)
}

template <typename PixelType>
void Image3<PixelType>::GetNeighborhood(int x, int y, int z,
    std::array<size_t, 6>& nbr) const
{
    size_t index = GetIndex(x, y, z);
    std::array<int, 6> inbr;
    GetNeighborhood(inbr);
    for (int i = 0; i < 6; ++i)
    {
        nbr[i] = index + inbr[i];
    }
}

template <typename PixelType>
void Image3<PixelType>::GetNeighborhood(int x, int y, int z,
    std::array<size_t, 18>& nbr) const
{
    size_t index = GetIndex(x, y, z);
    std::array<int, 18> inbr;
    GetNeighborhood(inbr);
    for (int i = 0; i < 18; ++i)
    {
        nbr[i] = index + inbr[i];
    }
}

template <typename PixelType>
void Image3<PixelType>::GetNeighborhood(int x, int y, int z,
    std::array<size_t, 26>& nbr) const
{
    size_t index = GetIndex(x, y, z);
    std::array<int, 26> inbr;
    GetNeighborhood(inbr);
    for (int i = 0; i < 26; ++i)
    {
        nbr[i] = index + inbr[i];
    }
}

template <typename PixelType>
void Image3<PixelType>::GetCorners(int x, int y, int z,
    std::array<size_t, 8>& nbr) const
{
    size_t index = GetIndex(x, y, z);
    std::array<int, 8> inbr;
    GetCorners(inbr);
    for (int i = 0; i < 8; ++i)
    {
        nbr[i] = index + inbr[i];
    }
}

template <typename PixelType>
void Image3<PixelType>::GetFull(int x, int y, int z,
    std::array<size_t, 27>& nbr) const
{
    size_t index = GetIndex(x, y, z);
    std::array<int, 27> inbr;
    GetFull(inbr);
    for (int i = 0; i < 27; ++i)
    {
        nbr[i] = index + inbr[i];
    }
}

template <typename PixelType>
void Image3<PixelType>::GetNeighborhood(
    std::array<std::array<int, 3>, 6>& nbr) const
{
    nbr[0] = { { -1, 0, 0 } };
    nbr[1] = { { +1, 0, 0 } };
    nbr[2] = { { 0, -1, 0 } };
    nbr[3] = { { 0, +1, 0 } };
    nbr[4] = { { 0, 0, -1 } };
    nbr[5] = { { 0, 0, +1 } };
}

template <typename PixelType>
void Image3<PixelType>::GetNeighborhood(
    std::array<std::array<int, 3>, 18>& nbr) const
{
    nbr[0] = { { -1, 0, 0 } };
    nbr[1] = { { +1, 0, 0 } };
    nbr[2] = { { 0, -1, 0 } };
    nbr[3] = { { 0, +1, 0 } };
    nbr[4] = { { 0, 0, -1 } };
    nbr[5] = { { 0, 0, +1 } };
    nbr[6] = { { -1, -1, 0 } };
    nbr[7] = { { +1, -1, 0 } };
    nbr[8] = { { -1, +1, 0 } };
    nbr[9] = { { +1, +1, 0 } };
    nbr[10] = { { -1, 0, +1 } };
    nbr[11] = { { +1, 0, +1 } };
    nbr[12] = { { 0, -1, +1 } };
    nbr[13] = { { 0, +1, +1 } };
    nbr[14] = { { -1, 0, -1 } };
    nbr[15] = { { +1, 0, -1 } };
    nbr[16] = { { 0, -1, -1 } };
    nbr[17] = { { 0, +1, -1 } };
}

template <typename PixelType>
void Image3<PixelType>::GetNeighborhood(
    std::array<std::array<int, 3>, 26>& nbr) const
{
    nbr[0] = { { -1, 0, 0 } };
    nbr[1] = { { +1, 0, 0 } };
    nbr[2] = { { 0, -1, 0 } };
    nbr[3] = { { 0, +1, 0 } };
    nbr[4] = { { 0, 0, -1 } };
    nbr[5] = { { 0, 0, +1 } };
    nbr[6] = { { -1, -1, 0 } };
    nbr[7] = { { +1, -1, 0 } };
    nbr[8] = { { -1, +1, 0 } };
    nbr[9] = { { +1, +1, 0 } };
    nbr[10] = { { -1, 0, +1 } };
    nbr[11] = { { +1, 0, +1 } };
    nbr[12] = { { 0, -1, +1 } };
    nbr[13] = { { 0, +1, +1 } };
    nbr[14] = { { -1, 0, -1 } };
    nbr[15] = { { +1, 0, -1 } };
    nbr[16] = { { 0, -1, -1 } };
    nbr[17] = { { 0, +1, -1 } };
    nbr[18] = { { -1, -1, -1 } };
    nbr[19] = { { +1, -1, -1 } };
    nbr[20] = { { -1, +1, -1 } };
    nbr[21] = { { +1, +1, -1 } };
    nbr[22] = { { -1, -1, +1 } };
    nbr[23] = { { +1, -1, +1 } };
    nbr[24] = { { -1, +1, +1 } };
    nbr[25] = { { +1, +1, +1 } };
}

template <typename PixelType>
void Image3<PixelType>::GetCorners(
    std::array<std::array<int, 3>, 8>& nbr) const
{
    nbr[0] = { { 0, 0, 0 } };
    nbr[1] = { { 1, 0, 0 } };
    nbr[2] = { { 0, 1, 0 } };
    nbr[3] = { { 1, 1, 0 } };
    nbr[4] = { { 0, 0, 1 } };
    nbr[5] = { { 1, 0, 1 } };
    nbr[6] = { { 0, 1, 1 } };
    nbr[7] = { { 1, 1, 1 } };
}

template <typename PixelType>
void Image3<PixelType>::GetFull(std::array<std::array<int, 3>, 27>& nbr) const
{
    nbr[0] = { { -1, -1, -1 } };
    nbr[1] = { { 0, -1, -1 } };
    nbr[2] = { { +1, -1, -1 } };
    nbr[3] = { { -1, 0, -1 } };
    nbr[4] = { { 0, 0, -1 } };
    nbr[5] = { { +1, 0, -1 } };
    nbr[6] = { { -1, +1, -1 } };
    nbr[7] = { { 0, +1, -1 } };
    nbr[8] = { { +1, +1, -1 } };
    nbr[9] = { { -1, -1, 0 } };
    nbr[10] = { { 0, -1, 0 } };
    nbr[11] = { { +1, -1, 0 } };
    nbr[12] = { { -1, 0, 0 } };
    nbr[13] = { { 0, 0, 0 } };
    nbr[14] = { { +1, 0, 0 } };
    nbr[15] = { { -1, +1, 0 } };
    nbr[16] = { { 0, +1, 0 } };
    nbr[17] = { { +1, +1, 0 } };
    nbr[18] = { { -1, -1, +1 } };
    nbr[19] = { { 0, -1, +1 } };
    nbr[20] = { { +1, -1, +1 } };
    nbr[21] = { { -1, 0, +1 } };
    nbr[22] = { { 0, 0, +1 } };
    nbr[23] = { { +1, 0, +1 } };
    nbr[24] = { { -1, +1, +1 } };
    nbr[25] = { { 0, +1, +1 } };
    nbr[26] = { { +1, +1, +1 } };
}

template <typename PixelType>
void Image3<PixelType>::GetNeighborhood(int x, int y, int z,
    std::array<std::array<size_t, 3>, 6>& nbr) const
{
    std::array<std::array<int, 3>, 6> inbr;
    GetNeighborhood(inbr);
    for (int i = 0; i < 6; ++i)
    {
        nbr[i][0] = (size_t)x + inbr[i][0];
        nbr[i][1] = (size_t)y + inbr[i][1];
        nbr[i][2] = (size_t)z + inbr[i][2];
    }
}

template <typename PixelType>
void Image3<PixelType>::GetNeighborhood(int x, int y, int z,
    std::array<std::array<size_t, 3>, 18>& nbr) const
{
    std::array<std::array<int, 3>, 18> inbr;
    GetNeighborhood(inbr);
    for (int i = 0; i < 18; ++i)
    {
        nbr[i][0] = (size_t)x + inbr[i][0];
        nbr[i][1] = (size_t)y + inbr[i][1];
        nbr[i][2] = (size_t)z + inbr[i][2];
    }
}

template <typename PixelType>
void Image3<PixelType>::GetNeighborhood(int x, int y, int z,
    std::array<std::array<size_t, 3>, 26>& nbr) const
{
    std::array<std::array<int, 3>, 26> inbr;
    GetNeighborhood(inbr);
    for (int i = 0; i < 26; ++i)
    {
        nbr[i][0] = (size_t)x + inbr[i][0];
        nbr[i][1] = (size_t)y + inbr[i][1];
        nbr[i][2] = (size_t)z + inbr[i][2];
    }
}

template <typename PixelType>
void Image3<PixelType>::GetCorners(int x, int y, int z,
    std::array<std::array<size_t, 3>, 8>& nbr) const
{
    std::array<std::array<int, 3>, 8> inbr;
    GetCorners(inbr);
    for (int i = 0; i < 8; ++i)
    {
        nbr[i][0] = (size_t)x + inbr[i][0];
        nbr[i][1] = (size_t)y + inbr[i][1];
        nbr[i][2] = (size_t)z + inbr[i][2];
    }
}

template <typename PixelType>
void Image3<PixelType>::GetFull(int x, int y, int z,
    std::array<std::array<size_t, 3>, 27>& nbr) const
{
    std::array<std::array<int, 3>, 27> inbr;
    GetFull(inbr);
    for (int i = 0; i < 27; ++i)
    {
        nbr[i][0] = (size_t)x + inbr[i][0];
        nbr[i][1] = (size_t)y + inbr[i][1];
        nbr[i][2] = (size_t)z + inbr[i][2];
    }
}

template <typename PixelType>
bool Image3<PixelType>::Resize(int dimension0, int dimension1,
    int dimension2)
{
    if (Image::Resize(3, dimension0, dimension1, dimension2))
    {
        AllocatePointers();
        return true;
    }
    return false;
}

template <typename PixelType>
void Image3<PixelType>::SetAllPixels(PixelType const& value)
{
    if (mPixels)
    {
        for (size_t i = 0; i < mNumPixels; ++i)
        {
            mPixels[0][0][i] = value;
        }
    }
}

template <typename PixelType>
bool Image3<PixelType>::Load(std::string const& filename,
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
void Image3<PixelType>::AllocatePointers()
{
    if (mPixels)
    {
        delete[] mPixels[0];
        delete[] mPixels;
    }

    mPixels = new PixelType**[mDimensions[2]];
    mPixels[0] = new PixelType*[mDimensions[1] * mDimensions[2]];
    mPixels[0][0] = reinterpret_cast<PixelType*>(mRawPixels);
    for (int i2 = 0; i2 < mDimensions[2]; ++i2)
    {
        int j1 = mDimensions[1] * i2;  // = mDimensions[1]*(i2 + j2), j2 = 0
        mPixels[i2] = &mPixels[0][j1];
        for (int i1 = 0; i1 < mDimensions[1]; ++i1)
        {
            int j0 = mDimensions[0] * (i1 + j1);
            mPixels[i2][i1] = &mPixels[0][0][j0];
        }
    }
}

template <typename PixelType>
void Image3<PixelType>::DeallocatePointers()
{
    if (mPixels)
    {
        delete[] mPixels[0];
        delete[] mPixels;
        mPixels = nullptr;
    }
}


}
