// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <LowLevel/GteLogger.h>
#include <fstream>
#include <string>
#include <vector>

// Testing for out-of-range indices in the image classes Image, Image1,
// Image2, and Image3.
//#define GTE_IMAGICS_ASSERT_ON_INVALID_INDEX

namespace gte
{

// The pixel data type must be "plain old data" (POD) and cannot have side
// effects from construction or destruction.  Classes derived from Image
// cannot have construction or destruction side effects.  These requirements
// allow the loading of an image independent of the pixel data type and
// derived class via the Load function.  The loaded image and raw pixels can
// be typecast without danger of an incorrect memory layout.

class GTE_IMPEXP Image
{
public:
    // Construction and destruction.
    virtual ~Image();
    Image();
    Image(Image const& image);
    Image(std::string const& pixelType, size_t pixelSize, int numDimensions,
        ...);

    // Assignment.
    Image& operator=(Image const& image);

    // Access to image data.
    inline std::string GetPixelType() const;
    inline size_t GetPixelSize() const;
    inline int GetNumDimensions() const;
    inline int const* GetDimensions() const;
    inline int GetDimension(int d) const;
    inline size_t const* GetOffsets() const;
    inline size_t GetOffset(int d) const;
    inline size_t GetNumPixels() const;
    inline char const* GetRawPixels() const;
    inline char* GetRawPixels();
    inline bool IsNullImage() const;

    // The input array must have the correct number of pixels as determined by
    // the image parameters.  Use at your own risk, because we cannot verify
    // the compatibility.  The base class sets mRawPixels.  The derived
    // classes must additionally process the multidimensional array pointers.
    virtual void SetRawPixels(char* rawPixels);

    // Access to metadata for tools that process images.  Per-pixel metadata
    // is allocated on demand (when SetPixelMetaData is called the first
    // time).  The metadata is streamed.  The ClearPixelMetaData function
    // sets all pixel metadata to "" (when the metadata exists).
    inline void SetImageMetaData(std::string const& metadata);
    inline std::string const& GetImageMetaData() const;
    void CreatePixelMetaData();
    void DestroyPixelMetaData();
    inline bool HasPixelMetaData() const;
    inline void SetPixelMetaData(size_t i, std::string const& metadata);
    inline std::string const& GetPixelMetaData(size_t i) const;
    void ClearPixelMetaData();

    // Conversions between n-dim and 1-dim structures.  The 'coord' arrays
    // must have GetNumDimensions() elements.
    size_t GetIndex(int const* coord) const;
    void GetCoordinates(size_t index, int* coord) const;

    // Set all image data to zero (using a memset).  The metadata is
    // unaffected by this call.
    void ClearPixels();

    // Streaming.  The LoadHeader function is useful for tools that can
    // load an image of any dimension and pixel type, and then want to
    // typecast to ImageD<T>, where D is the number of dimensions and T is the
    // pixel type.
    static bool LoadHeader(std::ifstream& input, std::string& pixelType,
        size_t& pixelSize, int& numDimensions);

    // Generic load of an image file.  You may specify a list of required
    // dimensions and/or a list of required pixel types.  If null pointers are
    // passed, no restrictions occur for that parameter.  This is useful for
    // tools that must load an image generically but want restrictions on the
    // type of image loaded.  An ImageD<T> object calls Load with the
    // appropriate requirements for number of dimensions and pixel types.
    bool Load(std::string const& name,
        std::vector<int> const* requiredNumDimensions,
        std::vector<std::string> const* requiredPixelTypes);

    bool Save(std::string const& name) const;

protected:
    // Deallocate arrays and initialize all data members to zero.
    void CreateNullImage();

    // Reallocate the image to match the requested dimensions.  If the image
    // is compatible with the inputs, nothing is modified and the image data
    // remains intact; otherwise, reallocations occur and all image data and
    // metadata information is lost by the operation.  The return value is
    // 'true' whenever the image is resized (the image is not compatible with
    // the inputs).
    bool Resize(int numDimensions, ...);

    // Two images are compatible when they have the same number of dimensions,
    // the same dimensions, and the same pixel type.
    bool IsCompatible(Image const& image) const;

    // Support for copying that avoids reallocation when possible.
    bool Copy(Image const& image);

    // Information needed for typecasting to ImageD<T>, where D is the number
    // of dimensions and T is the pixel type.
    std::string mPixelType;
    size_t mPixelSize;
    int mNumDimensions;
    int* mDimensions;
    size_t* mOffsets;
    size_t mNumPixels;

    // Pixel data.
    char* mRawPixels;
    bool mOwnerRawPixels;

    // Per-image metadata.
    std::string mImageMetaData;

    // Per-pixel metadata.
    std::string* mPixelMetaData;
    std::string mDefaultPixelMetaData;
};


inline std::string Image::GetPixelType() const
{
    return mPixelType;
}

inline size_t Image::GetPixelSize() const
{
    return mPixelSize;
}

inline int Image::GetNumDimensions() const
{
    return mNumDimensions;
}

inline int const* Image::GetDimensions() const
{
    return mDimensions;
}

inline int Image::GetDimension(int d) const
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (0 <= d && d < mNumDimensions)
    {
        return mDimensions[d];
    }
    LogError("Invalid index " + std::to_string(d) + ".");
    return 0;
#else
    return mDimensions[d];
#endif
}

inline size_t const* Image::GetOffsets() const
{
    return mOffsets;
}

inline size_t Image::GetOffset(int d) const
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (0 <= d && d < mNumDimensions)
    {
        return mOffsets[d];
    }
    LogError("Invalid index " + std::to_string(d) + ".");
    return 0;
#else
    return mOffsets[d];
#endif
}

inline size_t Image::GetNumPixels() const
{
    return mNumPixels;
}

inline char const* Image::GetRawPixels() const
{
    return mRawPixels;
}

inline char* Image::GetRawPixels()
{
    return mRawPixels;
}

inline bool Image::IsNullImage() const
{
    // There are many tests for a null image, all equivalent.
    return mRawPixels == nullptr;
}

inline void Image::SetImageMetaData(std::string const& metadata)
{
    mImageMetaData = metadata;
}

inline std::string const& Image::GetImageMetaData() const
{
    return mImageMetaData;
}

inline bool Image::HasPixelMetaData() const
{
    return mPixelMetaData != nullptr;
}

inline void Image::SetPixelMetaData(size_t i, std::string const& metadata)
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (mPixelMetaData && i < mNumPixels)
    {
        mPixelMetaData[i] = metadata;
        return;
    }
    LogError("No metadata or invalid index.");
#else
    mPixelMetaData[i] = metadata;
#endif
}

inline std::string const& Image::GetPixelMetaData(size_t i) const
{
#if defined(GTE_IMAGICS_ASSERT_ON_INVALID_INDEX)
    if (mPixelMetaData && i < mNumPixels)
    {
        return mPixelMetaData[i];
    }
    LogError("No metadata or invalid index.");
    return mDefaultPixelMetaData;
#else
    return mPixelMetaData[i];
#endif
}


}
