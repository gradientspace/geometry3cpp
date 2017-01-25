// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#include <GTEnginePCH.h>
#include <LowLevel/GteLogger.h>
#include <LowLevel/GteWrapper.h>
#include <Imagics/GteImage.h>
#include <cstdarg>
using namespace gte;


Image::~Image()
{
    delete[] mDimensions;
    delete[] mOffsets;
    delete[] mPixelMetaData;

    if (mOwnerRawPixels)
    {
        delete[] mRawPixels;
    }
}

Image::Image()
    :
    mPixelType(""),
    mPixelSize(0),
    mNumDimensions(0),
    mDimensions(nullptr),
    mOffsets(nullptr),
    mNumPixels(0),
    mRawPixels(nullptr),
    mOwnerRawPixels(true),
    mImageMetaData(""),
    mPixelMetaData(nullptr),
    mDefaultPixelMetaData("")
{
}

Image::Image(Image const& image)
    :
    mPixelType(image.mPixelType),
    mPixelSize(image.mPixelSize),
    mNumDimensions(image.mNumDimensions),
    mNumPixels(image.mNumPixels),
    mOwnerRawPixels(true),
    mImageMetaData(image.mImageMetaData),
    mDefaultPixelMetaData(image.mDefaultPixelMetaData)
{
    mDimensions = new int[mNumDimensions];
    mOffsets = new size_t[mNumDimensions];
    for (int d = 0; d < mNumDimensions; ++d)
    {
        mDimensions[d] = image.mDimensions[d];
        mOffsets[d] = image.mOffsets[d];
    }

    size_t numBytes = mNumPixels * mPixelSize;
    mRawPixels = new char[numBytes];
    Memcpy(mRawPixels, image.mRawPixels, numBytes);

    if (image.mPixelMetaData)
    {
        mPixelMetaData = new std::string[mNumPixels];
        for (size_t i = 0; i < mNumPixels; ++i)
        {
            mPixelMetaData[i] = image.mPixelMetaData[i];
        }
    }
    else
    {
        mPixelMetaData = nullptr;
    }
}

Image::Image(std::string const& pixelType, size_t pixelSize,
    int numDimensions, ...)
    :
    mPixelType(pixelType),
    mPixelSize(pixelSize),
    mNumDimensions(numDimensions),
    mDimensions(nullptr),
    mOffsets(nullptr),
    mNumPixels(0),
    mRawPixels(nullptr),
    mOwnerRawPixels(true),
    mImageMetaData(""),
    mPixelMetaData(nullptr),
    mDefaultPixelMetaData("")
{
    if (mNumDimensions > 0)
    {
        mDimensions = new int[mNumDimensions];
        mOffsets = new size_t[mNumDimensions];

        va_list arguments;
        va_start(arguments, numDimensions);
        mNumPixels = 1;
        for (int d = 0; d < mNumDimensions; ++d)
        {
            mDimensions[d] = va_arg(arguments, int);
            mNumPixels *= mDimensions[d];
        }
        va_end(arguments);

        mOffsets[0] = 1;
        for (int d = 1; d < mNumDimensions; ++d)
        {
            mOffsets[d] = mDimensions[d - 1]*mOffsets[d - 1];
        }

        mRawPixels = new char[mNumPixels*mPixelSize];
    }
}

Image& Image::operator=(Image const& image)
{
    Copy(image);
    return *this;
}

void Image::SetRawPixels(char* rawPixels)
{
    if (mOwnerRawPixels)
    {
        delete[] mRawPixels;
    }

    mRawPixels = rawPixels;
    mOwnerRawPixels = false;
}

void Image::CreatePixelMetaData()
{
    if (!mPixelMetaData && mNumPixels > 0)
    {
        mPixelMetaData = new std::string[mNumPixels];
        return;
    }
    LogError("No metadata.");
}

void Image::DestroyPixelMetaData()
{
    if (mPixelMetaData)
    {
        delete[] mPixelMetaData;
        mPixelMetaData = nullptr;
        return;
    }
    LogError("No metadata.");
}

void Image::ClearPixelMetaData()
{
    if (mPixelMetaData)
    {
        for (size_t i = 0; i < mNumPixels; ++i)
        {
            mPixelMetaData[i] = "";
        }
    }
}

size_t Image::GetIndex(int const* coord) const
{
    // assert:  coord is array of mNumDimensions elements
    size_t index = coord[0];
    for (int d = 1; d < mNumDimensions; ++d)
    {
        index += mOffsets[d]*coord[d];
    }
    return index;
}

void Image::GetCoordinates(size_t index, int* coord) const
{
    // assert:  coord is array of mNumDimensions elements
    for (int d = 0; d < mNumDimensions; ++d)
    {
        coord[d] = index % mDimensions[d];
        index /= mDimensions[d];
    }
}

void Image::ClearPixels()
{
    if (mRawPixels)
    {
        memset(mRawPixels, 0, mNumPixels*mPixelSize);
    }

    if (mPixelMetaData)
    {
        for (size_t i = 0; i < mNumPixels; ++i)
        {
            mPixelMetaData[i] = "";
        }
    }
}

bool Image::LoadHeader(std::ifstream& input, std::string& pixelType,
    size_t& pixelSize, int& numDimensions)
{
    // Read the length of the pixel type string.
    int length;
    if (input.read((char*)&length, sizeof(int)).bad())
    {
        LogError("Failed read length.");
        pixelType = "";
        pixelSize = 0;
        numDimensions = 0;
        return false;
    }

    // Read the pixel type string.
    char* temp = new char[length + 1];
    if (input.read(temp, length + 1).bad() || temp[length] != 0)
    {
        LogError("Failed read pixelType.");
        delete[] temp;
        pixelType = "";
        pixelSize = 0;
        numDimensions = 0;
        return false;
    }
    pixelType = std::string(temp);
    delete[] temp;

    // Read the pixel size.
    if (input.read((char*)&pixelSize, sizeof(size_t)).bad())
    {
        LogError("Failed read pixelSize.");
        pixelType = "";
        pixelSize = 0;
        numDimensions = 0;
        return false;
    }

    // Read the number of dimensions.
    if (input.read((char*)&numDimensions, sizeof(int)).bad())
    {
        LogError("Failed read numDimensions.");
        pixelType = "";
        pixelSize = 0;
        numDimensions = 0;
        return false;
    }

    return true;
}

bool Image::Load(std::string const& name,
    std::vector<int> const* requiredNumDimensions,
    std::vector<std::string> const* requiredPixelTypes)
{
    CreateNullImage();

    std::ifstream input(name, std::ios::in | std::ios::binary);
    if (!input)
    {
        LogError("Failed to open file " + name + " for reading.");
        return false;
    }

    if (!LoadHeader(input, mPixelType, mPixelSize, mNumDimensions))
    {
        input.close();
        return false;
    }

    if (requiredNumDimensions)
    {
        int const numElements = (int)requiredNumDimensions->size();
        int j;
        for (j = 0; j < numElements; ++j)
        {
            if (mNumDimensions == (*requiredNumDimensions)[j])
            {
                break;
            }
        }
        if (j == numElements)
        {
            CreateNullImage();
            input.close();
            return false;
        }
    }

    if (requiredPixelTypes)
    {
        int const numElements = (int)requiredPixelTypes->size();
        int j;
        for (j = 0; j < numElements; ++j)
        {
            if (mPixelType == (*requiredPixelTypes)[j])
            {
                break;
            }
        }
        if (j == numElements)
        {
            CreateNullImage();
            input.close();
            return false;
        }
    }

    // Read the dimensions.
    mDimensions = new int[mNumDimensions];
    if (input.read((char*)mDimensions, mNumDimensions*sizeof(int)).bad())
    {
        LogError("Failed read dimensions.");
        CreateNullImage();
        input.close();
        return false;
    }

    // Read the offsets.
    mOffsets = new size_t[mNumDimensions];
    if (input.read((char*)mOffsets, mNumDimensions*sizeof(size_t)).bad())
    {
        LogError("Failed read offsets.");
        CreateNullImage();
        input.close();
        return false;
    }

    // Read the number of pixels.
    if (input.read((char*)&mNumPixels, sizeof(size_t)).bad())
    {
        LogError("Failed read number of pixels.");
        CreateNullImage();
        input.close();
        return false;
    }

    // Read the pixels.
    size_t numBytes = mNumPixels*mPixelSize;
    mRawPixels = new char[numBytes];
    if (input.read(mRawPixels, numBytes).bad())
    {
        LogError("Failed read pixels.");
        CreateNullImage();
        input.close();
        return false;
    }

    // Read the length of the image metadata string.
    int length;
    if (input.read((char*)&length, sizeof(int)).bad())
    {
        LogError("Failed read length(image metadata).");
        CreateNullImage();
        input.close();
        return false;
    }

    // Read the image metadata string.
    numBytes = length + 1;
    char* temp = new char[numBytes];
    if (input.read(temp, numBytes).bad() || temp[length] != 0)
    {
        LogError("Failed read image metadata.");
        delete[] temp;
        CreateNullImage();
        input.close();
        return false;
    }
    mImageMetaData = std::string(temp);
    delete[] temp;

    // Read the existence flag for pixel metadata.
    int existsPixelMetaData = 0;
    if (input.read((char*)&existsPixelMetaData, sizeof(int)).bad())
    {
        LogError("Failed read pixel metadata existence.");
        CreateNullImage();
        input.close();
        return false;
    }

    if (existsPixelMetaData)
    {
        CreatePixelMetaData();

        size_t maxNumBytes = 128;
        temp = new char[maxNumBytes];
        for (size_t i = 0; i < mNumPixels; ++i)
        {
            // Read the length of the pixel metadata string.
            if (input.read((char*)&length, sizeof(int)).bad())
            {
                LogError("Failed read length.");
                delete[] temp;
                CreateNullImage();
                input.close();
                return false;
            }

            // Resize the temporary storage if necessary.
            numBytes = length + 1;
            if (numBytes > maxNumBytes)
            {
                maxNumBytes = numBytes;
                delete[] temp;
                temp = new char[maxNumBytes];
            }

            // Read the pixel metadata string.
            if (input.read(temp, numBytes).bad() || temp[length] != 0)
            {
                LogError("Failed read pixel metadata.");
                delete[] temp;
                CreateNullImage();
                input.close();
                return false;
            }
            mPixelMetaData[i] = std::string(temp);
        }
    }

    input.close();
    return true;
}

bool Image::Save(std::string const& name) const
{
    if (IsNullImage())
    {
        LogError("Cannot save a null image.");
        return false;
    }

    std::ofstream output(name, std::ios::out | std::ios::binary);
    if (!output)
    {
        LogError("Failed to open file " + name + " for writing.");
        return false;
    }

    // Write the length of the pixel type.
    int length = (int)mPixelType.length();
    if (output.write((char const*)&length, sizeof(int)).bad())
    {
        LogError("Failed write length(pixelType).");
        output.close();
        return false;
    }

    // Write the pixel type string.
    size_t numBytes = length + 1;
    if (output.write(mPixelType.c_str(), numBytes).bad())
    {
        LogError("Failed write pixelType.");
        output.close();
        return false;
    }

    // Write the pixel size.
    if (output.write((char const*)&mPixelSize, sizeof(size_t)).bad())
    {
        LogError("Failed write pixelSize.");
        output.close();
        return false;
    }

    // Write the number of dimensions.
    if (output.write((char const*)&mNumDimensions, sizeof(int)).bad())
    {
        LogError("Failed write numDimensions.");
        output.close();
        return false;
    }

    // Write the dimensions.
    numBytes = mNumDimensions*sizeof(int);
    if (output.write((char const*)mDimensions, numBytes).bad())
    {
        LogError("Failed write dimensions.");
        output.close();
        return false;
    }

    // Write the offsets.
    numBytes = mNumDimensions*sizeof(size_t);
    if (output.write((char const*)mOffsets, numBytes).bad())
    {
        LogError("Failed write offsets.");
        output.close();
        return false;
    }

    // Write the number of pixels.
    if (output.write((char const*)&mNumPixels, sizeof(size_t)).bad())
    {
        LogError("Failed write numPixels.");
        output.close();
        return false;
    }

    // Write the pixels.
    numBytes = mNumPixels*mPixelSize;
    if (output.write(mRawPixels, numBytes).bad())
    {
        LogError("Failed write pixels.");
        output.close();
        return false;
    }

    // Write the length of the image metadata string.
    length = (int)mImageMetaData.length();
    if (output.write((char const*)&length, sizeof(int)).bad())
    {
        LogError("Failed write length(image metadata).");
        output.close();
        return false;
    }

    // Write the image metadata string.
    numBytes = length + 1;
    if (output.write(mImageMetaData.c_str(), numBytes).bad())
    {
        LogError("Failed write image metadata.");
        output.close();
        return false;
    }

    // Write the existence flag for pixel metadata.
    int existsPixelMetaData = (mPixelMetaData ? 1 : 0);
    if (output.write((char const*)&existsPixelMetaData, sizeof(int)).bad())
    {
        LogError("Failed write pixel metadata existence.");
        output.close();
        return false;
    }

    if (mPixelMetaData)
    {
        for (size_t i = 0; i < mNumPixels; ++i)
        {
            // Write the length of the pixel metadata string.
            std::string const& metadata = mPixelMetaData[i];
            length = (int)metadata.length();
            if (output.write((char const*)&length, sizeof(int)).bad())
            {
                LogError("Failed write length(pixel metadata).");
                output.close();
                return false;
            }

            // Write the pixel metadata string.
            numBytes = length + 1;
            if (output.write(metadata.c_str(), numBytes).bad())
            {
                LogError("Failed write pixel metadata.");
                output.close();
                return false;
            }
        }
    }

    output.close();
    return true;
}

void Image::CreateNullImage()
{
    delete[] mDimensions;
    delete[] mOffsets;
    delete[] mPixelMetaData;
    if (mOwnerRawPixels)
    {
        delete[] mRawPixels;
    }

    mPixelType = "";
    mPixelSize = 0;
    mNumDimensions = 0;
    mDimensions = nullptr;
    mOffsets = nullptr;
    mNumPixels = 0;
    mRawPixels = nullptr;
    mOwnerRawPixels = true;
    mImageMetaData = "";
    mPixelMetaData = nullptr;
    mDefaultPixelMetaData = "";
}

bool Image::Resize(int numDimensions, ...)
{
    if (numDimensions < 0)
    {
        LogError("Invalid number of dimensions.");
        numDimensions = 0;
    }
    else if (numDimensions > 0 && numDimensions == mNumDimensions)
    {
        // Test for compatibility.  If they are, no resizing is necessary.
        va_list arguments;
        va_start(arguments, numDimensions);
        int d;
        for (d = 0; d < mNumDimensions; ++d)
        {
            int inputDimension = va_arg(arguments, int);
            if (mDimensions[d] != inputDimension)
            {
                break;
            }
        }
        va_end(arguments);
        if (d == mNumDimensions)
        {
            // The images are compatible.
            return false;
        }
    }

    std::string savePixelType = mPixelType;
    size_t savePixelSize = mPixelSize;
    bool existsPixelMetaData = (mPixelMetaData != nullptr);
    CreateNullImage();
    mPixelType = savePixelType;
    mPixelSize = savePixelSize;
    mNumDimensions = numDimensions;

    if (mNumDimensions > 0)
    {
        mDimensions = new int[mNumDimensions];
        mOffsets = new size_t[mNumDimensions];

        va_list arguments;
        va_start(arguments, numDimensions);
        mNumPixels = 1;
        int d;
        for (d = 0; d < mNumDimensions; ++d)
        {
            mDimensions[d] = va_arg(arguments, int);
            mNumPixels *= mDimensions[d];
        }
        va_end(arguments);

        mOffsets[0] = 1;
        for (d = 1; d < mNumDimensions; ++d)
        {
            mOffsets[d] = mDimensions[d - 1]*mOffsets[d - 1];
        }

        mRawPixels = new char[mNumPixels*mPixelSize];

        if (existsPixelMetaData)
        {
            mPixelMetaData = new std::string[mNumPixels];
        }
    }

    return true;
}

bool Image::IsCompatible(Image const& image) const
{
    if (mNumDimensions != image.mNumDimensions
        || mPixelType != image.mPixelType)
    {
        return false;
    }

    for (int d = 0; d < mNumDimensions; ++d)
    {
        if (mDimensions[d] != image.mDimensions[d])
        {
            return false;
        }
    }

    return true;
}

bool Image::Copy(Image const& image)
{
    bool compatible = IsCompatible(image);
    if (!compatible)
    {
        CreateNullImage();

        mPixelType = image.mPixelType;
        mPixelSize = image.mPixelSize;
        mNumDimensions = image.mNumDimensions;
        mNumPixels = image.mNumPixels;
        mImageMetaData = image.mImageMetaData;
        mDefaultPixelMetaData = image.mDefaultPixelMetaData;

        if (mNumDimensions > 0)
        {
            mDimensions = new int[mNumDimensions];
            mOffsets = new size_t[mNumDimensions];
            for (int d = 0; d < mNumDimensions; ++d)
            {
                mDimensions[d] = image.mDimensions[d];
                mOffsets[d] = image.mOffsets[d];
            }

            mRawPixels = new char[mNumPixels*mPixelSize];

            if (image.mPixelMetaData)
            {
                mPixelMetaData = new std::string[mNumPixels];
            }
        }
        else
        {
            mDimensions = nullptr;
            mOffsets = nullptr;
            mRawPixels = nullptr;
            mPixelMetaData = nullptr;
        }
    }

    if (mRawPixels)
    {
        Memcpy(mRawPixels, image.mRawPixels, mNumPixels*mPixelSize);
    }

    if (mPixelMetaData && image.mPixelMetaData)
    {
        for (size_t i = 0; i < mNumPixels; ++i)
        {
            mPixelMetaData[i] = image.mPixelMetaData[i];
        }
    }

    return compatible;
}

