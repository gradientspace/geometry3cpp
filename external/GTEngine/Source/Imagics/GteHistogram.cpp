// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#include <GTEnginePCH.h>
#include <Imagics/GteHistogram.h>
#include <Imagics/GteImage2.h>
using namespace gte;


Histogram::~Histogram ()
{
    delete[] mBuckets;
}

Histogram::Histogram (int numBuckets, int numSamples, int const* samples,
    bool noRescaling)
    :
    mNumBuckets(numBuckets),
    mBuckets(nullptr),
    mExcessLess(0),
    mExcessGreater(0)
{
    if (mNumBuckets <= 0 || numSamples <= 0 || !samples)
    {
        LogError("Invalid inputs.");
        mNumBuckets = 0;
        return;
    }

    mBuckets = new int[mNumBuckets];
    memset(mBuckets, 0, mNumBuckets*sizeof(int));

    int i, value;

    if (noRescaling)
    {
        // Map to the buckets, also counting out-of-range pixels.
        for (i = 0; i < numSamples; ++i)
        {
            value = samples[i];
            if (0 <= value)
            {
                if (value < mNumBuckets)
                {
                    ++mBuckets[value];
                }
                else
                {
                    ++mExcessGreater;
                }
            }
            else
            {
                ++mExcessLess;
            }
        }
    }
    else
    {
        // Compute the extremes.
        int minValue = samples[0], maxValue = minValue;
        for (i = 1; i < numSamples; ++i)
        {
            value = samples[i];
            if (value < minValue)
            {
                minValue = value;
            }
            else if (value > maxValue)
            {
                maxValue = value;
            }
        }

        // Map to the buckets.
        if (minValue < maxValue)
        {
            // The image is not constant.
            double numer = (double)(mNumBuckets - 1);
            double denom = (double)(maxValue - minValue);
            double mult = numer/denom;
            for (i = 0; i < numSamples; ++i)
            {
                int index = (int)(mult*(double)(samples[i] - minValue));
                ++mBuckets[index];
            }
        }
        else
        {
            // The image is constant.
            mBuckets[0] = numSamples;
        }
    }
}

Histogram::Histogram (int numBuckets, int numSamples, float const* samples)
    :
    mNumBuckets(numBuckets),
    mBuckets(nullptr),
    mExcessLess(0),
    mExcessGreater(0)
{
    if (mNumBuckets <= 0 || numSamples <= 0 || !samples)
    {
        LogError("Invalid inputs.");
        mNumBuckets = 0;
        return;
    }

    mBuckets = new int[mNumBuckets];
    memset(mBuckets, 0, mNumBuckets*sizeof(int));

    // Compute the extremes.
    float minValue = samples[0], maxValue = minValue;
    for (int i = 1; i < numSamples; ++i)
    {
        float value = samples[i];
        if (value < minValue)
        {
            minValue = value;
        }
        else if (value > maxValue)
        {
            maxValue = value;
        }
    }

    // Map to the buckets.
    if (minValue < maxValue)
    {
        // The image is not constant.
        double numer = (double)(mNumBuckets - 1);
        double denom = (double)(maxValue - minValue);
        double mult = numer/denom;
        for (int i = 0; i < numSamples; ++i)
        {
            int index = (int)(mult*(double)(samples[i] - minValue));
            ++mBuckets[index];
        }
    }
    else
    {
        // The image is constant.
        mBuckets[0] = numSamples;
    }
}

Histogram::Histogram (int numBuckets, int numSamples, double const* samples)
    :
    mNumBuckets(numBuckets),
    mBuckets(nullptr),
    mExcessLess(0),
    mExcessGreater(0)
{
    if (mNumBuckets <= 0 || numSamples <= 0 || !samples)
    {
        LogError("Invalid inputs.");
        mNumBuckets = 0;
        return;
    }

    mBuckets = new int[mNumBuckets];
    memset(mBuckets, 0, mNumBuckets*sizeof(int));

    // Compute the extremes.
    double minValue = samples[0], maxValue = minValue;
    for (int i = 1; i < numSamples; ++i)
    {
        double value = samples[i];
        if (value < minValue)
        {
            minValue = value;
        }
        else if (value > maxValue)
        {
            maxValue = value;
        }
    }

    // Map to the buckets.
    if (minValue < maxValue)
    {
        // The image is not constant.
        double numer = (double)(mNumBuckets - 1);
        double denom = maxValue - minValue;
        double mult = numer/denom;
        for (int i = 0; i < numSamples; ++i)
        {
            int index = (int)(mult*(samples[i] - minValue));
            ++mBuckets[index];
        }
    }
    else
    {
        // The image is constant.
        mBuckets[0] = numSamples;
    }
}

Histogram::Histogram (int numBuckets)
    :
    mNumBuckets(numBuckets),
    mBuckets(nullptr),
    mExcessLess(0),
    mExcessGreater(0)
{
    if (mNumBuckets <= 0)
    {
        LogError("Invalid inputs.");
        mNumBuckets = 0;
        return;
    }

    mBuckets = new int[mNumBuckets];
    memset(mBuckets, 0, mNumBuckets*sizeof(int));
}

void Histogram::InsertCheck (int value)
{
    if (0 <= value)
    {
        if (value < mNumBuckets)
        {
            ++mBuckets[value];
        }
        else
        {
            ++mExcessGreater;
        }
    }
    else
    {
        ++mExcessLess;
    }
}

int Histogram::GetLowerTail (double tailAmount)
{
    int hSum = 0;
    for (int i = 0; i < mNumBuckets; ++i)
    {
        hSum += mBuckets[i];
    }

    int hTailSum = (int)(tailAmount*hSum);

    int hLowerSum = 0;
    int lower;
    for (lower = 0; lower < mNumBuckets; ++lower)
    {
        hLowerSum += mBuckets[lower];
        if (hLowerSum >= hTailSum)
        {
            break;
        }
    }
    return lower;
}

int Histogram::GetUpperTail (double tailAmount)
{
    int hSum = 0;
    for (int i = 0; i < mNumBuckets; ++i)
    {
        hSum += mBuckets[i];
    }

    int hTailSum = (int)(tailAmount*hSum);

    int hUpperSum = 0;
    int upper;
    for (upper = mNumBuckets - 1; upper >= 0; --upper)
    {
        hUpperSum += mBuckets[upper];
        if (hUpperSum >= hTailSum)
        {
            break;
        }
    }
    return upper;
}

void Histogram::GetTails (double tailAmount, int& lower, int& upper)
{
    int hSum = 0;
    for (int i = 0; i < mNumBuckets; ++i)
    {
        hSum += mBuckets[i];
    }

    int hTailSum = (int)(0.5*tailAmount*hSum);

    int hLowerSum = 0;
    for (lower = 0; lower < mNumBuckets; ++lower)
    {
        hLowerSum += mBuckets[lower];
        if (hLowerSum >= hTailSum)
        {
            break;
        }
    }

    int hUpperSum = 0;
    for (upper = mNumBuckets - 1; upper >= 0; --upper)
    {
        hUpperSum += mBuckets[upper];
        if (hUpperSum >= hTailSum)
        {
            break;
        }
    }
}

void Histogram::SaveAsText (std::string const& name)
{
    std::ofstream output(name);
    if (!output)
    {
        LogError("Failed to open file " + name + ".");
        return;
    }

    for (int i = 0; i < mNumBuckets; ++i)
    {
        output << i << ' ' << mBuckets[i] << std::endl;
    }

    output.close();
}

void Histogram::SaveAsImage (std::string const& name, int dimension0)
{
    Image2<int> image(dimension0, mNumBuckets);

    // Get the maximum bucket value.  This is used to scale the histogram to
    // fit within the image bounds.
    int hMax = mBuckets[0];
    int y;
    for (y = 1; y < mNumBuckets; ++y)
    {
        int hValue = mBuckets[y];
        if (hValue > hMax)
        {
            hMax = hValue;
        }
    }

    // Map the histogram to the image.
    double mult = ((double)dimension0 - 1.0)/(double)hMax;
    for (y = 0; y < mNumBuckets; ++y)
    {
        int x = (int)(mult*mBuckets[y]);
        while (x >= 0)
        {
            image(x, y) = 1;
            --x;
        }
    }

    image.Save(name);
}

