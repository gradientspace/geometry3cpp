// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#include <GTEnginePCH.h>
#include <Imagics/GteHistogram.h>
#include <Imagics/GteImage2.h>
#include <algorithm>
using namespace gte;

Histogram::Histogram(int numBuckets, int numSamples, int const* samples, bool noRescaling)
    :
    mBuckets(numBuckets),
    mExcessLess(0),
    mExcessGreater(0)
{
    if (numBuckets <= 0 || numSamples <= 0 || !samples)
    {
        LogError("Invalid inputs.");
        return;
    }

    std::fill(mBuckets.begin(), mBuckets.end(), 0);

    int i, value;

    if (noRescaling)
    {
        // Map to the buckets, also counting out-of-range pixels.
        for (i = 0; i < numSamples; ++i)
        {
            value = samples[i];
            if (0 <= value)
            {
                if (value < numBuckets)
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
            double numer = static_cast<double>(numBuckets - 1);
            double denom = static_cast<double>(maxValue - minValue);
            double mult = numer / denom;
            for (i = 0; i < numSamples; ++i)
            {
                int index = static_cast<int>(mult * static_cast<double>(samples[i] - minValue));
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

Histogram::Histogram(int numBuckets, int numSamples, float const* samples)
    :
    mBuckets(numBuckets),
    mExcessLess(0),
    mExcessGreater(0)
{
    if (numBuckets <= 0 || numSamples <= 0 || !samples)
    {
        LogError("Invalid inputs.");
        return;
    }

    std::fill(mBuckets.begin(), mBuckets.end(), 0);

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
        double numer = static_cast<double>(numBuckets - 1);
        double denom = static_cast<double>(maxValue - minValue);
        double mult = numer / denom;
        for (int i = 0; i < numSamples; ++i)
        {
            int index = static_cast<int>(mult * static_cast<double>(samples[i] - minValue));
            ++mBuckets[index];
        }
    }
    else
    {
        // The image is constant.
        mBuckets[0] = numSamples;
    }
}

Histogram::Histogram(int numBuckets, int numSamples, double const* samples)
    :
    mBuckets(numBuckets),
    mExcessLess(0),
    mExcessGreater(0)
{
    if (numBuckets <= 0 || numSamples <= 0 || !samples)
    {
        LogError("Invalid inputs.");
        return;
    }

    std::fill(mBuckets.begin(), mBuckets.end(), 0);

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
        double numer = static_cast<double>(numBuckets - 1);
        double denom = maxValue - minValue;
        double mult = numer/denom;
        for (int i = 0; i < numSamples; ++i)
        {
            int index = static_cast<int>(mult * (samples[i] - minValue));
            ++mBuckets[index];
        }
    }
    else
    {
        // The image is constant.
        mBuckets[0] = numSamples;
    }
}

Histogram::Histogram(int numBuckets)
    :
    mBuckets(numBuckets),
    mExcessLess(0),
    mExcessGreater(0)
{
    if (numBuckets <= 0)
    {
        LogError("Invalid inputs.");
        return;
    }

    std::fill(mBuckets.begin(), mBuckets.end(), 0);
}

void Histogram::InsertCheck(int value)
{
    if (0 <= value)
    {
        if (value < static_cast<int>(mBuckets.size()))
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

int Histogram::GetLowerTail(double tailAmount)
{
    int const numBuckets = static_cast<int>(mBuckets.size());
    int hSum = 0;
    for (int i = 0; i < numBuckets; ++i)
    {
        hSum += mBuckets[i];
    }

    int hTailSum = static_cast<int>(tailAmount * hSum);
    int hLowerSum = 0;
    int lower;
    for (lower = 0; lower < numBuckets; ++lower)
    {
        hLowerSum += mBuckets[lower];
        if (hLowerSum >= hTailSum)
        {
            break;
        }
    }
    return lower;
}

int Histogram::GetUpperTail(double tailAmount)
{
    int const numBuckets = static_cast<int>(mBuckets.size());
    int hSum = 0;
    for (int i = 0; i < numBuckets; ++i)
    {
        hSum += mBuckets[i];
    }

    int hTailSum = static_cast<int>(tailAmount * hSum);
    int hUpperSum = 0;
    int upper;
    for (upper = numBuckets - 1; upper >= 0; --upper)
    {
        hUpperSum += mBuckets[upper];
        if (hUpperSum >= hTailSum)
        {
            break;
        }
    }
    return upper;
}

void Histogram::GetTails(double tailAmount, int& lower, int& upper)
{
    int const numBuckets = static_cast<int>(mBuckets.size());
    int hSum = 0;
    for (int i = 0; i < numBuckets; ++i)
    {
        hSum += mBuckets[i];
    }

    int hTailSum = static_cast<int>(0.5 * tailAmount * hSum);
    int hLowerSum = 0;
    for (lower = 0; lower < numBuckets; ++lower)
    {
        hLowerSum += mBuckets[lower];
        if (hLowerSum >= hTailSum)
        {
            break;
        }
    }

    int hUpperSum = 0;
    for (upper = numBuckets - 1; upper >= 0; --upper)
    {
        hUpperSum += mBuckets[upper];
        if (hUpperSum >= hTailSum)
        {
            break;
        }
    }
}
