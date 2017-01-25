// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <Imagics/GteImage1.h>

namespace gte
{

class GTE_IMPEXP Histogram
{
public:
    // Construction and destruction.  The first three constructors are for
    // computing the histogram of the entire image at one time.
    ~Histogram ();

    // In the constructor with input 'int const* samples', set noRescaling to
    // 'true' when you want the sample values mapped directly to the buckets.
    // Typically, you know that the sample values are in the set of numbers
    // {0,1,...,numBuckets-1}, but in the event of out-of-range values, the
    // histogram stores a count for those numbers smaller than 0 and those
    // numbers larger or equal to numBuckets.
    Histogram (int numBuckets, int numSamples, int const* samples,
        bool noRescaling);
    Histogram (int numBuckets, int numSamples, float const* samples);
    Histogram (int numBuckets, int numSamples, double const* samples);

    // Construction where you plan on updating the histogram incrementally.
    // The incremental update is implemented only for integer samples and
    // no rescaling.
    Histogram (int numBuckets);

    // This function is called when you have used the Histogram(int)
    // constructor.  No bounds checking is used; you must ensure that the
    // input value is in {0,...,numBuckets-1}.
    inline void Insert (int value);

    // This function is called when you have used the Histogram(int)
    // constructor.  Bounds checking is used.
    void InsertCheck (int value);

    // Member access.
    inline int GetNumBuckets () const;
    inline int const* GetBuckets () const;
    inline int GetExcessLess () const;
    inline int GetExcessGreater () const;

    // In the following, define cdf(V) = sum_{i=0}^{V} bucket[i], where
    // 0 <= V < B and B is the number of buckets.  Define N = cdf(B-1),
    // which must be the number of pixels in the image.

    // Get the lower tail of the histogram.  The returned index L has the
    // properties:  cdf(L-1)/N < tailAmount and cdf(L)/N >= tailAmount.
    int GetLowerTail (double tailAmount);

    // Get the upper tail of the histogram.  The returned index U has the
    // properties:  cdf(U)/N >= 1-tailAmount and cdf(U+1) < 1-tailAmount.
    int GetUpperTail (double tailAmount);

    // Get the lower and upper tails of the histogram.  The returned indices
    // are L and U and have the properties:
    // cdf(L-1)/N < tailAmount/2, cdf(L)/N >= tailAmount/2,
    // cdf(U)/N >= 1-tailAmount/2, and cdf(U+1) < 1-tailAmount/2.
    void GetTails (double tailAmount, int& lower, int& upper);

    // Save the histogram to a text file with two columns.  The first
    // column is the bucket index and the second column is the number of
    // items in the bucket.
    void SaveAsText (std::string const& name);

    // Save the histogram to an integer-valued image with dimension0 columns
    // and numBuckets rows.  The image is 32-bit, but only black and white
    // colors are used.
    void SaveAsImage (std::string const& name, int dimension0);

private:
    int mNumBuckets;
    int* mBuckets;
    int mExcessLess, mExcessGreater;
};


inline void Histogram::Insert(int value)
{
    ++mBuckets[value];
}

inline int Histogram::GetNumBuckets() const
{
    return mNumBuckets;
}

inline int const* Histogram::GetBuckets() const
{
    return mBuckets;
}

inline int Histogram::GetExcessLess() const
{
    return mExcessLess;
}

inline int Histogram::GetExcessGreater() const
{
    return mExcessGreater;
}


}
