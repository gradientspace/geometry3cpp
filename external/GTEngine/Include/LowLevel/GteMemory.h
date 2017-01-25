// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <GTEngineDEF.h>
#include <cstddef>

// The allocation and deallocation functions minimize the number of new and
// delete calls to represent an N-dimensional array.  If 'data' is the
// returned pointer, the actual T objects are stored in contiguous arrays
// pointed to by data[0] for Allocate2, data[0][0] for Allocate3, and
// data[0][0][0] for Allocate4.  The contiguous arrays are dynamically
// allocated, thus owned in a sense by 'data'.  The DeallocateN functions
// will delete all memory allocate in the AllocateN functions.
//
// Sometimes it is convenient to take an already existing contiguous array
// of T objects and wrap it for multidimensional access.  In this case, you
// must provide the contiguous array as input to AllocateMapN and pair that
// call with a DeallocateMapN (if you want to retain ownership of the T
// objects) or DeallocateN (if you want to transfer ownership of the T
// objects).  It is your responsibility to ensure that the 'objects' arrays
// have the correct number of elements.

namespace gte
{

// For 2D arrays:  data[bound1][bound0]
template <typename T>
T** Allocate2(size_t const bound0, size_t const bound1);

template <typename T>
T** AllocateMap2(size_t const bound0, size_t const bound1, T* objects);

// For 3D arrays:  data[bound2][bound1][bound0]
template <typename T>
T*** Allocate3(size_t const bound0, size_t const bound1, size_t const bound2);

template <typename T>
T*** AllocateMap3(size_t const bound0, size_t const bound1,
    size_t const bound2, T* objects);

// For 4D arrays:  data[bound3][bound2][bound1][bound0]
template <typename T>
T**** Allocate4(size_t const bound0, size_t const bound1, size_t const bound2,
    size_t const bound3);

template <typename T>
T**** AllocateMap4(size_t const bound0, size_t const bound1,
    size_t const bound2, size_t const bound3, T* objects);

// For 2D arrays:  data[bound1][bound0]
template <typename T>
void Deallocate2(T**& data);

template <typename T>
void DeallocateMap2(T**& data);

// For 3D arrays:  data[bound2][bound1][bound0]
template <typename T>
void Deallocate3(T***& data);

template <typename T>
void DeallocateMap3(T***& data);

// For 4D arrays:  data[bound3][bound2][bound1][bound0]
template <typename T>
void Deallocate4(T****& data);

template <typename T>
void DeallocateMap4(T****& data);


template <typename T>
T** Allocate2(size_t const bound0, size_t const bound1)
{
    size_t const bound01 = bound0*bound1;
    T** data = new T*[bound1];
    data[0] = new T[bound01];

    for (size_t i1 = 1; i1 < bound1; ++i1)
    {
        size_t j0 = bound0*i1;  // = bound0*(i1 + j1) where j1 = 0
        data[i1] = &data[0][j0];
    }
    return data;
}

template <typename T>
T** AllocateMap2(size_t const bound0, size_t const bound1, T* objects)
{
    T** data = new T*[bound1];
    data[0] = objects;

    for (size_t i1 = 1; i1 < bound1; ++i1)
    {
        size_t j0 = bound0*i1;  // = bound0*(i1 + j1) where j1 = 0
        data[i1] = &data[0][j0];
    }
    return data;
}

template <typename T>
T*** Allocate3(size_t const bound0, size_t const bound1, size_t const bound2)
{
    size_t const bound12 = bound1*bound2;
    size_t const bound012 = bound0*bound12;
    T*** data = new T**[bound2];
    data[0] = new T*[bound12];
    data[0][0] = new T[bound012];

    for (size_t i2 = 0; i2 < bound2; ++i2)
    {
        size_t j1 = bound1*i2;  // = bound1*(i2 + j2) where j2 = 0
        data[i2] = &data[0][j1];
        for (size_t i1 = 0; i1 < bound1; ++i1)
        {
            size_t j0 = bound0*(i1 + j1);
            data[i2][i1] = &data[0][0][j0];
        }
    }
    return data;
}

template <typename T>
T*** AllocateMap3(size_t const bound0, size_t const bound1,
    size_t const bound2, T* objects)
{
    size_t const bound12 = bound1*bound2;
    T*** data = new T**[bound2];
    data[0] = new T*[bound12];
    data[0][0] = objects;

    for (size_t i2 = 0; i2 < bound2; ++i2)
    {
        size_t j1 = bound1*i2;  // = bound1*(i2 + j2) where j2 = 0
        data[i2] = &data[0][j1];
        for (size_t i1 = 0; i1 < bound1; ++i1)
        {
            size_t j0 = bound0*(i1 + j1);
            data[i2][i1] = &data[0][0][j0];
        }
    }
    return data;
}

template <typename T>
T**** Allocate4(size_t const bound0, size_t const bound1, size_t const bound2,
    size_t const bound3)
{
    size_t const bound23 = bound2*bound3;
    size_t const bound123 = bound1*bound23;
    size_t const bound0123 = bound0*bound123;
    T**** data = new T***[bound3];
    data[0] = new T**[bound23];
    data[0][0] = new T*[bound123];
    data[0][0][0] = new T[bound0123];

    for (size_t i3 = 0; i3 < bound3; ++i3)
    {
        size_t j2 = bound2*i3;  // = bound2*(i3 + j3) where j3 = 0
        data[i3] = &data[0][j2];
        for (size_t i2 = 0; i2 < bound2; ++i2)
        {
            size_t j1 = bound1*(i2 + j2);
            data[i3][i2] = &data[0][0][j1];
            for (size_t i1 = 0; i1 < bound1; ++i1)
            {
                size_t j0 = bound0*(i1 + j1);
                data[i3][i2][i1] = &data[0][0][0][j0];
            }
        }
    }
    return data;
}

template <typename T>
T**** AllocateMap4(size_t const bound0, size_t const bound1,
    size_t const bound2, size_t const bound3, T* objects)
{
    size_t const bound23 = bound2*bound3;
    size_t const bound123 = bound1*bound23;
    T**** data = new T***[bound3];
    data[0] = new T**[bound23];
    data[0][0] = new T*[bound123];
    data[0][0][0] = objects;

    for (size_t i3 = 0; i3 < bound3; ++i3)
    {
        size_t j2 = bound2*i3;  // = bound2*(i3 + j3) where j3 = 0
        data[i3] = &data[0][j2];
        for (size_t i2 = 0; i2 < bound2; ++i2)
        {
            size_t j1 = bound1*(i2 + j2);
            data[i3][i2] = &data[0][0][j1];
            for (size_t i1 = 0; i1 < bound1; ++i1)
            {
                size_t j0 = bound0*(i1 + j1);
                data[i3][i2][i1] = &data[0][0][0][j0];
            }
        }
    }
    return data;
}

template <typename T>
void Deallocate2(T**& data)
{
    if (data)
    {
        delete[] data[0];
        delete[] data;
        data = nullptr;
    }
}

template <typename T>
void DeallocateMap2(T**& data)
{
    if (data)
    {
        // data[0] owned by the caller
        delete[] data;
        data = nullptr;
    }
}

template <typename T>
void Deallocate3(T***& data)
{
    if (data)
    {
        delete[] data[0][0];
        delete[] data[0];
        delete[] data;
        data = nullptr;
    }
}

template <typename T>
void DeallocateMap3(T***& data)
{
    if (data)
    {
        // data[0][0] owned by the caller
        delete[] data[0];
        delete[] data;
        data = nullptr;
    }
}

template <typename T>
void Deallocate4(T****& data)
{
    if (data)
    {
        delete[] data[0][0][0];
        delete[] data[0][0];
        delete[] data[0];
        delete[] data;
        data = nullptr;
    }
}

template <typename T>
void DeallocateMap4(T****& data)
{
    if (data)
    {
        // data[0][0][0] owned by the caller
        delete[] data[0][0];
        delete[] data[0];
        delete[] data;
        data = nullptr;
    }
}


}
