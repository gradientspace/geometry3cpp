// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#include <GTEnginePCH.h>
#include <Mathematics/GteTriangleKey.h>
#include <algorithm>

namespace gte
{


template<>
TriangleKey<true>::TriangleKey(int v0, int v1, int v2)
{
    if (v0 < v1)
    {
        if (v0 < v2)
        {
            // v0 is minimum
            V[0] = v0;
            V[1] = v1;
            V[2] = v2;
        }
        else
        {
            // v2 is minimum
            V[0] = v2;
            V[1] = v0;
            V[2] = v1;
        }
    }
    else
    {
        if (v1 < v2)
        {
            // v1 is minimum
            V[0] = v1;
            V[1] = v2;
            V[2] = v0;
        }
        else
        {
            // v2 is minimum
            V[0] = v2;
            V[1] = v0;
            V[2] = v1;
        }
    }
}

template<>
TriangleKey<false>::TriangleKey(int v0, int v1, int v2)
{
    if (v0 < v1)
    {
        if (v0 < v2)
        {
            // v0 is minimum
            V[0] = v0;
            V[1] = std::min(v1, v2);
            V[2] = std::max(v1, v2);
        }
        else
        {
            // v2 is minimum
            V[0] = v2;
            V[1] = std::min(v0, v1);
            V[2] = std::max(v0, v1);
        }
    }
    else
    {
        if (v1 < v2)
        {
            // v1 is minimum
            V[0] = v1;
            V[1] = std::min(v2, v0);
            V[2] = std::max(v2, v0);
        }
        else
        {
            // v2 is minimum
            V[0] = v2;
            V[1] = std::min(v0, v1);
            V[2] = std::max(v0, v1);
        }
    }
}


}
