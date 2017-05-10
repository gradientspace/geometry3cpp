// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#include <GTEnginePCH.h>
#include <Mathematics/GteEdgeKey.h>

namespace gte
{


template<>
EdgeKey<true>::EdgeKey(int v0, int v1)
{
    V[0] = v0;
    V[1] = v1;
}

template<>
EdgeKey<false>::EdgeKey(int v0, int v1)
{
    if (v0 < v1)
    {
        // v0 is minimum
        V[0] = v0;
        V[1] = v1;
    }
    else
    {
        // v1 is minimum
        V[0] = v1;
        V[1] = v0;
    }
}


}
