// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <Mathematics/GteFeatureKey.h>

namespace gte
{

template <bool Ordered>
class EdgeKey : public FeatureKey<2, Ordered>
{
public:
    // An ordered edge has (V[0],V[1]) = (v0,v1).  An unordered edge has
    // (V[0],V[1]) = (min(V[0],V[1]),max(V[0],V[1])).
    EdgeKey(int v0 = -1, int v1 = -1);
};

}
