// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <GTEngineDEF.h>

namespace gte
{

// Test-intersection queries.

template <typename Real, typename Type0, typename Type1>
class TIQuery
{
public:
    struct Result
    {
        // A TIQuery-base class B must define a B::Result struct with member
        // 'bool intersect'.  A TIQuery-derived class D must also derive a
        // D::Result from B:Result but may have no members.  The member
        // 'intersect' is 'true' iff the primitives intersect.  The operator()
        // is non-const to allow TIQuery to store and modify private state
        // that supports the query.
    };
    Result operator()(Type0 const& primitive0, Type1 const& primitive1);
};

}
