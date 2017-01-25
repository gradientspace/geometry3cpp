// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <GTEngineDEF.h>
#include <cstddef>

// Wrappers around platform-specific low-level library calls.

namespace gte
{

void Memcpy(void* target, void const* source, size_t count);
void Memcpy(wchar_t* target, wchar_t const* source, size_t count);

}
