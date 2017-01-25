// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.2 (2015/11/21)

#ifndef WM5DX9RENDERERLIB_H
#define WM5DX9RENDERERLIB_H

#include "Wm5GraphicsLIB.h"

// Disable the 'min' and 'max' macros that are sucked in by the indirect
// inclusion of windows.h in the DirectX header files.  These conflict
// with std::numeric_limits<type>::max().
#ifndef NOMINMAX
#define NOMINMAX
#endif

#include <d3d9.h>
#include <d3dx9.h>

// NOTE:  When trapping compiler and linker problems with MSVS 2015,
// the DebugDLL|x64 configuration generated a linker error for dxerr.lib,
// complaining about a _vsnprintf reference in that library.  The function
// appears to be deprecated.  Removing the hr-to-string diagnostics.
//
//#include <dxerr.h>
#define DXGetErrorString(hr) "<unknown>"

#endif
