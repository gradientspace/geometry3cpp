// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#if defined(WIN32)
#if !defined(_MSC_VER)
#error Microsoft Visual Studio 2013 or later is required.
#endif
#if _MSC_VER < 1800
#error Microsoft Visual Studio 2013 or later is required.
#endif

// Debug build values (choose_your_value is 0, 1, or 2)
// 0:  Disables checked iterators and disables iterator debugging.
// 1:  Enables checked iterators and disables iterator debugging.
// 2:  (default) Enables iterator debugging; checked iterators are not relevant.
//
// Release build values (choose_your_value is 0 or 1)
// 0:  (default) Disables checked iterators.
// 1:  Enables checked iterators; iterator debugging is not relevant.
//
// #define _ITERATOR_DEBUG_LEVEL choose_your_value

#endif

// TODO: Windows DLL configurations have not yet been added to the project,
// but these defines are required to support them (when we do add them).
//
// Add GTE_EXPORT to project preprocessor options for dynamic library
// configurations to export their symbols.
#if defined(GTE_EXPORT)
    // For the dynamic library configurations.
    #define GTE_IMPEXP __declspec(dllexport)
#else
    // For a client of the dynamic library or for the static library
    // configurations.
    #define GTE_IMPEXP
#endif

// Expose exactly one of these.
#define GTE_USE_ROW_MAJOR
//#define GTE_USE_COL_MAJOR

// Expose exactly one of these.
#define GTE_USE_MAT_VEC
//#define GTE_USE_VEC_MAT

#if (defined(GTE_USE_ROW_MAJOR) && defined(GTE_USE_COL_MAJOR)) || (!defined(GTE_USE_ROW_MAJOR) && !defined(GTE_USE_COL_MAJOR))
#error Exactly one storage order must be specified.
#endif

#if (defined(GTE_USE_MAT_VEC) && defined(GTE_USE_VEC_MAT)) || (!defined(GTE_USE_MAT_VEC) && !defined(GTE_USE_VEC_MAT))
#error Exactly one multiplication convention must be specified.
#endif
