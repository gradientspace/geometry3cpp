#pragma once

#ifndef g3External
#if defined(WIN32) && ! defined(G3_STATIC_LIB)
#ifdef GEOMETRY3_DLL_EXPORT
#define g3External   __declspec( dllexport )
#else
#define g3External   __declspec( dllimport )
#endif
#else
#define g3External
#endif
#endif   // ifndef g3External
