#pragma once


// [TODO] move this declaration into a header that g3Types includes!
#ifndef g3External
#ifdef WIN32
#ifdef GEOMETRY3_DLL_EXPORT
#define g3External   __declspec( dllexport )
#else
#define g3External   __declspec( dllimport )
#endif
#else
#define g3External
#endif
#endif



namespace g3
{

g3External void g3_testAssert(bool b);

g3External void g3_devAssert(bool b);

g3External void g3_debugPrint(const char * fmt, ...);
g3External void g3_debugPrint(const wchar_t * fmt, ...);

// must pass-by-copy here because of limitation of varargs va_start macro
g3External void g3_debugPrint(std::string fmt, ...);
g3External void g3_debugPrint(std::wstring fmt, ...);

}


#define gDevAssert g3_devAssert

#define gDebugPrint g3_debugPrint

#define gDevAssertReturnOnFail(x, retval) if (!(x)) { gDevAssert(false); return (retval); }

