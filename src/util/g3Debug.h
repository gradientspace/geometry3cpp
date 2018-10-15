#pragma once

#include <g3Config.h>


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

// TODO: https://stackoverflow.com/questions/173618/is-there-a-portable-equivalent-to-debugbreak-debugbreak
#define gBreakToDebugger __debugbreak