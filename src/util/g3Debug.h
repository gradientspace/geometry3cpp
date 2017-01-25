#pragma once

namespace g3
{

void g3_testAssert(bool b);

void g3_devAssert(bool b);

void g3_debugPrint(const char * fmt, ...);
void g3_debugPrint(const wchar_t * fmt, ...);

// must pass-by-copy here because of limitation of varargs va_start macro
void g3_debugPrint(std::string fmt, ...);
void g3_debugPrint(std::wstring fmt, ...);

}


#define gDevAssert g3_devAssert

#define gDebugPrint g3_debugPrint

#define gDevAssertReturnOnFail(x, retval) if (!(x)) { gDevAssert(false); return (retval); }

