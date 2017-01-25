#include <geometry3PCH.h>
#include "g3Debug.h"

#include <stdarg.h>
#include <stdio.h>
#include <wchar.h>

// for OutputDebugString
#ifdef WIN32
#include <windows.h>
#endif

using namespace g3;

void g3::g3_testAssert(bool b)
{
    if (!b)
#ifdef __APPLE__
        __asm__("int $3\n" : : );
#else
    __debugbreak();
#endif
}

void g3::g3_devAssert(bool b)
{
	if (!b)
#ifdef __APPLE__
        __asm__("int $3\n" : : );        
#else
		__debugbreak();
#endif
}

void g3::g3_debugPrint(std::string fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	g3_debugPrint(fmt.c_str(), args);
	va_end(args);
}
void g3::g3_debugPrint(const char * fmt, ...)
{
	va_list args;
	va_start(args,fmt);
	char buf[2048];
	vsprintf(buf, fmt, args);
	va_end(args);
    
#ifdef WIN32
	OutputDebugString(buf);
#else
	fprintf(stderr, "%s", buf);
#endif		
}

void g3::g3_debugPrint(std::wstring fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	g3_debugPrint(fmt.c_str(), args);
	va_end(args);
}
void g3::g3_debugPrint(const wchar_t * fmt, ...)
{
	va_list args;
	va_start(args,fmt);
	wchar_t buf[2048];
	vswprintf(buf, 2048, fmt, args);
	va_end(args);
    
#ifdef WIN32
	OutputDebugStringW(buf);
#else
	fwprintf(stderr, L"%s", buf);
#endif		
}