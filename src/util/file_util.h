#pragma once


#ifdef _WIN32
#include <io.h> 
#define access    _access_s			// FileExists
#define FILE_SEPARATOR "\\"
#else
#include <unistd.h>					// FileExists
#define FILE_SEPARATOR "/"
#endif


namespace g3
{
class FileUtil
{
public:
	FileUtil() = delete;


	static bool FileExists(const std::string &Filename)
	{
		return access(Filename.c_str(), 0) == 0;
	}


	static std::string PathCombine(const std::string & s1, const std::string & s2)
	{
		// AAAHHH NOOOO
		return s1 + FILE_SEPARATOR + s2;
	}


};



}