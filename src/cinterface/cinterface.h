#pragma once


#ifdef WIN32
#define g3ExternalC __declspec(dllexport)
#else
#define g3ExternalC
#endif


// types
typedef void * CMeshHandle;
