#
# Top-level CMake file for geometry3 library
# Ryan Schmidt - rms@rms80.com
#

cmake_minimum_required (VERSION 3.4.0)
project (geometry3)


set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})


# PCH macros:
#  - SET_TARGET_PRECOMPILED_HEADER ( target pch.h pch.cpp )
include ("cmake/PCHSupportV3.cmake")

# Frame3 custom cmake macros
#include ("cmake/F3UseThirdPartyShared.cmake")
#include ("cmake/F3CopyThirdPartyFiles.cmake")
#include ("cmake/F3CopyThirdPartyFolder.cmake")
#include ("cmake/F3AddForcedInclude.cmake")


if (APPLE)
	# tell xcode to use C++14 language & std library
	set (CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++14")
	set (CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")	
endif ()


# set up 3P stuff - root, platform, architecture
#set ( F3ThirdPartyRoot ${CMAKE_CURRENT_SOURCE_DIR}/../3P )

set ( G3Arch "64" )
if(MSVC)
	# this is the path that we find 3P/s in
	set ( G3Platform "windows" )
else()
	# OSX
	set ( G3Platform "osx" )
endif()
message ( "geometry3 Platform is ${G3Platform} Arch is ${G3Arch}" )


# boost third-party (cmake find_package works for this)
#set ( F3Boost_Root ${F3ThirdPartyRoot}/boost/${F3Platform} )
#if(MSVC)
#	set ( F3Boost_Lib ${F3Boost_Root}/lib64 )
#	set ( F3Boost_Bin ${F3Boost_Root}/lib64 )
#else ()
#	set ( F3Boost_Lib ${F3Boost_Root}/lib )
#endif ()
#set ( F3Boost_Include ${F3Boost_Root}/include )	



# tinylibs third-party
#set ( tinylibs_Root ${F3ThirdPartyRoot}/tinylibs/${F3Platform} )
#set ( tinylibs_Include ${F3ThirdPartyRoot}/tinylibs/include )


# Intel TBB third-party
set ( F3TBB_Root ${F3ThirdPartyRoot}/tbb/${F3Platform} )
if (MSVC) 
   set ( F3TBB_Lib ${F3TBB_Root}/lib64 )
   set ( F3TBB_Bin ${F3TBB_Root}/bin64 )
else ()
   set ( F3TBB_Lib ${F3TBB_Root}/lib )
endif ()
set ( F3TBB_Include ${F3TBB_Root}/include )	




#add_subdirectory (GTEngine)
add_subdirectory (src)
add_subdirectory (external)
add_subdirectory (test)
