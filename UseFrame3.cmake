#
# [RMS] set up frame3 library
#   Requires the following variables be set:
#     FRAME3_SOURCE_DIR 
#     THREEP_DIR
#


# PCH macros:
#  - SET_TARGET_PRECOMPILED_HEADER ( target pch.h pch.cpp )
include ("${FRAME3_SOURCE_DIR}/cmake/PCHSupportV3.cmake")

# Frame3 custom cmake macros
include ("${FRAME3_SOURCE_DIR}/cmake/F3UseThirdPartyShared.cmake")
include ("${FRAME3_SOURCE_DIR}/cmake/F3UseThirdPartyStatic.cmake")
include ("${FRAME3_SOURCE_DIR}/cmake/F3CopyThirdPartyFiles.cmake")
include ("${FRAME3_SOURCE_DIR}/cmake/F3CopyThirdPartyFolder.cmake")
include ("${FRAME3_SOURCE_DIR}/cmake/F3AddForcedInclude.cmake")


if (APPLE)
	# tell xcode to use C++11 language & std library
	set (CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++11")
	set (CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")	
endif ()


# set up 3P stuff - root, platform, architecture
set ( F3ThirdPartyRoot ${THREEP_DIR} )

set ( F3Arch "64" )
if(MSVC)
	# this is the path that we find 3P/s in
	set ( F3Platform "windows" )
else()
	# OSX
	set ( F3Platform "osx" )
endif()
message ( "[Frame3] Platform is ${F3Platform} Arch is ${F3Arch}" )


# boost third-party (cmake find_package works for this)
set ( F3Boost_Root ${F3ThirdPartyRoot}/boost/${F3Platform} )
if(MSVC)
	set ( F3Boost_Lib ${F3Boost_Root}/lib64 )
	set ( F3Boost_Bin ${F3Boost_Root}/lib64 )
else ()
	set ( F3Boost_Lib ${F3Boost_Root}/lib )
endif ()
set ( F3Boost_Include ${F3Boost_Root}/include )	
message ( "[Frame3] F3Boost_Include is ${F3Boost_Include}" )


# glm third-party
set ( glm_Root ${F3ThirdPartyRoot}/glm/${F3Platform} )
set ( glm_Include ${glm_Root}/include )


# tinylibs third-party
set ( tinylibs_Root ${F3ThirdPartyRoot}/tinylibs/${F3Platform} )
set ( tinylibs_Include ${F3ThirdPartyRoot}/tinylibs/include )

	
# [RMS] we have to put win64 and win32 in separate folders because otherwise
#  CMake's Qt support cannot find files (it uses relative paths generated during build)
if (MSVC) 
	if ( ${F3Arch} EQUAL "64" )
		set ( Qt_Root ${F3ThirdPartyRoot}/qt/win64 )
	else ()
		set ( Qt_Root ${F3ThirdPartyRoot}/qt/win32 )
	endif ()
	set ( Qt_Lib ${Qt_Root}/lib )
	set ( Qt_Bin ${Qt_Root}/bin )
	set ( Qt_Plugins ${Qt_Root}/plugins )
	set ( Qt_Include ${Qt_Root}/include )
else()
	set ( Qt_Root ${F3ThirdPartyRoot}/qt/osx )
	set ( Qt_Lib ${Qt_Root}/lib )
	set ( Qt_Bin ${Qt_Root}/bin )
	set ( Qt_Plugins ${Qt_Root}/plugins )
	set ( Qt_Include ${Qt_Root}/include )
endif ()
message ( "[Frame3] Qt_Lib is ${Qt_Lib}" )
message ( "[Frame3] Qt_Bin is ${Qt_Bin}" )


# Intel TBB third-party
set ( F3TBB_Root ${F3ThirdPartyRoot}/tbb/${F3Platform} )
if (MSVC) 
   set ( F3TBB_Lib ${F3TBB_Root}/lib64 )
   set ( F3TBB_Bin ${F3TBB_Root}/bin64 )
else ()
   set ( F3TBB_Lib ${F3TBB_Root}/lib )
endif ()
set ( F3TBB_Include ${F3TBB_Root}/include )	


# python third-party
set ( F3Python_Root ${F3ThirdPartyRoot}/python/${F3Platform} )
if (MSVC) 
   set ( F3Python_Lib ${F3Python_Root}/lib64 )
   set ( F3Python_Bin ${F3Python_Root}/lib64 )
else ()
   set ( F3Python_Lib ${F3Python_Root}/lib )
endif ()
set ( F3Python_Include ${F3Python_Root}/include )	