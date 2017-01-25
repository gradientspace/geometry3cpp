
set ( LIBCORENAME "WmlCore" )

file ( GLOB wml_root "*.cpp" "*.h" )
file ( GLOB wml_io "IO/*.*" )
file ( GLOB wml_utility "Utility/*.*" )

add_library ( ${LIBCORENAME} ${wml_root} ${wml_io} ${wml_utility} )

add_definitions (-DWM5_USE_PRECOMPILED_HEADERS)
SET_TARGET_PRECOMPILED_HEADER ( ${LIBCORENAME} Wm5CorePCH.h Wm5CorePCH.cpp )

target_include_directories ( ${LIBCORENAME} PUBLIC "." )
target_include_directories ( ${LIBCORENAME} PUBLIC "IO" )
target_include_directories ( ${LIBCORENAME} PUBLIC "Utility" )


# [RMS] configure source groups for visual studio project

source_group ( "IO" FILES ${wml_io} )
source_group ( "Utility" FILES ${wml_utility} )


