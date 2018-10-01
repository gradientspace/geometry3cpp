

file ( GLOB g3_geomroot "*.cpp" "*.h" )
file ( GLOB g3_geom "geometry/*.*" )
file ( GLOB g3_util "util/*.*" )
file ( GLOB g3_mesh "mesh/*.*" )
file ( GLOB g3_cinterface "cinterface/*.*" )


add_library ( geometry3 SHARED ${g3_geomroot} ${g3_geom} ${g3_util} ${g3_mesh} ${g3_cinterface} )
SET_TARGET_PRECOMPILED_HEADER ( geometry3 geometry3PCH.h geometry3PCH.cpp )



target_compile_definitions(geometry3 PRIVATE GEOMETRY3_DLL_EXPORT)

target_include_directories ( geometry3 PRIVATE . )
target_include_directories ( geometry3 PUBLIC "geometry" )
target_include_directories ( geometry3 PUBLIC "util" )
target_include_directories ( geometry3 PUBLIC "mesh" )
target_include_directories ( geometry3 PUBLIC "cinterface" )

target_link_libraries ( geometry3 GTEngine )
target_link_libraries ( geometry3 WmlCore )
target_link_libraries ( geometry3 WmlMathematics )

target_link_libraries ( geometry3 Eigen3::Eigen)


source_group ( "util" FILES ${g3_util} )
source_group ( "mesh" FILES ${g3_mesh} )
source_group ( "geometry" FILES ${g3_geom} )
source_group ( "core" FILES ${g3_geomroot} )
source_group ( "cinterface" FILES ${g3_cinterface} )



# [RMS] disabled for now
# include boost
#target_include_directories ( geometry3 PUBLIC ${F3Boost_Include})

# [RMS] disabled for now
# include TBB
#target_include_directories ( geometry3 PRIVATE ${F3TBB_Include} )
