

file ( GLOB code "*.cpp" "*.h" )

add_executable ( geometry3_tests ${code}  )

target_compile_definitions(geometry3_tests PRIVATE GEOMETRY3_DLL_IMPORT)

target_link_libraries ( geometry3_tests geometry3 )

