

file ( GLOB code "*.cpp" "*.h" )

add_executable ( geometry3_tests ${code}  )

# geometry3 project now "exports" these so they are automatically picked up here
#target_compile_definitions(geometry3_tests PRIVATE GEOMETRY3_DLL_IMPORT)
#target_compile_definitions(geometry3_tests PRIVATE G3_STATIC_LIB)

target_link_libraries ( geometry3_tests geometry3 )

