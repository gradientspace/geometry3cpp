
set ( LIBMATHNAME "WmlMathematics" )

file ( GLOB math_root "*.cpp" "*.h" )
file ( GLOB approx "Approximation/*.*" )
file ( GLOB compgeom "ComputationalGeometry/*.*" )
file ( GLOB contain "Containment/*.*" )
file ( GLOB core "Core/*.*" )
file ( GLOB curvsurf "CurvesSurfacesVolumes/*.*" )
file ( GLOB dist "Distance/*.*" )
file ( GLOB interp "Interpolation/*.*" )
file ( GLOB inters "Intersection/*.*" )
file ( GLOB numerics "NumericalAnalysis/*.*" )
file ( GLOB objects "Objects/*.*" )


add_library ( ${LIBMATHNAME} ${math_root} ${approx} ${compgeom} ${contain} 
			  ${core} ${curvsurf} ${dist} ${interp} ${inters} ${numerics} ${objects} )

add_definitions (-DWM5_USE_PRECOMPILED_HEADERS)
SET_TARGET_PRECOMPILED_HEADER ( ${LIBMATHNAME} Wm5MathematicsPCH.h Wm5MathematicsPCH.cpp )

target_link_libraries ( ${LIBMATHNAME} WmlCore )


target_include_directories ( ${LIBMATHNAME} PUBLIC "." )
target_include_directories ( ${LIBMATHNAME} PUBLIC "Approximation" )
target_include_directories ( ${LIBMATHNAME} PUBLIC "ComputationalGeometry" )
target_include_directories ( ${LIBMATHNAME} PUBLIC "Containment" )
target_include_directories ( ${LIBMATHNAME} PUBLIC "Core" )
target_include_directories ( ${LIBMATHNAME} PUBLIC "CurvesSurfacesVolumes" )
target_include_directories ( ${LIBMATHNAME} PUBLIC "Distance" )
target_include_directories ( ${LIBMATHNAME} PUBLIC "Interpolation" )
target_include_directories ( ${LIBMATHNAME} PUBLIC "Intersection" )
target_include_directories ( ${LIBMATHNAME} PUBLIC "NumericalAnalysis" )
target_include_directories ( ${LIBMATHNAME} PUBLIC "Objects" )

# [RMS] configure source groups for visual studio project

source_group ( "Core" FILES ${core} )
source_group ( "Approximation" FILES ${approx} )
source_group ( "ComputationalGeometry" FILES ${compgeom} )
source_group ( "Containment" FILES ${contain} )
source_group ( "CurvesSurfacesVolumes" FILES ${curvsurf} )
source_group ( "Distance" FILES ${dist} )
source_group ( "Interpolation" FILES ${interp} )
source_group ( "Intersection" FILES ${inters} )
source_group ( "NumericalAnalysis" FILES ${numerics} )
source_group ( "Objects" FILES ${objects} )



