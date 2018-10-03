
set ( LIBMATHNAME "WmlMathematics" )

file ( GLOB math_root "*.cpp" "*.h" )
file ( GLOB algebra "Algebra/*.*" )
file ( GLOB approx "Approximation/*.*" )
file ( GLOB base  "Base/*.*" )
file ( GLOB compgeom "ComputationalGeometry/*.*" )
file ( GLOB contain "Containment/*.*" )
file ( GLOB curvsurf "CurvesSurfacesVolumes/*.*" )
file ( GLOB dist "Distance/*.*" )
file ( GLOB interp "Interpolation/*.*" )
file ( GLOB inters "Intersection/*.*" )
file ( GLOB misc "Miscellaneous/*.*" )
file ( GLOB numerics "NumericalAnalysis/*.*" )
file ( GLOB objects2D "Objects2D/*.*" )
file ( GLOB objects3D "Objects3D/*.*" )
file ( GLOB query "Query/*.*" )
file ( GLOB rational "Rational/*.*" )


add_library ( ${LIBMATHNAME} ${math_root} ${algebra} ${approx} ${base} ${compgeom} ${contain} 
			  ${misc} ${curvsurf} ${dist} ${interp} ${inters} ${numerics} ${objects2D} ${objects3D} 
			  ${query} ${rational} )

add_definitions (-DWM5_USE_PRECOMPILED_HEADERS)
SET_TARGET_PRECOMPILED_HEADER ( ${LIBMATHNAME} Wm5MathematicsPCH.h Wm5MathematicsPCH.cpp )

target_link_libraries ( ${LIBMATHNAME} WmlCore )


target_include_directories ( ${LIBMATHNAME} PUBLIC "." )
target_include_directories ( ${LIBMATHNAME} PUBLIC "Algebra" )
target_include_directories ( ${LIBMATHNAME} PUBLIC "Approximation" )
target_include_directories ( ${LIBMATHNAME} PUBLIC "Base" )
target_include_directories ( ${LIBMATHNAME} PUBLIC "ComputationalGeometry" )
target_include_directories ( ${LIBMATHNAME} PUBLIC "Containment" )
target_include_directories ( ${LIBMATHNAME} PUBLIC "CurvesSurfacesVolumes" )
target_include_directories ( ${LIBMATHNAME} PUBLIC "Distance" )
target_include_directories ( ${LIBMATHNAME} PUBLIC "Interpolation" )
target_include_directories ( ${LIBMATHNAME} PUBLIC "Intersection" )
target_include_directories ( ${LIBMATHNAME} PUBLIC "Miscellaneous" )
target_include_directories ( ${LIBMATHNAME} PUBLIC "NumericalAnalysis" )
target_include_directories ( ${LIBMATHNAME} PUBLIC "Objects2D" )
target_include_directories ( ${LIBMATHNAME} PUBLIC "Objects3D" )
target_include_directories ( ${LIBMATHNAME} PUBLIC "Query" )
target_include_directories ( ${LIBMATHNAME} PUBLIC "Rational" )

# [RMS] configure source groups for visual studio project

source_group ( "Algebra" FILES ${algebra} )
source_group ( "Approximation" FILES ${approx} )
source_group ( "Base" FILES ${base} )
source_group ( "ComputationalGeometry" FILES ${compgeom} )
source_group ( "Containment" FILES ${contain} )
source_group ( "CurvesSurfacesVolumes" FILES ${curvsurf} )
source_group ( "Distance" FILES ${dist} )
source_group ( "Interpolation" FILES ${interp} )
source_group ( "Intersection" FILES ${inters} )
source_group ( "Miscellaneous" FILES ${core} )
source_group ( "NumericalAnalysis" FILES ${numerics} )
source_group ( "Objects2D" FILES ${objects2D} )
source_group ( "Objects3D" FILES ${objects3D} )
source_group ( "Query" FILES ${query} )
source_group ( "Rational" FILES ${rational} )



