

file ( GLOB HDR_BASE "Include/*.h" )
file ( GLOB HDR_MATH "Include/Mathematics/*.h" )
file ( GLOB HDR_IMG "Include/Imagics/*.h" )
file ( GLOB HDR_LOW "Include/LowLevel/*.h" )

file ( GLOB_RECURSE SRC_MATH "Source/Mathematics/*.cpp" )
file ( GLOB_RECURSE SRC_IMG "Source/Imagics/*.cpp" )

add_library ( GTEngine ${HDR_BASE} ${HDR_MATH} ${HDR_IMG} ${HDR_LOW} ${SRC_MATH} ${SRC_IMG} )

target_include_directories ( GTEngine PUBLIC "Include" )
target_include_directories ( GTEngine PUBLIC "Include/Mathematics" )
target_include_directories ( GTEngine PUBLIC "Include/Imagics" )
target_include_directories ( GTEngine PUBLIC "Include/LowLevel" )


# [RMS] configure source groups for visual studio project

source_group ( "Imagics" FILES ${SRC_IMG} ${HDR_IMG} )

set(GTEMATH "*/Mathematics")


file ( GLOB GTE_APPR "Include/Mathematics/GteAppr*.h"  )
file ( GLOB GTE_DIST "Include/Mathematics/GteDist*.h"  )
file ( GLOB GTE_INTP "Include/Mathematics/GteIntp*.h"  )
file ( GLOB GTE_INTR "Include/Mathematics/GteIntr*.h"  )
file ( GLOB GTE_CONT "Include/Mathematics/GteCont*.h"  )
file ( GLOB GTE_CONT "Include/Mathematics/GteCont*.h"  )
file ( GLOB GTE_SOLV "${GTEMATH}/GteMin*.h" 
					 "${GTEMATH}/GteOde*.h" 
					 "${GTEMATH}/GteRoots*.h"
					 "${GTEMATH}/Gte*solver*.h"
					 "${GTEMATH}/GteSingular*.h"
					 "${GTEMATH}/GteSingular*.h"
					 "${GTEMATH}/Gte*Elim*.h"
					 "${GTEMATH}/Gte*Linear*.h"
					 )
file ( GLOB GTE_COMP "${GTEMATH}/Gte*Delaunay*.h" 
					 "${GTEMATH}/Gte*Hull*.h" 
					 "${GTEMATH}/Gte*Query*.h" 
					 "${GTEMATH}/Gte*Key*.*" 
					 "${GTEMATH}/Gte*Triangulate*.*" 
					 )					 
file ( GLOB GTE_CURV "${GTEMATH}/Gte*Mesh*.*" 
					 "${GTEMATH}/Gte*Spline*.h" 
					 "${GTEMATH}/Gte*NURBS*.h" 
					 "${GTEMATH}/Gte*Bezier*.h" 
					 "${GTEMATH}/Gte*Frame*.h" 
					 "${GTEMATH}/Gte*Curve.h" 
					 "${GTEMATH}/Gte*Surface.h" 
					 )					 

source_group ( "Mathematics" FILES ${SRC_MATH} ${HDR_MATH} )
source_group ( "Approximation" FILES ${GTE_APPR} )
source_group ( "Distance" FILES ${GTE_DIST} )
source_group ( "Interpolation" FILES ${GTE_INTP} )
source_group ( "Intersection" FILES ${GTE_INTR} )
source_group ( "Containment" FILES ${GTE_CONT} )
source_group ( "Solving" FILES ${GTE_SOLV} )
source_group ( "CompGeometry" FILES ${GTE_COMP} )
source_group ( "CurvesSurfaces" FILES ${GTE_CURV} )


file ( GLOB GTE_JUNK "${GTEMATH}/GteIEEE*.*" 
					 "${GTEMATH}/GteUInteger*.*"
					 "${GTEMATH}/GteIntel*.*"
					 "${GTEMATH}/GteBitHacks*.*"
					 "${GTEMATH}/Gte*Estimate*.h"
					 )
source_group ( "LowLevel" FILES ${HDR_LOW} ${GTE_JUNK} )
