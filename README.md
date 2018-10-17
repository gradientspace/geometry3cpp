# geometry3cpp

Open-Source (Boost-license) C++ library for geometric computing. 

geometry3cpp is an in-progress port of [geometry3Sharp](https://github.com/gradientspace/geometry3Sharp), the gradientspace C# library for geometric computing. Except the parts of that library that were ported from the C++ WildMagic/GTEngine libraries, which are included in this one. It's a bit confusing.

**WORK IN PROGRESS - HARDLY TESTED - BUYER BEWARE**


# Goals

g3cpp is intended to be a general-purpose high-level geometric computing package, with a focus on triangle mesh processing. It includes *other* math libraries which provide most of the low-level vector math stuff, solvers, etc. 

I would like the library to be a header-only library. *However* the WildMagic5 dependency is currently structured to produce a compiled DLL, and GTEngine is also not fully header-only, there are few .cpp files. My intention is to eventually refactor these libraries so that they are fully header-only as well. 

*Note: I would welcome comments about whether a header-only approach is suitable/desirable for production use.*

# Current State

The dependencies contain an enormous amount of functionality, much more than geometry3Sharp, at the lower level. At the mesh level the following classes have been ported

* **DMesh3** - dynamic mesh, fully ported
* **DMeshAABBTree3** - AABB bounding volume hierarchy for DMesh3, only nearest-point and generic traversal currently ported
* **Remesher** - majority ported, no parallel smoothing/projection currently (avoiding tbb dependency for now)
* **OBJReader** and **OBJWriter** - functional but support for materials/texture maps is not fully ported

**Old Code** This repository has evolved from an earlier attempt (pre-geometry3Sharp) at a mesh processing library. This previous code has not been fully deleted/updated, but will be. 

**goemetry3_tests** project has some basic sample code in **main.cpp**. The CMake generators add this as a separate project, and it generates a standalone executable linking to g3Cpp. If you want to just try something out, I would recommend starting here.


# Dependencies

All dependencies are included in the repository, for convenience.

1) [**Eigen** ](https://eigen.tuxfamily.org/), the C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms (according to their website). **MPL2 License**. Header-only, source included in */external/Eigen* for ease of compiling, but you could probably get CMake to look somewhere else. The Vector types (Vector2d, Vector3d, etc) in g3cpp are Eigen vectors.

2) **WildMagic5** from [GeometricTools](https://www.geometrictools.com/), written by David Eberly. **Boost License**. Only the LibCore and LibMathematics components. Vector math, Geometric intersection and distance tests in 2D and 3D, containment fitters, geometric approximations fitters, Computational Geometry algorithms, Numerical methods, rational number types, 1/2/3D interpolation methods. It's amazing, I've been using WildMagic for 10+ years, since version 2. The source is included in the */external/* subdirectory. This library is no longer maintained and so I have made various local changes to ease porting from the C# version and make it easier to pass vector types between Eigen and Wm5. 

3) **GTEngine** also from [GeometricTools](https://www.geometrictools.com/) and written by David Eberly. **Boost License**. Only the Math parts, again. This is the "next version" of WildMagic5, however quite a few of the classes from Wm5 were not ported. However it has some new things like Constrained Delaunay. So both are included. The source is included in */external/GTEngine*. Namespace for this one is "gte". *(Currently not actually called by the g3cpp code, but it will be)



# Building

**Currently only definitely working on Windows, but the code should be portable..."should"...***

1) Run **bin\get_cmake.bat** to download and unzip the required version of CMake. This will be unzipped in bin\cmake_win, it will not affect your system installation of CMake (if installed). 

2) Run **external\cmake_install_eigen.bat**. This does some cmake stuff so that the Eigen CMake module can be found. Again, it doesn't affect any other Eigen installations/etc. *(Not entirely clear why this is necessary, but CMake cannot find Eigen without it)*

3) Run top-level **StartVS2017_Debug.bat** or **StartVS2017_Release.bat** to run CMake and automatically open the resulting VS2017 solution. 

If you are actively working on the g3cpp code, UpdateVS_Debug_2017.bat will regenerate the solution file, which will cause VS2017 to prompt you to reload it. 


# Packaging / Build Output

Currently a dll, output to visual studio build folders, ie *build\Win64_Debug* or *_Release*. 

You can switch to a static .lib by editing *src\CMakeLists.txt* and commenting/uncommenting the relevant block (see comments in file)

No packaging as of yet. You will have to add the various header source folders to other projects. Sorry!! (happy to take PRs here!)


# Testing

There isn't any, yet.
