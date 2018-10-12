# geometry3 / g3cpp

Open-Source (Boost-license) C++ library for geometric computing. 

geometry3 is an in-progress port of [geometry3Sharp](https://github.com/gradientspace/geometry3Sharp), the gradientspace C# library for geometric computing. Except the parts of that library that were ported from the C++ WildMagic/GTEngine libraries, which are included in this one. It's a bit confusing.

**WORK IN PROGRESS - HARDLY TESTED - BUYER BEWARE**


# Goals

g3cpp is intended to be a general-purpose high-level geometric computing package, with a focus on triangle mesh processing. It includes *other* math libraries which provide most of the low-level vector math stuff, solvers, etc. 

I would like the library to be a header-only library. *However* the WildMagic5 dependency is currently structured to produce a compiled DLL, and GTEngine is also not fully header-only, there are few .cpp files. My intention is to eventually refactor these libraries so that they are fully header-only as well. 



# Current State

The dependencies contain an enormous amount of functionality, much more than geometry3Sharp, at the lower level. At the mesh level the following classes have been ported

* **DMesh3** - dynamic mesh, fully ported
* **DMeshAABBTree3** - AABB bounding volume hierarchy for DMesh3, only nearest-point and generic traversal currently ported
* **Remesher** - majority ported, no parallel smoothing/projection currently (avoiding tbb dependency for now)
* **OBJReader** and **OBJWriter** - functional but support for materials/texture maps is not fully ported

**Old Code** This repository has evolved from an earlier attempt (pre-geometry3Sharp) at a mesh processing library. This previous code has not been fully deleted/updated, but will be. 

**goemetry3_tests** project has some basic sample code. It is included in the CMake generators but is a separate executable


# Dependencies

All dependencies are included in the repository, for convenience.

1) [**Eigen** ](https://eigen.tuxfamily.org/), the C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms (according to their website). **MPL2 License**. Header-only, source included in */external/Eigen* for ease of compiling, but you could probably get CMake to look somewhere else. The Vector types (Vector2d, Vector3d, etc) in g3cpp are Eigen vectors.

2) **WildMagic5** from [GeometricTools](https://www.geometrictools.com/), written by David Eberly. **Boost License**. Only the LibCore and LibMathematics components. Vector math, Geometric intersection and distance tests in 2D and 3D, containment fitters, geometric approximations fitters, Computational Geometry algorithms, Numerical methods, rational number types, 1/2/3D interpolation methods. It's amazing, I've been using WildMagic for 10+ years, since version 2. The source is included in the */external/* subdirectory. This library is no longer maintained and so I have made various local changes to ease porting from the C# version and make it easier to pass vector types between Eigen and Wm5. 

3) **GTEngine** also from [GeometricTools](https://www.geometrictools.com/) and written by David Eberly. **Boost License**. Only the Math parts, again. This is the "next version" of WildMagic5, however quite a few of the classes from Wm5 were not ported. However it has some new things like Constrained Delaunay. So both are included. The source is included in */external/GTEngine*. Namespace for this one is "gte". *(Currently not actually called by the g3cpp code, but it will be)



# Building

**Currently only definitely working on Windows, but the code should be portable..."should"...***

a CMake build system is provided. Run top-level StartVS2017_Debug.bat or StartVS2017_Release.bat to automatically run CMake and open the resulting VS2017 solution. 

If you are actively working on the g3cpp code, UpdateVS_Debug_2017.bat will regenerate the solution file, which will cause VS2017 to prompt you to reload it. 


# Build Output

Currently a dll.


# Testing

There isn't any, yet.
