#pragma once

/*
 * [RMS] this file defines and/or pre-declares many core types for the geometry3 library
 *   Note that order-of-include matters, because of nested includes, so some of the .h files
 *   are only included further down...
 */

#define EIGEN_MATRIXBASE_PLUGIN "eigen_MatrixBase_AddOns.h"
#include <Eigen/Core>
#include <Eigen/Geometry>   // required for MatrixBase.cross()  !!

#include <Wm5Core.h>
#include <Wm5Mathematics.h>

#include <memory>
#include <limits>

#ifdef WIN32
    #ifdef GEOMETRY3_DLL_EXPORT
    #define g3External   __declspec( dllexport )
    #else
    #define g3External   __declspec( dllimport )
    #endif
#else
    #define g3External
#endif



namespace Wml = Wm5;

namespace g3 
{
	static constexpr int InvalidID = -1;
	static constexpr int NonManifoldID = -2;
	static constexpr int MarkerID1 = -10;
	static constexpr int MarkerID2 = -11;
	static constexpr int MarkerID3 = -12;
	typedef int VertexID;
	typedef int TriangleID;
	typedef int EdgeID;

	typedef unsigned int GroupID;
	static constexpr GroupID InvalidGroupID = (1<<30);		// currently GMesh group IDs are only 24 bits

	template<typename T> using Math = Wml::Math<T>;
	typedef Wml::Math<float> Mathf;
	typedef Wml::Math<double> Mathd;

	typedef Wml::Interval1d Interval1d;

	template <typename T, int D> using Vector = Eigen::Matrix<T, D, 1>;
	template <typename T, int A> using SquareMatrix = Eigen::Matrix<T, A, A>;

	template<typename T> using Vector2 = Vector<T, 2>;
	typedef Vector2<float> Vector2f;
	typedef Vector2<double> Vector2d;

	template<typename T> using Vector3 = Vector<T, 3>;
	typedef Vector3<float> Vector3f;
	typedef Vector3<double> Vector3d;

	template<typename T> using Vector4 = Vector<T, 4>;
	typedef Vector4<float> Vector4f;
	typedef Vector4<double> Vector4d;

	template<typename T> using Matrix2 = SquareMatrix<T, 2>;
	typedef Matrix2<float> Matrix2f;
	typedef Matrix2<double> Matrix2d;

	template<typename T> using Matrix3 = SquareMatrix<T, 3>;
	typedef Matrix3<float> Matrix3f;
	typedef Matrix3<double> Matrix3d;

	template<typename T> using Matrix4 = SquareMatrix<T, 4>;
	typedef Matrix4<float> Matrix4f;
	typedef Matrix4<double> Matrix4d;

	typedef Vector3<float> Color3f;
	typedef Vector3<double> Color3d;
	typedef Vector4<float> Color4f;
	typedef Vector4<double> Color4d;

	typedef Vector2<int> Vector2i;
	typedef Vector3<int> Vector3i;
	typedef Vector4<int> Vector4i;

	typedef Vector2<int> Index2i;
	typedef Vector3<int> Index3i;
	typedef Vector4<int> Index4i;


	typedef Wml::Color4b Color4b;



	template<typename T> using AxisAlignedBox2 = Wml::AxisAlignedBox2<T>;
	typedef Wml::AxisAlignedBox2<float> AxisAlignedBox2f;
	typedef Wml::AxisAlignedBox2<double> AxisAlignedBox2d;

	template<typename T> using AxisAlignedBox3 = Wml::AxisAlignedBox3<T>;
	typedef Wml::AxisAlignedBox3<float> AxisAlignedBox3f;
	typedef Wml::AxisAlignedBox3<double> AxisAlignedBox3d;

	template<typename T> using Box2 = Wml::Box2<T>;
	typedef Wml::Box2<float> Box2f;
	typedef Wml::Box2<double> Box2d;

	template<typename T> using Box3 = Wml::Box3<T>;
	typedef Wml::Box3<float> Box3f;
	typedef Wml::Box3<double> Box3d;


	class GeometryAssembly;
	typedef std::shared_ptr<GeometryAssembly> GeometryAssemblyPtr;

	class DMesh3;
	typedef std::shared_ptr<DMesh3> DMesh3Ptr;

}

// interface files
#include <SpatialQueryTypes.h>


// [RMS] why is this here? doesn't seem like a good idea...
#include <Frame3.h>


std::ostream& operator<<(std::ostream& os, const g3::Vector2f & v);
std::ostream& operator<<(std::ostream& os, const g3::Vector2d & v);
std::ostream& operator<<(std::ostream& os, const g3::Vector3f & v);
std::ostream& operator<<(std::ostream& os, const g3::Vector3d & v);
std::ostream& operator<<(std::ostream& os, const g3::Vector4f & v);
std::ostream& operator<<(std::ostream& os, const g3::Vector4d & v);
std::ostream& operator<<(std::ostream& os, const g3::Vector2i & v);
std::ostream& operator<<(std::ostream& os, const g3::Vector3i & v);
std::ostream& operator<<(std::ostream& os, const g3::Vector4i & v);

std::ostream& operator<<(std::ostream& os, const g3::AxisAlignedBox2f & b);
std::ostream& operator<<(std::ostream& os, const g3::AxisAlignedBox2d & b);
std::ostream& operator<<(std::ostream& os, const g3::AxisAlignedBox3f & b);
std::ostream& operator<<(std::ostream& os, const g3::AxisAlignedBox3d & b);

