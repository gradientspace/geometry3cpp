#pragma once

/*
 * [RMS] this file defines and/or pre-declares many core types for the geometry3 library
 *   Note that order-of-include matters, because of nested includes, so some of the .h files
 *   are only included further down...
 */

#include <Wm5Core.h>
#include <Wm5Mathematics.h>
#include <memory>
#include <limits>


#ifdef GEOMETRY3_DLL_EXPORT
#define g3External   __declspec( dllexport )  
#else
#define g3External   __declspec( dllimport )  
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


	template<typename T> using Vector2 = Wml::Vector2<T>;
	typedef Wml::Vector2<float> Vector2f;
	typedef Wml::Vector2<double> Vector2d;

	template<typename T> using Vector3 = Wml::Vector3<T>;
	typedef Wml::Vector3<float> Vector3f;
	typedef Wml::Vector3<double> Vector3d;
	typedef Wml::Vector3<float> Color3f;
	typedef Wml::Vector3<double> Color3d;

	template<typename T> using Vector4 = Wml::Vector4<T>;
	typedef Wml::Vector4<float> Vector4f;
	typedef Wml::Vector4<double> Vector4d;
	typedef Wml::Vector4<float> Color4f;
	typedef Wml::Vector4<double> Color4d;

	typedef Wml::IVector2 Vector2i;
	typedef Wml::IVector3 Vector3i;
	typedef Wml::IVector4 Vector4i;
	typedef Wml::Color4b Color4b;

	template<typename T> using Matrix2 = Wml::Matrix2<T>;
	typedef Wml::Matrix2<float> Matrix2f;
	typedef Wml::Matrix2<double> Matrix2d;

	template<typename T> using Matrix3 = Wml::Matrix3<T>;
	typedef Wml::Matrix3<float> Matrix3f;
	typedef Wml::Matrix3<double> Matrix3d;

	template<typename T> using Matrix4 = Wml::Matrix4<T>;
	typedef Wml::Matrix4<float> Matrix4f;
	typedef Wml::Matrix4<double> Matrix4d;

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
