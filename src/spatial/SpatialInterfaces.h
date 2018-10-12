#pragma once

#include <g3types.h>

namespace g3
{


class ISpatial
{
	virtual bool SupportsPointContainment() = 0;

	/// <summary> return true if query point is inside object </summary>
	virtual bool IsInside(const Vector3d & p) = 0;
};



class IMeshSpatial : public ISpatial
{
public:
	virtual bool SupportsNearestTriangle() = 0;

	/// <summary>
	/// Find id of triangle nearest to p within distance fMaxDist, or return DMesh3.InvalidID if not found
	/// </summary>
	virtual int FindNearestTriangle(const Vector3d & p, double & fNearestDistSqr, double fMaxDist = std::numeric_limits<double>::max()) = 0;


	virtual bool SupportsTriangleRayIntersection() = 0;

	/// <summary>
	/// Find id of triangle intersected by ray, where intersection point is within distance fMaxDist, or return DMesh3.InvalidID if not found
	/// </summary>
	virtual int FindNearestHitTriangle(const Ray3d & ray, double fMaxDist = std::numeric_limits<double>::max()) = 0;
};


class IProjectionTarget
{
public:
	virtual Vector3d Project(const Vector3d & vPoint, int identifier = -1) = 0;
};

class IOrientedProjectionTarget : public IProjectionTarget
{
public:
	virtual Vector3d Project(const Vector3d & vPoint, Vector3d & vProjectNormal, int identifier = -1) = 0;
};

class IIntersectionTarget
{
public:
	virtual bool HasNormal() = 0;
	virtual bool RayIntersect(const Ray3d & ray, Vector3d & vHit, Vector3d & vHitNormal) = 0;
};


}