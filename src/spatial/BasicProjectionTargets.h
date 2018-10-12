#pragma once

#include <g3types.h>
#include <SpatialInterfaces.h>
#include <DMeshAABBTree3.h>

namespace g3
{


/// <summary>
/// MeshProjectionTarget provides an IProjectionTarget interface to a mesh + spatial data structure.
/// Use to project points to mesh surface.
/// </summary>
class MeshProjectionTarget : public IOrientedProjectionTarget
{
public:
	DMesh3Ptr Mesh;
	IMeshSpatialPtr Spatial;

	MeshProjectionTarget() { }

	MeshProjectionTarget(DMesh3Ptr mesh, IMeshSpatialPtr spatial)
	{
		Mesh = mesh;
		Spatial = spatial;
		if (Spatial == nullptr)
			Spatial = std::make_shared<DMeshAABBTree3>(mesh, true);
	}

	MeshProjectionTarget(DMesh3Ptr mesh)
	{
		Mesh = mesh;
		Spatial = std::make_shared<DMeshAABBTree3>(mesh, true);
	}

	virtual Vector3d Project(const Vector3d & vPoint, int identifier = -1) override
	{
		double fDistSqr;
		int tNearestID = Spatial->FindNearestTriangle(vPoint, fDistSqr);
		Vector3d v0, v1, v2;
		Mesh->GetTriVertices(tNearestID, v0, v1, v2);

		Wml::DistPoint3Triangle3d dist(vPoint, Triangle3d(v0, v1, v2));
		double distSqr = dist.GetSquared();
		Vector3d vProj = dist.GetClosestPoint1();;
		return vProj;
	}

	virtual Vector3d Project(const Vector3d & vPoint, Vector3d & vProjectNormal, int identifier = -1)
	{
		double fDistSqr;
		int tNearestID = Spatial->FindNearestTriangle(vPoint, fDistSqr);
		Vector3d v0, v1, v2;
		Mesh->GetTriVertices(tNearestID, v0, v1, v2);

		vProjectNormal = Normal(v0, v1, v2);

		Wml::DistPoint3Triangle3d dist(vPoint, Triangle3d(v0, v1, v2));
		double distSqr = dist.Get();
		return dist.GetClosestPoint1();
	}

	/// <summary>
	/// Automatically construct fastest projection target for mesh
	/// </summary>
	static MeshProjectionTarget Auto(DMesh3Ptr mesh, bool bForceCopy = true)
	{
		if (bForceCopy)
			return MeshProjectionTarget(std::make_shared<DMesh3>(*mesh, false, MeshComponents::None));
		else
			return MeshProjectionTarget(mesh);
	}
	static MeshProjectionTargetPtr AutoPtr(DMesh3Ptr mesh, bool bForceCopy = true)
	{
		if (bForceCopy)
			return std::make_shared<MeshProjectionTarget>(std::make_shared<DMesh3>(*mesh, false, MeshComponents::None));
		else
			return std::make_shared<MeshProjectionTarget>(mesh);
	}

	static MeshProjectionTarget Auto(const DMesh3 & mesh)
	{
		return MeshProjectionTarget(std::make_shared<DMesh3>(mesh, false, MeshComponents::None));
	}
	static MeshProjectionTargetPtr AutoPtr(const DMesh3 & mesh)
	{
		return std::make_shared<MeshProjectionTarget>(std::make_shared<DMesh3>(mesh, false, MeshComponents::None));
	}


	/// <summary>
	/// Automatically construct fastest projection target for region of mesh
	/// </summary>
	//static MeshProjectionTarget Auto(DMesh3 mesh, IEnumerable<int> triangles, int nExpandRings = 5)
	//{
	//	MeshFaceSelection targetRegion = new MeshFaceSelection(mesh);
	//	targetRegion.Select(triangles);
	//	targetRegion.ExpandToOneRingNeighbours(nExpandRings);
	//	DSubmesh3 submesh = new DSubmesh3(mesh, targetRegion);
	//	return new MeshProjectionTarget(submesh.SubMesh);
	//}
};



//
///// <summary>
///// Extension of MeshProjectionTarget that allows the target to have a transformation
///// relative to it's internal space. Call SetTransform(), or initialize the transforms yourself
///// </summary>
//class TransformedMeshProjectionTarget : MeshProjectionTarget
//{
//	TransformSequence SourceToTargetXForm;
//	TransformSequence TargetToSourceXForm;
//
//	TransformedMeshProjectionTarget() { }
//	TransformedMeshProjectionTarget(DMesh3 mesh, ISpatial spatial) : base(mesh, spatial)
//	{
//	}
//	TransformedMeshProjectionTarget(DMesh3 mesh) : base(mesh)
//	{
//	}
//
//	void SetTransform(TransformSequence sourceToTargetX)
//	{
//		SourceToTargetXForm = sourceToTargetX;
//		TargetToSourceXForm = SourceToTargetXForm.MakeInverse();
//	}
//
//	override Vector3d Project(Vector3d vPoint, int identifier = -1)
//	{
//		Vector3d vTargetPt = SourceToTargetXForm.TransformP(vPoint);
//		Vector3d vTargetProj = base.Project(vTargetPt, identifier);
//		return TargetToSourceXForm.TransformP(vTargetProj);
//	}
//
//
//	override Vector3d Project(Vector3d vPoint, out Vector3d vProjectNormal, int identifier = -1)
//	{
//		Vector3d vTargetPt = SourceToTargetXForm.TransformP(vPoint);
//		Vector3d vTargetProjNormal;
//		Vector3d vTargetProj = base.Project(vTargetPt, out vTargetProjNormal, identifier);
//		vProjectNormal = TargetToSourceXForm.TransformV(vTargetProjNormal).Normalized;
//		return TargetToSourceXForm.TransformP(vTargetProj);
//	}
//};
//
//
//
//
//
//
//
//
//
//class PlaneProjectionTarget : IProjectionTarget
//{
//	Vector3d Origin;
//	Vector3d Normal;
//
//	Vector3d Project(Vector3d vPoint, int identifier = -1)
//	{
//		Vector3d d = vPoint - Origin;
//		return Origin + (d - d.Dot(Normal) * Normal);
//	}
//};
//
//
//
//
//class CircleProjectionTarget : IProjectionTarget
//{
//	Circle3d Circle;
//
//	Vector3d Project(Vector3d vPoint, int identifier = -1)
//	{
//		DistPoint3Circle3 d = new DistPoint3Circle3(vPoint, Circle);
//		d.GetSquared();
//		return d.CircleClosest;
//	}
//};
//
//
//
//class CylinderProjectionTarget : IProjectionTarget
//{
//	Cylinder3d Cylinder;
//
//	Vector3d Project(Vector3d vPoint, int identifer = -1)
//	{
//		DistPoint3Cylinder3 d = new DistPoint3Cylinder3(vPoint, Cylinder);
//		d.GetSquared();
//		return d.CylinderClosest;
//	}
//};
//
//
//
//
//class SequentialProjectionTarget : IProjectionTarget
//{
//	IProjectionTarget[] Targets{ get; set; }
//
//	SequentialProjectionTarget() { }
//	SequentialProjectionTarget(params IProjectionTarget[] targets)
//	{
//		Targets = targets;
//	}
//
//	Vector3d Project(Vector3d vPoint, int identifier = -1)
//	{
//		Vector3d vCur = vPoint;
//		for (int i = 0; i < Targets.Length; ++i) {
//			vCur = Targets[i].Project(vCur, identifier);
//		}
//		return vCur;
//	}
//};


}