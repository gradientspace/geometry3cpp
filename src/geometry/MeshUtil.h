#pragma once

#include <g3types.h>
#include <GeometryInterfaces.h>
#include <DMesh3.h>


namespace g3
{

	double Area( const IPackedMesh * pMesh );
	double Volume( const IPackedMesh * pMesh );

	Box3f AxisBoundingBox( const IPackedMesh * pMesh );
	Box3f OrientedBoundingBox( const IPackedMesh * pMesh, bool bFast = true );

	void Translate( IPackedMesh * pMesh, const Vector3f & vTranslation );
	void Scale( IPackedMesh * pMesh, const Vector3f & vScale);
	inline void Scale( IPackedMesh * pMesh, float fScale) 
		{ Scale(pMesh, Vector3f(fScale,fScale,fScale)); }



inline Vector3d UniformSmooth(const DMesh3 & mesh, int vID, double t)
{
	Vector3d v = mesh.GetVertex(vID);
	Vector3d c;
	mesh.VtxOneRingCentroid(vID, c);
	return (1.0-t)*v + (t)*c;
}


}  // namespace g3



