#include <geometry3PCH.h>
#include <CTriangulation.h>
#include "cinterface_internal.h"

#include <GteTriangulateEC.h>
#include <GteBSRational.h>
#include <GteUIntegerAP32.h>

CMeshHandle triangulate_polygon( int nVertices, double * vertexBuffer )
{
	typedef gte::BSRational<gte::UIntegerAP32> gtnumber;
	typedef gte::Vector2<double> Vector2d;

	std::vector<Vector2d> vertices(nVertices);
	for ( int i = 0; i < nVertices; ++i ) 
		vertices[i] = { vertexBuffer[2*i], vertexBuffer[2*i+1] };

	gte::TriangulateEC<double, gtnumber> t(vertices);
	bool bOK = t();
	if ( bOK == false )
		return nullptr;

	PackedMesh * pMesh = new PackedMesh();
	pMesh->ResizeVertices(nVertices, false, false);
	for (int i = 0; i < nVertices; ++i) {
		pMesh->SetPosition(i, (float)vertexBuffer[2*i], (float)vertexBuffer[2*i+1], 0);
	}

	auto triangles = t.GetTriangles();
	unsigned int NT = (unsigned int)triangles.size();
	pMesh->ResizeTriangles(NT);
	for ( unsigned int i = 0; i < NT; ++i ) {
		auto t = triangles[i];
		pMesh->SetTriangle(i, t[0],t[1],t[2]);
	}

	return new MeshHandle(pMesh);
}