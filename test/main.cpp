
//#include <cinterface.h>

//#include <CTriangulation.h>
//#include <CMeshIO.h>

//#include <cinterface_internal.h>

//#include <vector>
//#include <DMesh3.h>

#include <iostream>

#include <VectorUtil.h>
#include <refcount_vector.h>
#include <small_list_set.h>
#include <DMesh3.h>
using namespace g3;

int main(int argc, char ** argv) 
{
	refcount_vector rcvec1;
	small_list_set smalls;
	
	smalls.Resize(100);

	smalls.AllocateAt(7);
	//for (int k = 0; k <= 12; ++k)
	//	smalls.Insert(7, k);

	for (int value : smalls.values(7, [](int a) { return a + 10; }) ) {
		std::cout << "value: " << value << std::endl;
	}

	DMesh3 mesh;
	mesh.AppendVertex(Vector3d(0, 0, 0));
	mesh.AppendVertex(Vector3d(1, 0, 0));
	mesh.AppendVertex(Vector3d(0, 1, 0));
	mesh.AppendVertex(Vector3d(1, 1, 0));
	int t0 = mesh.AppendTriangle(Vector3i(0, 1, 2));
	int t1 = mesh.AppendTriangle(Vector3i(0, 3, 1));

	bool bValid = mesh.CheckValidity();
	mesh.CompactInPlace();
	bValid = mesh.CheckValidity();

	for (Vector3d v : mesh.Vertices()) {
		std::cout << v[0] << " " << v[1] << " " << v[2] << std::endl;
	}

	for (Index3i t : mesh.Triangles()) {
		std::cout << t[0] << " " << t[1] << " " << t[2] << std::endl;
	}

	for (int eid : mesh.BoundaryEdgeIndices()) {
		Index4i ev = mesh.GetEdge(eid);
		if (ev[3] != InvalidID)
			std::cout << "BAD BOUNDARY EDGE ITR!!";
	}

	for (int vid : mesh.VertexIndices()) {
		std::cout << "vertex " << vid << " tris: ";
		for (int tid : mesh.VtxTrianglesItr(vid))
			std::cout << tid << " ";
		std::cout << std::endl;
	}


	Vector3d axis(0,0,1);
	Matrix3d matrix;
	g3::ComputeAlignZAxisMatrix(axis, matrix, false);

	Vector3d r1 = matrix.row(0);
	Vector3d r2 = matrix.row(1);
	Vector3d r3 = matrix.row(2);

	//std::vector<int> temp; temp.push_back(7); temp.push_back(5);
	//CMeshHandle garbage = (void *)&temp;
	//int valid = is_mesh(garbage);

	//std::vector<double> poly = { 0,0,  1,0,   1,1,   0,1 };
	//CMeshHandle handle = triangulate_polygon(4, &poly[0]);

	//int nv = get_vertex_count(handle);
	//int nt = get_triangle_count(handle);

	//std::vector<double> verts(nv*3);
	//get_vertices_double(handle, nv*3, &verts[0]);

	//std::vector<int> tris(nt*3);
	//get_triangles(handle, nt*3, (unsigned int *)&tris[0]);

}
