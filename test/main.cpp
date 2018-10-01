
//#include <cinterface.h>

//#include <CTriangulation.h>
//#include <CMeshIO.h>

//#include <cinterface_internal.h>

//#include <vector>
//#include <DMesh3.h>

#include <VectorUtil.h>
using namespace g3;

int main(int argc, char ** argv) 
{
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
