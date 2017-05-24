#include <geometry3PCH.h>
#include <CMeshIO.h>
#include "cinterface_internal.h"


CMeshHandle create_packed_mesh() {
	return new MeshHandle(new PackedMesh);
}

CMeshHandle create_dynamic_mesh() {
	return new MeshHandle(new DMesh3<double>());
}

void release_mesh( CMeshHandle handle ) {
	if (is_mesh( handle )) {
		MeshHandle * p = (MeshHandle *)handle;
		p->Release();
		delete p;
	}
}


int is_mesh( CMeshHandle handle ) {
	MeshHandle * h = (MeshHandle *)handle;
	if ( h->valid_token != MeshHandle::VALID )
		return 0;
	return ((MeshHandle *)handle)->is_valid();
}

int get_vertex_count( CMeshHandle handle ) {
	return is_mesh(handle) ? ((MeshHandle*)handle)->vertex_count() : -1;
}
int get_triangle_count(CMeshHandle handle) {
	return is_mesh(handle) ? ((MeshHandle*)handle)->triangle_count() : -1;
}

g3ExternalC int get_vertices_float(CMeshHandle handle, int buffer_size, float * buffer ) {
	if ( is_mesh(handle) == 0)
		return -1;
	return ((MeshHandle*)handle)->get_vertices_float(buffer_size, buffer);
}

g3ExternalC int get_vertices_double(CMeshHandle handle, int buffer_size, double * buffer ) {
	if ( is_mesh(handle) == 0)
		return -1;
	return ((MeshHandle*)handle)->get_vertices_double(buffer_size, buffer);
}

g3ExternalC int get_triangles(CMeshHandle handle, int buffer_size, unsigned int * buffer ) {
	if ( is_mesh(handle) == 0)
		return -1;
	return ((MeshHandle*)handle)->get_triangles(buffer_size, buffer);
}