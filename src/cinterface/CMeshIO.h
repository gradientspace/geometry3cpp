#pragma once


#include "cinterface.h"

#ifdef __cplusplus
extern "C"
{
#endif


// type = 0: PackedMesh 1: DMesh3
g3ExternalC CMeshHandle create_packed_mesh();
g3ExternalC CMeshHandle create_dynamic_mesh();
g3ExternalC void release_mesh(CMeshHandle handle);

g3ExternalC int is_mesh(CMeshHandle handle);
g3ExternalC int get_vertex_count(CMeshHandle handle);
g3ExternalC int get_triangle_count(CMeshHandle handle);

g3ExternalC int get_vertices_float(CMeshHandle handle, int buffer_size, float * buffer );
g3ExternalC int get_vertices_double(CMeshHandle handle, int buffer_size, double * buffer );
g3ExternalC int get_triangles(CMeshHandle handle, int buffer_size, unsigned int * buffer );

#ifdef __cplusplus
}
#endif
