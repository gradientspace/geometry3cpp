#pragma once

#define g3ExternalC __declspec(dllexport)

#ifdef __cplusplus
extern "C"
{
#endif

typedef void * CMeshHandle;

// type = 0: PackedMesh 1: DMesh3
g3ExternalC CMeshHandle create_packed_mesh();
g3ExternalC CMeshHandle create_dynamic_mesh();
g3ExternalC void release_mesh(CMeshHandle handle);


#ifdef __cplusplus
}
#endif