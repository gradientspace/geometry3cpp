#include <geometry3PCH.h>
#include <CMeshIO.h>

#include <g3types.h>
#include <PackedMesh.h>
#include <DMesh3.h>

using namespace g3;

struct MeshHandle
{
	PackedMesh * packed;
	DMesh3<double> * dynamic;
	int valid_token;
};


CMeshHandle create_packed_mesh()
{
	MeshHandle * p = new MeshHandle();
	p->packed = new PackedMesh();
	p->dynamic = nullptr;
	p->valid_token = 0xdeadbeef;
	return p;
}

CMeshHandle create_dynamic_mesh()
{
	MeshHandle * p = new MeshHandle();
	p->packed = nullptr;
	p->dynamic = new DMesh3<double>();
	p->valid_token = 0xdeadbeef;
	return p;
}


void release_mesh( CMeshHandle handle )
{
	MeshHandle * p = (MeshHandle * )handle;
	if ( p->valid_token != 0xdeadbeef )
		return;
	p->valid_token = 0xf00ff00f;
	if ( p->packed != nullptr )
		delete p->packed;
	if ( p->dynamic != nullptr )
		delete p->dynamic;
	delete p;
}