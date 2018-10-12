#include <geometry3PCH.h>

#include <DMesh3.h>
#include <DMeshAABBTree3.h>
#include <DMesh3Builder.h>
#include <OBJReader.h>
#include <OBJWriter.h>
#include <MeshConstraints.h>
#include <MeshRefinerBase.h>
#include <Remesher.h>

using namespace g3;

// [RMS] this function just instantiates many of the classes above, which are
// header-only, templates, etc. This helps us find compile errors.

static void test_mesh_classes()
{
	DMesh3Ptr pMesh = std::make_shared<DMesh3>();
	DMeshAABBTree3 test(pMesh, true);

	DMesh3Builder builder;

	OBJReader reader;
	OBJWriter writer;

	MeshConstraints mc;
	MeshRefinerBase refbase;

	Remesher remesher(std::make_shared<DMesh3>());
}