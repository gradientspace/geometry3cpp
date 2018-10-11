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

static void test_mesh_classes()
{
	DMesh3 mesh;

	DMesh3Builder builder;

	OBJReader reader;
	OBJWriter writer;

	MeshConstraints mc;
	MeshRefinerBase refbase;

	Remesher remesher(std::make_shared<DMesh3>());
}