#include <geometry3PCH.h>
#include <MeshSubdivider.h>

using namespace g3;


MeshSubdivider::MeshSubdivider()
{
}


void MeshSubdivider::Split1to4( DMesh3 & mesh )
{
	EdgeID nMaxEdge = mesh.MaxEdgeID();
	std::vector<VertexID> vNewV(nMaxEdge, InvalidID);

	for (auto eid : mesh.EdgeIndices())
		vNewV[eid] = MarkerID1;

	// [TODO] can we do this without saving triangles? Could we process each triangle
	//   sequentially somehow? The problem is that our original triangles go away
	//   each time we split...

	// save triangle edges & verts (do we need both?)
	static const Vector3i InvalidTri(-1,-1,-1);
	TriangleID nMaxInitialTri = mesh.MaxTriangleID();
	std::vector<Vector3i> vTris(nMaxInitialTri);
	std::vector<Vector3i> vTriEdges(nMaxInitialTri);
	for (int k = 0; k < nMaxInitialTri; ++k) {
		if (mesh.IsTriangle( k )) {
			vTris[k] = mesh.GetTriangle( k );
			vTriEdges[k] = mesh.GetTriEdges( k );
		} else {
			vTris[k] = InvalidTri;
			vTriEdges[k] = InvalidTri;
		}
	}

	// split all edges
	for (auto eid : mesh.EdgeIndices()) {
		if (vNewV[eid] == MarkerID1) {
			DMesh3::EdgeSplitInfo info;
			auto eResult = mesh.SplitEdge(eid, info);
			vNewV[eid] = (eResult == MeshResult::Ok) ? info.vNew : InvalidID;
		}
	}

	// for each triangle, if split edge doesn't exist, then it is flipped 
	//   and we need to rotate it one more time. Should happen once per triangle.
	for (int k = 0; k < nMaxInitialTri; ++k) {
		if (vTriEdges[k] == InvalidTri)
			continue;
		for (int j = 0; j < 3; ++j) {
			EdgeID e0 = vTriEdges[k][j];
			EdgeID e1 = vTriEdges[k][(j+1)%3];
			VertexID a = vNewV[e0];
			VertexID b = vNewV[e1];
			if (mesh.FindEdge( a, b ) == InvalidID) {
				VertexID o = vTris[k][(j+1)%3];
				VertexID c = vNewV[vTriEdges[k][(j+2)%3]];
				if ( c != InvalidID && mesh.FindEdge( o, c ) ) {
					DMesh3::EdgeFlipInfo info;
					auto eResult = mesh.FlipEdge( o, c, info );
				}
			}
		}
	}
}

