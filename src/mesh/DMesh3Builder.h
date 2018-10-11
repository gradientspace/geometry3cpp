#pragma once

#include <MaterialTypes.h>
#include <DMesh3.h>

namespace g3
{



class MeshMetadata
{
public:
	virtual ~MeshMetadata() {}
};



class DMesh3Builder
{
public:
	enum class AddTriangleFailBehaviors
	{
		DiscardTriangle = 0,
		DuplicateAllVertices = 1
	};

    /// <summary>
    /// What should we do when AddTriangle() fails because triangle is non-manifold?
    /// </summary>
    AddTriangleFailBehaviors NonManifoldTriBehavior = AddTriangleFailBehaviors::DuplicateAllVertices;

    /// <summary>
    /// What should we do when AddTriangle() fails because the triangle already exists?
    /// </summary>
    AddTriangleFailBehaviors DuplicateTriBehavior = AddTriangleFailBehaviors::DiscardTriangle;

	using PDMesh3 = std::shared_ptr<DMesh3>;
	using PMaterial = std::shared_ptr<GenericMaterial>;
	using PMetadata = std::shared_ptr<MeshMetadata>;

    std::vector<PDMesh3> Meshes;
	std::vector<PMaterial> Materials;

    // this is a map from index into Meshes to index into Materials (-1 if no material)
    //  (so, currently we can only have 1 material per mesh!)
    std::vector<int> MaterialAssignment;

	// [RMS] std::any is C++17 feature... 
	std::vector<std::map<std::string, PMetadata>> Metadata;

    int nActiveMesh;

    DMesh3Builder()
    {
        nActiveMesh = -1;
    }

    int AppendNewMesh(bool bHaveVtxNormals, bool bHaveVtxColors, bool bHaveVtxUVs, bool bHaveFaceGroups)
    {
        int index = (int)Meshes.size();
		PDMesh3 m = std::make_unique<DMesh3>(bHaveVtxNormals, bHaveVtxColors, bHaveVtxUVs, bHaveFaceGroups);
        //DMesh3 m = DMesh3(bHaveVtxNormals, bHaveVtxColors, bHaveVtxUVs, bHaveFaceGroups);
        Meshes.push_back(m);
        MaterialAssignment.push_back(-1);     // no material is known
		Metadata.push_back(std::map<std::string, PMetadata>());
        nActiveMesh = index;
        return index;
    }

	int AppendNewMesh(PDMesh3 mesh)
	{
		int index = (int)Meshes.size();
		Meshes.push_back(mesh);
		MaterialAssignment.push_back(-1);     // no material is known
		Metadata.push_back(std::map<std::string, PMetadata>());
		nActiveMesh = index;
		return index;
	}

    int AppendNewMesh(const DMesh3 & existingMesh)
    {
        int index = (int)Meshes.size();
        Meshes.push_back(std::make_unique<DMesh3>(existingMesh));
        MaterialAssignment.push_back(-1);     // no material is known
		Metadata.push_back(std::map<std::string, PMetadata>());
        nActiveMesh = index;
        return index;
    }

    void SetActiveMesh(int id)
    {
        if (id >= 0 && id < Meshes.size())
            nActiveMesh = id;
        else
            throw std::exception("active mesh id is out of range");
    }

    int AppendTriangle(int i, int j, int k)
    {
        return AppendTriangle(i, j, k, -1);
    }

    int AppendTriangle(int i, int j, int k, int g)
    {
        // [RMS] What to do here? We definitely do not want to add a duplicate triangle!!
        //   But is silently ignoring the right thing to do?
        int existing_tid = Meshes[nActiveMesh]->FindTriangle(i, j, k);
        if (existing_tid != InvalidID) {
            if (DuplicateTriBehavior == AddTriangleFailBehaviors::DuplicateAllVertices)
                return append_duplicate_triangle(i, j, k, g);
            else
                return existing_tid;
        }

        int tid = Meshes[nActiveMesh]->AppendTriangle(i, j, k, g);
        if ( tid == NonManifoldID ) {
            if (NonManifoldTriBehavior == AddTriangleFailBehaviors::DuplicateAllVertices)
                return append_duplicate_triangle(i, j, k, g);
            else
                return NonManifoldID;
        }
        return tid;
    }
    int append_duplicate_triangle(int i, int j, int k, int g)
    {
        NewVertexInfo vinfo = NewVertexInfo();
        Meshes[nActiveMesh]->GetVertex(i, vinfo, true, true, true);
        int new_i = Meshes[nActiveMesh]->AppendVertex(vinfo);
        Meshes[nActiveMesh]->GetVertex(j, vinfo, true, true, true);
        int new_j = Meshes[nActiveMesh]->AppendVertex(vinfo);
        Meshes[nActiveMesh]->GetVertex(k, vinfo, true, true, true);
        int new_k = Meshes[nActiveMesh]->AppendVertex(vinfo);
        return Meshes[nActiveMesh]->AppendTriangle(new_i, new_j, new_k, g);
    }



    int AppendVertex(double x, double y, double z)
    {
        return Meshes[nActiveMesh]->AppendVertex(Vector3d(x, y, z));
    }
    int AppendVertex(NewVertexInfo info)
    {
        return Meshes[nActiveMesh]->AppendVertex(info);
    }

    bool SupportsMetaData() {return true; }
    void AppendMetaData(const std::string & identifier, PMetadata data)
    {
		Metadata[nActiveMesh][identifier] = data;
    }


    // just store GenericMaterial object, we can't use it here
    int BuildMaterial(PMaterial m)
    {
        int id = (int)Materials.size();
        Materials.push_back(m);
        return id;
    }

    // do material assignment to mesh
    void AssignMaterial(int materialID, int meshID)
    {
        if (meshID >= MaterialAssignment.size() || materialID >= Materials.size())
            throw std::exception("[DMesh3Builder::AssignMaterial] meshID or materialID are out-of-range");
        MaterialAssignment[meshID] = materialID;
    }





    ////
    //// DMesh3 construction utilities
    ////

    ///// <summary>
    ///// ultimate generic mesh-builder, pass it arrays of floats/doubles, or std::vectors
    ///// of Vector3d, or anything in-between. Will figure out how to interpret
    ///// </summary>
    //static DMesh3 Build<VType,TType,NType>(IEnumerable<VType> Vertices,  
    //                                                IEnumerable<TType> Triangles, 
    //                                                IEnumerable<NType> Normals = null,
    //                                                IEnumerable<int> TriGroups = null)
    //{
    //    DMesh3 mesh = DMesh3(Normals != null, false, false, TriGroups != null);

    //    Vector3d[] v = BufferUtil.ToVector3d(Vertices);
    //    for (int i = 0; i < v.Length; ++i)
    //        mesh.AppendVertex(v[i]);

    //    if ( Normals != null ) {
    //        Vector3f[] n = BufferUtil.ToVector3f(Normals);
    //        if ( n.Length != v.Length )
    //            throw Exception("DMesh3Builder.Build: incorrect number of normals provided");
    //        for (int i = 0; i < n.Length; ++i)
    //            mesh.SetVertexNormal(i, n[i]);
    //    }

    //    Index3i[] t = BufferUtil.ToIndex3i(Triangles);
    //    for (int i = 0; i < t.Length; ++i)
    //        mesh.AppendTriangle(t[i]);

    //    if ( TriGroups != null ) {
    //        std::vector<int> groups = std::vector<int>(TriGroups);
    //        if (groups.size() != t.Length)
    //            throw new Exception("DMesh3Builder.Build: incorect number of triangle groups");
    //        for (int i = 0; i < t.Length; ++i)
    //            mesh.SetTriangleGroup(i, groups[i]);
    //    }

    //    return mesh;
    //}

         


};





}
