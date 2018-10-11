#pragma once

#include <MeshIO.h>
#include <DMesh3Builder.h>
#include <parse_util.h>
#include <file_util.h>
#include <string_util.h>

#include <iostream>
#include <iomanip>

namespace g3 {


class OBJWriter {
public:

    // stream-opener. Override to write to something other than a file.
	//std::function<std::ostream(std::string)> OpenStreamF = [](std::string filename) {
	//	return std::ofstream(filename);
 //   };
	//std::function<void(std::ostream&)> CloseStreamF = [](std::ostream & stream) {
	//	// do we need to explicitly close ofstream &&
 //   };

    std::string GroupNamePrefix = "mmGroup";   // default, compatible w/ meshmixer
    std::function<std::string(int)> GroupNameF = nullptr;  // use this to replace standard group names w/ your own

    //IOWriteResult Write(BinaryWriter writer, List<WriteMesh> vMeshes, WriteOptions options)
    //{
    //    // [RMS] not supported
    //    throw new NotImplementedException();
    //}

    IOWriteResult Write(std::ostream & stream, const std::vector<WriteMesh> & vMeshes, const WriteOptions & options)
    {
        if (! options.groupNamePrefix.empty())
            GroupNamePrefix = options.groupNamePrefix;
        if (options.GroupNameF != nullptr)
            GroupNameF = options.GroupNameF;

        int nAccumCountV = 1;       // OBJ indices always start at 1
        int nAccumCountUV = 1;

		// collect materials
		std::string sMaterialLib = "";
		int nHaveMaterials = 0;
		if ( options.bWriteMaterials && options.MaterialFilePath.size() > 0 ) {
			// [RMS] disabled for now
			//List<GenericMaterial> vMaterials = MeshIOUtil.FindUniqueMaterialList(vMeshes);
			//IOWriteResult ok = write_materials(vMaterials, options);
			//if ( ok.code == IOCode.Ok ) {
			//	sMaterialLib = Path.GetFileName(options.MaterialFilePath);
			//	nHaveMaterials = vMeshes.Count;
			//}
		}


		if (options.AsciiHeaderFunc != nullptr)
			stream << options.AsciiHeaderFunc() << std::endl;

		if (! sMaterialLib.empty())
			stream << "mtllib " << sMaterialLib << std::endl;

        for (int mi = 0; mi < vMeshes.size(); ++mi) {
            DMesh3 & mesh = * vMeshes[mi].Mesh;

            if (options.ProgressFunc != nullptr)
                options.ProgressFunc(mi, (int)vMeshes.size());

            bool bVtxColors = options.bPerVertexColors && mesh.HasVertexColors();
            bool bNormals = options.bPerVertexNormals && mesh.HasVertexNormals();

			// use separate UV set if we have it, otherwise write per-vertex UVs if we have those
			bool bVtxUVs = options.bPerVertexUVs && mesh.HasVertexUVs();
			if ( vMeshes[mi].UVs != nullptr)
				bVtxUVs = false;

			std::vector<int> mapV;
			mapV.resize(mesh.MaxVertexID());

			int digits = options.RealPrecisionDigits;
			stream << std::setprecision(digits);

			// write vertices for this mesh
            for ( int vi : mesh.VertexIndices() ) { 
				mapV[vi] = nAccumCountV++;
                Vector3d v = mesh.GetVertex(vi);
                if ( bVtxColors ) {
                    Vector3f c = mesh.GetVertexColor(vi);
					stream << "v " 
						<< std::setprecision(digits) << v[0] << " " << v[1] << " " << v[2] << " "
						<< std::setprecision(8) << c[0] << " " << c[1] << " " << c[2] << std::endl;
				} else {
					stream << "v " << v[0] << " " << v[1] << " " << v[2] << std::endl;
				}

				if ( bNormals ) {
                    Vector3f n = mesh.GetVertexNormal(vi);
					stream << "vn " << n[0] << " " << n[1] << " " << n[2] << std::endl;
                }

				if ( bVtxUVs ) {
					Vector2f uv = mesh.GetVertexUV(vi);
					stream << "vt " << uv[0] << " " << uv[1] << std::endl;
				}
            }

			// [RMS] disabled...

            // write independent UVs for this mesh, if we have them
			//IIndexMap mapUV = (bVtxUVs) ? new IdentityIndexMap() : nullptr;
   //         DenseUVMesh uvSet = nullptr;
   //         if ( vMeshes[mi].UVs != nullptr) {
   //             uvSet = vMeshes[mi].UVs;
   //             int nUV = uvSet.UVs.Length;
			//	IndexMap fullMap = new IndexMap(false, nUV);   // [TODO] do we really need a map here? is just integer shift, no?
   //             for (int ui = 0; ui < nUV; ++ui) {
   //                 writer.WriteLine("vt {0:F8} {1:F8}", uvSet.UVs[ui].x, uvSet.UVs[ui].y);
			//		fullMap[ui] = nAccumCountUV++;
   //             }
			//	mapUV = fullMap;
   //         }

			// check if we need to write usemtl lines for this mesh
			//bool bWriteMaterials = nHaveMaterials > 0 
   //                 && vMeshes[mi].TriToMaterialMap != nullptr
   //                 && vMeshes[mi].Materials != nullptr;
			bool bWriteMaterials = false;

			DenseUVMesh * uvSet = nullptr;
			std::map<int, int> * mapUV = nullptr;

			// various ways we can write triangles to minimize state changes...
			// [TODO] support writing materials when mesh has groups!!
            if (options.bWriteGroups && mesh.HasTriangleGroups())
                write_triangles_bygroup(stream, mesh, mapV, uvSet, mapUV, bNormals);
            else
				write_triangles_flat(stream, vMeshes[mi], mapV, uvSet, mapUV, bNormals, bWriteMaterials);

            if (options.ProgressFunc != nullptr)
                options.ProgressFunc(mi+1, (int)vMeshes.size());
        }


		return IOWriteResult::Ok();
    }



	// write triangles of mesh with re-ordering to minimize group changes
	// (note: this may mean lots of material changes, depending on mesh...)
	void write_triangles_bygroup(std::ostream & stream, 
		DMesh3 & mesh, std::vector<int> & mapV, 
		DenseUVMesh * uvSet, std::map<int, int> * mapUV,
		bool bNormals)
    {
        // This makes N passes over mesh indices, but doesn't use much extra memory.
        // would there be a faster way? could construct integer-pointer-list during initial
        // scan, this would need O(N) memory but then write is effectively O(N) instead of O(N*k)

		// TODO construct triangle-list pairs <tid,gid> and then sort by gid!

        bool bUVs = (mapUV != nullptr);

		std::set<int> vGroups;		// std::set is sorted!
        for (int ti : mesh.TriangleIndices())
            vGroups.insert(mesh.GetTriangleGroup(ti));

        for ( int g : vGroups ) {
            std::string group_name = GroupNamePrefix;
            if (GroupNameF != nullptr ) {
                group_name = GroupNameF(g);
            } else {
				group_name = GroupNamePrefix + std::to_string(g);
            }
			stream << "g " << group_name << std::endl;

            for (int ti : mesh.TriangleIndices() ) {
                if (mesh.GetTriangleGroup(ti) != g)
                    continue;

                Index3i t = mesh.GetTriangle(ti);
				t[0] = mapV[t[0]];
				t[1] = mapV[t[1]];
				t[2] = mapV[t[2]];

                if (bUVs) {
					Index3i tuv = (uvSet != nullptr) ? uvSet->TriangleUVs[ti] : t;
					tuv[0] = (*mapUV)[tuv[0]];
					tuv[1] = (*mapUV)[tuv[1]];
					tuv[2] = (*mapUV)[tuv[2]];
					write_tri(stream, t, bNormals, true, tuv);
                } else {
					write_tri(stream, t, bNormals, false, t);
                }

            }
        }
    }


	// sequential write of input mesh triangles. preserves triangle IDs up to constant shift.
	void write_triangles_flat(std::ostream & stream, 
		const WriteMesh & write_mesh, std::vector<int> & mapV, 
		DenseUVMesh * uvSet, std::map<int,int> * mapUV, 
		bool bNormals, bool bMaterials)
    {
        bool bUVs = (mapUV != nullptr);

		int cur_material = -1;

		const DMesh3 & mesh = * write_mesh.Mesh;
        for (int ti : mesh.TriangleIndices() ) { 
			if ( bMaterials )
				set_current_material(stream, ti, write_mesh, cur_material);

            Index3i t = mesh.GetTriangle(ti);
			t[0] = mapV[t[0]];
			t[1] = mapV[t[1]];
			t[2] = mapV[t[2]];

            if (bUVs) {
				Index3i tuv = (uvSet != nullptr) ? uvSet->TriangleUVs[ti] : t;
                tuv[0] = (*mapUV)[tuv[0]];
                tuv[1] = (*mapUV)[tuv[1]];
                tuv[2] = (*mapUV)[tuv[2]];
                write_tri(stream, t, bNormals, true, tuv);
            } else {
                write_tri(stream, t, bNormals, false, t);
            }
        }
    }

	// update material state if necessary
	void set_current_material(std::ostream & stream, int ti, const WriteMesh & mesh, int & cur_material) 
	{
		int mi = mesh.TriToMaterialMap.find(ti)->second;  // because of const!
		if ( mi != cur_material && mi >= 0 && mi < mesh.Materials.size() ) {
			stream << "usemtl " << mesh.Materials[mi]->name << std::endl;
			cur_material = mi;
		}
	}


	// actually write triangle line, in proper OBJ format
    void write_tri(std::ostream & stream, const Index3i & t, bool bNormals, bool bUVs, const Index3i & tuv)
    {
        if ( bNormals == false && bUVs == false ) {
			stream << "f " << t[0] << " " << t[1] << " " << t[2] << std::endl;
        } else if ( bNormals == true && bUVs == false ) {
			stream << "f " << t[0] << "//" << t[0] << " "
							<< t[1] << "//" << t[1] << " "
							<< t[2] << "//" << t[2] << std::endl;
        } else if ( bNormals == false && bUVs == true ) {
			stream << "f " << t[0] << "/" << tuv[0] << " "
							<< t[1] << "/" << tuv[1] << " "
							<< t[2] << "/" << tuv[2] << std::endl;
        } else {
			stream << "f " << t[0] << "/" << tuv[0] << "/" << t[0] << " "
							<< t[1] << "/" << tuv[1] << "/" << t[1] << " "
							<< t[2] << "/" << tuv[2] << "/" << t[2] << std::endl;
        }
    }



	// [RMS] disabled for now, not critical

	// write .mtl file
/*
	IOWriteResult write_materials(List<GenericMaterial> vMaterials, WriteOptions options) 
	{
        Stream stream = OpenStreamF(options.MaterialFilePath);
        if (stream == nullptr)
            return new IOWriteResult(IOCode.FileAccessError, "Could not open file " + options.MaterialFilePath + " for writing");


        try { 
            StreamWriter w = new StreamWriter(stream);

			foreach ( GenericMaterial gmat in vMaterials ) {
				if ( gmat is OBJMaterial == false )
					continue;
				OBJMaterial mat = gmat as OBJMaterial;

				w.WriteLine("newmtl {0}", mat.name);
				if ( mat.Ka != GenericMaterial.Invalid )
					w.WriteLine("Ka {0} {1} {2}", mat.Ka.x, mat.Ka.y, mat.Ka.z);
				if ( mat.Kd != GenericMaterial.Invalid)
					w.WriteLine("Kd {0} {1} {2}", mat.Kd.x, mat.Kd.y, mat.Kd.z);
				if ( mat.Ks != GenericMaterial.Invalid )
					w.WriteLine("Ks {0} {1} {2}", mat.Ks.x, mat.Ks.y, mat.Ks.z);
				if ( mat.Ke != GenericMaterial.Invalid )
					w.WriteLine("Ke {0} {1} {2}", mat.Ke.x, mat.Ke.y, mat.Ke.z);
				if ( mat.Tf != GenericMaterial.Invalid )
					w.WriteLine("Tf {0} {1} {2}", mat.Tf.x, mat.Tf.y, mat.Tf.z);
				if ( mat.d != Single.MinValue )
					w.WriteLine("d {0}", mat.d);
				if ( mat.Ns != Single.MinValue )
					w.WriteLine("Ns {0}", mat.Ns);
				if ( mat.Ni != Single.MinValue )
					w.WriteLine("Ni {0}", mat.Ni);
				if ( mat.sharpness != Single.MinValue )
					w.WriteLine("sharpness {0}", mat.sharpness);
				if ( mat.illum != -1 )
					w.WriteLine("illum {0}", mat.illum);

				if ( mat.map_Ka != null && mat.map_Ka != "" )
					w.WriteLine("map_Ka {0}", mat.map_Ka);
				if ( mat.map_Kd != null && mat.map_Kd != "" )
					w.WriteLine("map_Kd {0}", mat.map_Kd);
				if ( mat.map_Ks != null && mat.map_Ks != "" )
					w.WriteLine("map_Ks {0}", mat.map_Ks);
				if ( mat.map_Ke != null && mat.map_Ke != "" )
					w.WriteLine("map_Ke {0}", mat.map_Ke);
				if ( mat.map_d != null && mat.map_d != "" )
					w.WriteLine("map_d {0}", mat.map_d);
				if ( mat.map_Ns != null && mat.map_Ns != "" )
					w.WriteLine("map_Ns {0}", mat.map_Ns);
				
				if ( mat.bump != null && mat.bump != "" )
					w.WriteLine("bump {0}", mat.bump);
				if ( mat.disp != null && mat.disp != "" )
					w.WriteLine("disp {0}", mat.disp);				
				if ( mat.decal != null && mat.decal != "" )
					w.WriteLine("decal {0}", mat.decal);
				if ( mat.refl != null && mat.refl != "" )
					w.WriteLine("refl {0}", mat.refl);
			}

        } finally {
            CloseStreamF(stream);
        }

		return IOWriteResult.Ok;
	}
*/

};


}