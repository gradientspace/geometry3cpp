#pragma once

#include <MeshIO.h>
#include <DMesh3Builder.h>
#include <parse_util.h>
#include <file_util.h>
#include <string_util.h>


namespace g3 {














/// <summary>
/// gradientspace OBJ mesh format parser
/// 
/// Basic structure is:
///   1) parse OBJ into internal data structures that represent OBJ exactly
///   2) convert to mesh objects based on options/etc
/// 
/// [TODO] major current limitation is that we do not support multiple UVs per-vertex
///   (this is a limitation of DMesh3). Similarly no multiple normals per-vertex. So,
///   in the current code, such vertices are duplicated. See append_vertex() and use
///   of Index3i triplet (vi,ni,ui) to represent vertex in BuildMeshes_X() functions
/// 
/// [TODO] only a single material can be assigned to a mesh, this is a current limitation
///   of DMesh3Builder. So, in this case we are splitting the input mesh by material, IE
///   multiple meshes are returned for a single input mesh, each with one material.
/// 
/// 
/// </summary>


struct Triangle
{
    static constexpr int InvalidMaterialID = -1;
	static constexpr int InvalidGroupID = -1;

    Index3i vIndices;
    Index3i vNormals;
    Index3i vUVs;
    int nMaterialID;
    int nGroupID;

    void clear()
    {
        nMaterialID = InvalidMaterialID;
        nGroupID = InvalidGroupID;
        vIndices = vNormals = vUVs = Index3i(-1, -1, -1);
    }

    void set_vertex(int j, int vi, int ni = -1, int ui = -1)
    {
        vIndices[j] = vi;
        if (ni != -1) vNormals[j] = ni;
        if (ui != -1) vUVs[j] = ui;
    }

    void move_vertex(int jFrom, int jTo)
    {
        vIndices[jTo] = vIndices[jFrom];
        vNormals[jTo] = vNormals[jFrom];
        vUVs[jTo] = vUVs[jFrom];
    }

    bool is_complex()
    {
        for ( int j = 0; j < 3; ++j ) {
            if (vNormals[j] != -1 && vNormals[j] != vNormals[j])
                return true;
            if (vUVs[j] != -1 && vUVs[j] != vUVs[j])
                return true;
        }
        return false;
    }
};





class OBJReader
{
public:
    dvector<double> vPositions;
    dvector<float> vNormals;
    dvector<float> vUVs;
    dvector<float> vColors;
    dvector<Triangle> vTriangles;

	std::map<std::string, std::shared_ptr<OBJMaterial>> Materials;
    std::map<int, std::string> UsedMaterials;

    bool m_bOBJHasPerVertexColors;
    int m_nUVComponents;

    bool m_bOBJHasTriangleGroups;
    int m_nSetInvalidGroupsTo;

    int nWarningLevel = 0;      // 0 == no diagnostics, 1 == basic, 2 == crazy
    std::map<std::string, int> warningCount = std::map<std::string, int>();

    OBJReader()
    {
        MTLFileSearchPaths = std::vector<std::string>();
    }

	// you need to initialize this with paths if you want .MTL files to load
	std::vector<std::string> MTLFileSearchPaths;

    // connect to this to get warning messages
	std::function<void(const std::string &)> warningEvent = nullptr;

    bool HasPerVertexColors() { return m_bOBJHasPerVertexColors; }
    int UVDimension() { return m_nUVComponents; }

    bool HasTriangleGroups() { return m_bOBJHasTriangleGroups; }


    // if this is true, means during parsing we found vertices of faces that
    //  had different indices for vtx/normal/uv
	bool HasComplexVertices;


	std::vector<std::string> LINES;


    //IOReadResult Read(BinaryReader reader, ReadOptions options, DMesh3Builder & builder)
    //{
    //    throw NotImplementedException();
    //}

    IOReadResult Read(std::istream & reader, const ReadOptions & options, DMesh3Builder & builder)
    {
        Materials = std::map<std::string, std::shared_ptr<OBJMaterial>>();
        UsedMaterials = std::map<int, std::string>();
        HasComplexVertices = false;

        if (nWarningLevel >= 1)
            emit_warning("[OBJReader] starting parse");

        auto parseResult = ParseInput(reader, options);
        if (parseResult.code != IOCode::Ok)
            return parseResult;

        if (nWarningLevel >= 1)
            emit_warning("[OBJReader] completed parse. building.");

        auto buildResult = 
            (UsedMaterials.size()> 1 || HasComplexVertices) ?
                BuildMeshes_ByMaterial(options, builder) : BuildMeshes_Simple(options, builder);

        if (nWarningLevel >= 1)
            emit_warning("[OBJReader] build complete.");

        if (buildResult.code != IOCode::Ok)
            return buildResult;

        return IOReadResult(IOCode::Ok, "");
    }







	int append_vertex(DMesh3Builder & builder, const Index3i & vertIdx, bool bHaveNormals, bool bHaveColors, bool bHaveUVs )
    {
        int vi = 3 * vertIdx[0];
        if ( vertIdx[0] < 0 || vertIdx[0] >= vPositions.size()/3 ) {
            emit_warning("[OBJReader] append_vertex() referencing invalid vertex " + std::to_string(vertIdx[0]) );
            return -1;
        }

        if ( bHaveNormals == false && bHaveColors == false && bHaveUVs == false )
            return builder.AppendVertex(vPositions[vi], vPositions[vi + 1], vPositions[vi + 2]);

        NewVertexInfo vinfo = NewVertexInfo();
        vinfo.bHaveC = vinfo.bHaveN = vinfo.bHaveUV = false;
        vinfo.v = Vector3d(vPositions[vi], vPositions[vi + 1], vPositions[vi + 2]);
        if ( bHaveNormals ) {
            vinfo.bHaveN = true;
            int ni = 3 * vertIdx[1];
            vinfo.n = Vector3f(vNormals[ni], vNormals[ni + 1], vNormals[ni + 2]);
        }
        if ( bHaveColors ) {
            vinfo.bHaveC = true;
            vinfo.c = Vector3f(vColors[vi], vColors[vi + 1], vColors[vi + 2]);
        }
        if ( bHaveUVs ) {
            vinfo.bHaveUV = true;
            int ui = 2 * vertIdx[2];
            vinfo.uv = Vector2f(vUVs[ui], vUVs[ui + 1]);
        }

        return builder.AppendVertex(vinfo);
    }



    int append_triangle(DMesh3Builder & builder, int nTri, const std::vector<int> mapV)
    {
        Triangle t = vTriangles[nTri];
        int v0 = mapV[t.vIndices[0] - 1];
        int v1 = mapV[t.vIndices[1] - 1];
        int v2 = mapV[t.vIndices[2] - 1];
        if ( v0 == -1 || v1 == -1 || v2 == -1 ) {
            emit_warning( StringUtil::Format("[OBJReader] invalid triangle:  {0} {1} {2}  mapped to {3} {4} {5}",
                t.vIndices[0], t.vIndices[1], t.vIndices[2], v0, v1, v2));
            return -1;
        }
        int gid = (vTriangles[nTri].nGroupID == InvalidGroupID) ? 
            m_nSetInvalidGroupsTo : vTriangles[nTri].nGroupID;
        return builder.AppendTriangle(v0, v1, v2, gid);
    }
    int append_triangle(DMesh3Builder & builder, Triangle t)
    {
        if ( t.vIndices[0] < 0 || t.vIndices[1] < 0 || t.vIndices[2] < 0 ) {
            emit_warning( StringUtil::Format("[OBJReader] invalid triangle:  {0} {1} {2}",
                t.vIndices[0], t.vIndices[1], t.vIndices[2]));
            return -1;
        }
        int gid = (t.nGroupID == InvalidGroupID) ? 
            m_nSetInvalidGroupsTo : t.nGroupID;
        return builder.AppendTriangle(t.vIndices[0], t.vIndices[1], t.vIndices[2], gid);
    }


    IOReadResult BuildMeshes_Simple(const ReadOptions & options, DMesh3Builder & builder)
    {
        if (vPositions.size() == 0)
            return IOReadResult(IOCode::GarbageDataError, "No vertices in file");
        if (vTriangles.size() == 0)
            return IOReadResult(IOCode::GarbageDataError, "No triangles in file");

        // [TODO] support non-per-vertex normals/colors
        bool bHaveNormals = (vNormals.size() == vPositions.size());
        bool bHaveColors = (vColors.size() == vPositions.size());
        bool bHaveUVs = (vUVs.size()/2 == vPositions.size()/3);

        int nVertices = (int)vPositions.size() / 3;
		std::vector<int> mapV(nVertices);

        int meshID = builder.AppendNewMesh(bHaveNormals, bHaveColors, bHaveUVs, m_bOBJHasTriangleGroups);
        for (int k = 0; k < nVertices; ++k) {
			Index3i vk = Index3i(k,k,k);
            mapV[k] = append_vertex(builder, vk, bHaveNormals, bHaveColors, bHaveUVs);
        }

        // [TODO] this doesn't handle missing vertices...
        for (int k = 0; k < vTriangles.size(); ++k)
            append_triangle(builder, k, mapV);

        if ( UsedMaterials.size() == 1 ) {       // [RMS] should not be in here otherwise
			int material_id = UsedMaterials.begin()->first;
            std::string sMatName = UsedMaterials[material_id];
            auto useMat = Materials[sMatName];
            int matID = builder.BuildMaterial( std::dynamic_pointer_cast<GenericMaterial>(useMat) );
            builder.AssignMaterial(matID, meshID);
        }

        return IOReadResult(IOCode::Ok, "");
    }




	IOReadResult BuildMeshes_ByMaterial(ReadOptions options, DMesh3Builder & builder)
	{
		throw std::exception("TODO");
	}
   // IOReadResult BuildMeshes_ByMaterial(ReadOptions options, DMesh3Builder & builder)
   // {
   //     if (vPositions.size() == 0)
   //         return IOReadResult(IOCode.GarbageDataError, "No vertices in file");
   //     if (vTriangles.size() == 0)
   //         return IOReadResult(IOCode.GarbageDataError, "No triangles in file");

   //     bool bHaveNormals = (vNormals.size() > 0);
   //     bool bHaveColors = (vColors.size() > 0);
   //     bool bHaveUVs = (vUVs.size() > 0);

   //     std::vector<int> usedMaterialIDs = std::vector<int>(UsedMaterials.Keys);
   //     usedMaterialIDs.push_back(Triangle.InvalidMaterialID);
   //     foreach ( int material_id in usedMaterialIDs) {
   //         int matID = Triangle.InvalidMaterialID;
   //         if (material_id != Triangle.InvalidMaterialID) {
   //             string sMatName = UsedMaterials[material_id];
   //             OBJMaterial useMat = Materials[sMatName];
   //             matID = builder.BuildMaterial(useMat);
   //         }
   //         bool bMatHaveUVs = (material_id == Triangle.InvalidMaterialID) ? false : bHaveUVs;

			//// don't append mesh until we actually see triangles
			//int meshID = -1;

   //         std::map<Index3i, int> mapV = std::map<Index3i, int>();

   //         for ( int k = 0; k < vTriangles.size(); ++k ) {

   //             Triangle t = vTriangles[k];
   //             if (t.nMaterialID == material_id) {

			//		if ( meshID == -1 )
			//			meshID = builder.AppendNewMesh(bHaveNormals, bHaveColors, bMatHaveUVs, false);
			//			
   //                 Triangle t2 = Triangle();
   //                 for (int j = 0; j < 3; ++j) {
   //                     Index3i vk = Index3i(
			//				t.vIndices[j] - 1, t.vNormals[j] - 1, t.vUVs[j] - 1 );

   //                     int use_vtx = -1;
   //                     if (mapV.ContainsKey(vk) == false) {
   //                         use_vtx = append_vertex(builder, vk, bHaveNormals, bHaveColors, bMatHaveUVs);
   //                         mapV[vk] = use_vtx;
   //                     } else
   //                         use_vtx = mapV[vk];

   //                     t2.vIndices[j] = use_vtx;
   //                 }
   //                 append_triangle(builder, t2);
   //             }
   //         }

   //         if ( matID != Triangle.InvalidMaterialID )
   //             builder.AssignMaterial(matID, meshID);
   //     }

   //     return IOReadResult(IOCode.Ok, "");
   // }





    IOReadResult ParseInput(std::istream & reader, ReadOptions options)
    {
        vPositions = dvector<double>();
        vNormals = dvector<float>();
        vUVs = dvector<float>();
        vColors = dvector<float>();
        vTriangles = dvector<Triangle>();

        bool bVerticesHaveColors = false;
        int nMaxUVLength = 0;
        std::shared_ptr<OBJMaterial> activeMaterial = nullptr;

        std::map<std::string, int> GroupNames = std::map<std::string, int>();
        int nGroupCounter = 0;
        int nActiveGroup = InvalidGroupID;

		std::vector<std::string> tokens, partsBuf;
		int nLines = 0;
		std::string line;
		while (std::getline(reader, line)) {
		//for ( int k = 0; k < LINES.size(); ++k ) {
		//	std::string & line = LINES[k];
			
            nLines++;
			ParseUtil::Split(line, ' ', tokens, ParseUtil::Options::RemoveEmpty);
            if (tokens.size() == 0)
                continue;

            // [RMS] this will hang VS on large models...
            //if (nWarningLevel >= 2)
            //    emit_warning("Parsing line " + line);
            try {

                if (tokens[0][0] == 'v') {
                    if (tokens[0].size() == 1) {
                        if (tokens.size() == 7) {
                            vPositions.push_back(ParseUtil::ToDouble(tokens[1]));
                            vPositions.push_back(ParseUtil::ToDouble(tokens[2]));
                            vPositions.push_back(ParseUtil::ToDouble(tokens[3]));

                            vColors.push_back(ParseUtil::ToFloat(tokens[4]));
                            vColors.push_back(ParseUtil::ToFloat(tokens[5]));
                            vColors.push_back(ParseUtil::ToFloat(tokens[6]));
                            bVerticesHaveColors = true;
                        } else if (tokens.size() >= 4) {
                            vPositions.push_back(ParseUtil::ToDouble(tokens[1]));
                            vPositions.push_back(ParseUtil::ToDouble(tokens[2]));
                            vPositions.push_back(ParseUtil::ToDouble(tokens[3]));

                        } 
                        if ( tokens.size() != 4 && tokens.size() != 7)
                            emit_warning("[OBJReader] vertex has unknown format: " + line);

                    } else if (tokens[0][1] == 'n') {
                        if (tokens.size() >= 4) {
                            vNormals.push_back(ParseUtil::ToFloat(tokens[1]));
                            vNormals.push_back(ParseUtil::ToFloat(tokens[2]));
                            vNormals.push_back(ParseUtil::ToFloat(tokens[3]));
                        } 
                        if (tokens.size() != 4)
                            emit_warning("[OBJReader] normal has more than 3 coordinates: " + line);

                    } else if (tokens[0][1] == 't') {
                        if (tokens.size() >= 3) {
                            vUVs.push_back(ParseUtil::ToFloat(tokens[1]));
                            vUVs.push_back(ParseUtil::ToFloat(tokens[2]));
                            nMaxUVLength = std::max(nMaxUVLength, (int)tokens.size());
                        } 
                        if ( tokens.size() != 3 )
                            emit_warning("[OBJReader] UV has unknown format: " + line);
                    }


                } else if (tokens[0][0] == 'f') {
                    if ( tokens.size() < 4 ) {
                        emit_warning("[OBJReader] degenerate face specified : " + line);
                    } else if (tokens.size() == 4) {
                        Triangle tri = Triangle();
                        parse_triangle(tokens, partsBuf, tri);

                        tri.nGroupID = nActiveGroup;

                        if (activeMaterial != nullptr) {
                            tri.nMaterialID = activeMaterial->id;
                            UsedMaterials[activeMaterial->id] = activeMaterial->name;
                        }

                        vTriangles.push_back(tri);
                        if (tri.is_complex())
                            HasComplexVertices = true;

                    } else {
                        append_face(tokens, partsBuf, activeMaterial, nActiveGroup);
                    }

                } else if (tokens[0][0] == 'g') {
                    std::string sGroupName = (tokens.size() == 2) ? tokens[1] : line.substr( line.find(tokens[1]), std::string::npos );
                    if ( ContainsKey( GroupNames, sGroupName) ) {
                        nActiveGroup = GroupNames[sGroupName];
                    } else {
                        nActiveGroup = nGroupCounter;
                        GroupNames[sGroupName] = nGroupCounter++;
                    }

                } else if (tokens[0][0] == 'o') {
                    // TODO multi-object support

                } else if (tokens[0] == "mtllib" && options.ReadMaterials) {
                    if (MTLFileSearchPaths.size() == 0)
                        emit_warning("Materials requested but Material Search Paths not initialized!");
					std::string sMTLPathString = (tokens.size() == 2) ? tokens[1] :
						line.substr(line.find(tokens[1]), std::string::npos);
					std::string sFile = FindMTLFile(sMTLPathString);
                    if ( ! sFile.empty() ) {
                        IOReadResult result = ReadMaterials(sFile);
                        if (result.code != IOCode::Ok)
                            emit_warning("error parsing " + sFile + " : " + result.message);
                    } else
						emit_warning("material file " + sMTLPathString + " could not be found in material search paths");

                } else if (tokens[0] == "usemtl" && options.ReadMaterials) {
                    activeMaterial = find_material(tokens[1]);
                }

            } catch (const std::exception & e) {
                emit_warning("error parsing line " + std::to_string(nLines) + ": " + line + ", exception " + e.what());
            }

        }

        m_bOBJHasPerVertexColors = bVerticesHaveColors;
        m_bOBJHasTriangleGroups = (nActiveGroup != InvalidGroupID);
        m_nSetInvalidGroupsTo = nGroupCounter++;
        m_nUVComponents = nMaxUVLength;

        return IOReadResult(IOCode::Ok, "");
    }


    int parse_v(const std::string & sToken)
    {
        int vi = ParseUtil::ToInt(sToken);
        if (vi < 0)
            vi = ((int)vPositions.size() / 3) + vi + 1;
        return vi;
    }
    int parse_n(const std::string &  sToken)
    {
        int vi = ParseUtil::ToInt(sToken);
        if (vi < 0)
            vi = ((int)vNormals.size() / 3) + vi + 1;
        return vi;
    }
    int parse_u(const std::string &  sToken)
    {
        int vi = ParseUtil::ToInt(sToken);
        if (vi < 0)
            vi = ((int)vUVs.size() / 2) + vi + 1;
        return vi;
    }

    void append_face(const std::vector<std::string> & tokens,  
					 std::vector<std::string> & partsBuf,
					 std::shared_ptr<OBJMaterial> activeMaterial, int nActiveGroup)
    {
		int nMode = 0;
		if (ParseUtil::Contains(tokens[1], "//"))
			nMode = 1;
		else if (ParseUtil::Contains(tokens[1], '/'))
			nMode = 2;

        Triangle t = Triangle();
        t.clear();
        for ( int ti = 0; ti < tokens.size()-1; ++ti) {
            int j = (ti < 3) ? ti : 2;
            if (ti >= 3)
                t.move_vertex(2, 1);

            // parse next vertex
            if (nMode == 0) {
                // "f v1 v2 v3"
                t.set_vertex(j, parse_v(tokens[ti + 1]));

            } else if (nMode == 1) {
                // "f v1//vn1 v2//vn2 v3//vn3"
				ParseUtil::Split(tokens[ti + 1], '/', partsBuf, ParseUtil::Options::RemoveEmpty);
                t.set_vertex(j, parse_v(partsBuf[0]), parse_n(partsBuf[1]));

            } else if (nMode == 2) {
				ParseUtil::Split(tokens[ti + 1], '/', partsBuf, ParseUtil::Options::RemoveEmpty);
				if (partsBuf.size() == 2) {
                    // "f v1/vt1 v2/vt2 v3/vt3"
                    t.set_vertex(j, parse_v(partsBuf[0]), -1, parse_u(partsBuf[1]));
                } else if (partsBuf.size() == 3) {
                    // "f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3"
                    t.set_vertex(j, parse_v(partsBuf[0]), parse_n(partsBuf[2]), parse_u(partsBuf[1]));
                } else {
                    emit_warning("append_face unexpected face component " + tokens[j]);
                }
            }


            // do append
            if (ti >= 2) {
                if (activeMaterial != nullptr) {
                    t.nMaterialID = activeMaterial->id;
                    UsedMaterials[activeMaterial->id] = activeMaterial->name;
                }
                t.nGroupID = nActiveGroup;
                vTriangles.push_back(t);
                if (t.is_complex())
                    HasComplexVertices = true;
            }
        }
    }

    void parse_triangle(const std::vector<std::string> & tokens, 
						std::vector<std::string> & partsBuf, 
						Triangle & t ){
        int nMode = 0;
        if ( ParseUtil::Contains(tokens[1], "//") )
            nMode = 1;
        else if ( ParseUtil::Contains(tokens[1], '/') )
            nMode = 2;

        t.clear();

        for (int j = 0; j < 3; ++j) {
            if (nMode == 0) {
                // "f v1 v2 v3"
                t.set_vertex(j, parse_v(tokens[j + 1]));

            } else if (nMode == 1) {
                // "f v1//vn1 v2//vn2 v3//vn3"
				ParseUtil::Split(tokens[j+1], '/', partsBuf, ParseUtil::Options::RemoveEmpty);
                t.set_vertex(j, parse_v(partsBuf[0]), parse_n(partsBuf[1]));

            } else if (nMode == 2) {
				ParseUtil::Split(tokens[j + 1], '/', partsBuf, ParseUtil::Options::RemoveEmpty);
				if (partsBuf.size() == 2) {
                    // "f v1/vt1 v2/vt2 v3/vt3"
                    t.set_vertex(j, parse_v(partsBuf[0]), -1, parse_u(partsBuf[1]));
                } else if (partsBuf.size() == 3) {
                    // "f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3"
                    t.set_vertex(j, parse_v(partsBuf[0]), parse_n(partsBuf[2]), parse_u(partsBuf[1]));
                } else {
                    emit_warning("parse_triangle unexpected face component " + tokens[j]);
                }
            }
        }

    }


	std::string FindMTLFile(const std::string sMTLFilePath) {
		for ( const std::string & sPath : MTLFileSearchPaths ) {
			std::string sFullPath = FileUtil::PathCombine(sPath, sMTLFilePath);
			if ( FileUtil::FileExists(sFullPath) )
				return sFullPath;
		}
		return std::string();
	}

	IOReadResult ReadMaterials(const std::string & sPath)
	{
		throw std::exception("need to implement this...");
	}
	//IOReadResult ReadMaterials(const std::string & sPath)
	//{
 //       if (nWarningLevel >= 1)
 //           emit_warning("[OBJReader] ReadMaterials " + sPath);

 //       StreamReader reader;
 //       try {
 //           reader = StreamReader(sPath);
 //           if (reader.EndOfStream)
 //               return IOReadResult(IOCode.FileAccessError, "");
 //       } catch {
 //           return IOReadResult(IOCode.FileAccessError, "");
 //       }


 //       OBJMaterial curMaterial = null;

 //       while (reader.Peek() >= 0) {

 //           string line = reader.ReadLine();
 //           string[] tokens = line.Split((char[])null, StringSplitOptions.RemoveEmptyEntries);
 //           if (tokens.size() == 0)
 //               continue;

 //           if ( tokens[0][0] == '#' ) {
 //               continue;
 //           } else if (tokens[0] == "newmtl") {
 //               curMaterial = OBJMaterial();
 //               curMaterial.name = tokens[1];
 //               curMaterial.id = Materials.Count;

 //               if (Materials.ContainsKey(curMaterial.name))
 //                   emit_warning("Material file " + sPath + " / material " + curMaterial.name + " : already exists in Material set. Replacing.");
 //               if (nWarningLevel >= 1)
 //                   emit_warning("[OBJReader] parsing material " + curMaterial.name);

 //               Materials[curMaterial.name] = curMaterial;

 //           } else if (tokens[0] == "Ka") {
 //               if (curMaterial != null) curMaterial.Ka = parse_mtl_color(tokens);
 //           } else if (tokens[0] == "Kd") {
 //               if (curMaterial != null) curMaterial.Kd = parse_mtl_color(tokens);
 //           } else if (tokens[0] == "Ks") {
 //               if (curMaterial != null) curMaterial.Ks = parse_mtl_color(tokens);
 //           } else if (tokens[0] == "Ke") {
 //               if (curMaterial != null) curMaterial.Ke = parse_mtl_color(tokens);
 //           } else if (tokens[0] == "Tf") {
 //               if (curMaterial != null) curMaterial.Tf = parse_mtl_color(tokens);

 //           } else if (tokens[0] == "illum") {
 //               if (curMaterial != null) curMaterial.illum = int.Parse(tokens[1]);

 //           } else if (tokens[0] == "d") {
 //               if (curMaterial != null) curMaterial.d = Single.Parse(tokens[1]);
 //           } else if (tokens[0] == "Tr") {     // alternate to d/alpha, [Tr]ansparency is 1-d
 //               if (curMaterial != null) curMaterial.d = 1.0f - Single.Parse(tokens[1]);
 //           } else if (tokens[0] == "Ns") {
 //               if (curMaterial != null) curMaterial.Ns = Single.Parse(tokens[1]);
 //           } else if (tokens[0] == "sharpness") {
 //               if (curMaterial != null) curMaterial.sharpness = Single.Parse(tokens[1]);
 //           } else if (tokens[0] == "Ni") {
 //               if (curMaterial != null) curMaterial.Ni = Single.Parse(tokens[1]);

 //           } else if (tokens[0] == "map_Ka") {
	//			if (curMaterial != null) curMaterial.map_Ka = parse_mtl_path(line,tokens);
 //           } else if (tokens[0] == "map_Kd") {
	//			if (curMaterial != null) curMaterial.map_Kd = parse_mtl_path(line,tokens);
 //           } else if (tokens[0] == "map_Ks") {
	//			if (curMaterial != null) curMaterial.map_Ks = parse_mtl_path(line,tokens);
 //           } else if (tokens[0] == "map_Ke") {
	//			if (curMaterial != null) curMaterial.map_Ke = parse_mtl_path(line,tokens);
 //           } else if (tokens[0] == "map_d") {
	//			if (curMaterial != null) curMaterial.map_d = parse_mtl_path(line,tokens);
 //           } else if (tokens[0] == "map_Ns") {
	//			if (curMaterial != null) curMaterial.map_Ns = parse_mtl_path(line,tokens);

 //           } else if (tokens[0] == "bump" || tokens[0] == "map_bump") {
	//			if (curMaterial != null) curMaterial.bump = parse_mtl_path(line,tokens);
 //           } else if (tokens[0] == "disp") {
	//			if (curMaterial != null) curMaterial.disp = parse_mtl_path(line,tokens);
 //           } else if (tokens[0] == "decal") {
	//			if (curMaterial != null) curMaterial.decal = parse_mtl_path(line,tokens);
 //           } else if (tokens[0] == "refl") {
	//			if (curMaterial != null) curMaterial.refl = parse_mtl_path(line,tokens);
 //           } else {
 //               emit_warning("unknown material command " + tokens[0]);
 //           }

 //       }

 //       if (nWarningLevel >= 1)
 //           emit_warning("[OBJReader] ReadMaterials completed");

 //       return IOReadResult(IOCode.Ok, "ok");
	//}


	//string parse_mtl_path(const std::string & line, const std::vector<std::string> & tokens) {
	//	if ( tokens.size() == 2 )
	//		return tokens[1];
	//	else
	//		return line.Substring( line.IndexOf(tokens[1]) );
	//}

 //   Vector3f parse_mtl_color(const std::vector<std::string> & tokens)
 //   {
 //       if ( tokens[1] == "spectral" ) {
 //           emit_warning("OBJReader::parse_material_color : spectral color not supported!");
 //           return Vector3f(1, 0, 0);
 //       } else if (tokens[1] == "xyz" ) {
 //           emit_warning("OBJReader::parse_material_color : xyz color not supported!");
 //           return Vector3f(1, 0, 0);
 //       } else {
 //           float r = ParseUtil::ToFloat(tokens[1]);
 //           float g = ParseUtil::ToFloat(tokens[2]);
 //           float b = ParseUtil::ToFloat(tokens[3]);
 //           return Vector3f(r, g, b);
 //       }
 //   }



    std::shared_ptr<OBJMaterial> find_material(const std::string & sName)
    {
        if ( ContainsKey(Materials, sName))
            return Materials[sName];

		// [RMSC#] ??
        // try case-insensitive search
        //try {
        //    return Materials.First(x => String.Equals(x.Key, sName, StringComparison.OrdinalIgnoreCase)).Value;
        //} catch {
        //    // didn't work
        //}

        emit_warning("unknown material " + sName + " referenced");
        return nullptr;
    }




    void emit_warning(std::string sMessage)
    {
		auto sPrefix = sMessage.substr(0, 15);
        int nCount = ContainsKey(warningCount, sPrefix) ? warningCount[sPrefix] : 0;
        nCount++; 
		warningCount[sPrefix] = nCount;
        if (nCount > 10)
            return;
        else if (nCount == 10)
            sMessage += " (additional message surpressed)";

		if (warningEvent != nullptr)
			warningEvent(sMessage);
    }


};











}