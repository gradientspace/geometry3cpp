#pragma once

#include <g3types.h>

#include <string>
#include <map>

#include <dvector.h>
#include <refcount_vector.h>
#include <small_list_set.h>
#include <g3Debug.h>
#include <VectorUtil.h>
#include <index_util.h>
#include <iterator_util.h>



namespace g3
{

enum class MeshResult
{
	Ok = 0,
	Failed_NotAVertex = 1,
	Failed_NotATriangle = 2,
	Failed_NotAnEdge = 3,

	Failed_BrokenTopology = 10,
	Failed_HitValenceLimit = 11,

	Failed_IsBoundaryEdge = 20,
	Failed_FlippedEdgeExists = 21,
	Failed_IsBowtieVertex = 22,
	Failed_InvalidNeighbourhood = 23,       // these are all failures for CollapseEdge
	Failed_FoundDuplicateTriangle = 24,
	Failed_CollapseTetrahedron = 25,
	Failed_CollapseTriangle = 26,
	Failed_NotABoundaryEdge = 27,
	Failed_SameOrientation = 28,

	Failed_WouldCreateBowtie = 30,
	Failed_VertexAlreadyExists = 31,
	Failed_CannotAllocateVertex = 32,

	Failed_WouldCreateNonmanifoldEdge = 50,
	Failed_TriangleAlreadyExists = 51,
	Failed_CannotAllocateTriangle = 52
};


enum class MeshComponents
{
	None = 0,
	VertexNormals = 1,
	VertexColors = 2,
	VertexUVs = 4,
	FaceGroups = 8
};

enum class MeshHints
{
	None = 0,
	IsCompact = 1
};


/*
* Abstracts construction of meshes, so that we can construct different types, etc
*/
struct NewVertexInfo
{
	Vector3d v;
	Vector3f n, c;
	Vector2f uv;
	bool bHaveN, bHaveUV, bHaveC;

	NewVertexInfo() {
		this->v = Vector3d::Zero(); n = c = Vector3f::Zero(); uv = Vector2f::Zero();
		bHaveN = bHaveC = bHaveUV = false;
	}
	NewVertexInfo(Vector3d v) {
		this->v = v; n = c = Vector3f::Zero(); uv = Vector2f::Zero();
		bHaveN = bHaveC = bHaveUV = false;
	}
	NewVertexInfo(Vector3d v, Vector3f n) {
		this->v = v; this->n = n; c = Vector3f::Zero(); uv = Vector2f::Zero();
		bHaveN = true; bHaveC = bHaveUV = false;
	}
	NewVertexInfo(Vector3d v, Vector3f n, Vector3f c) {
		this->v = v; this->n = n; this->c = c; uv = Vector2f::Zero();
		bHaveN = bHaveC = true; bHaveUV = false;
	}
	NewVertexInfo(Vector3d v, Vector3f n, Vector3f c, Vector2f uv) {
		this->v = v; this->n = n; this->c = c; this->uv = uv;
		bHaveN = bHaveC = bHaveUV = true;
	}
};





//
// DMesh3 is a dynamic triangle mesh class. The mesh has has connectivity, 
//  is an indexed mesh, and allows for gaps in the index space.
//
// internally, all data is stored in POD-type buffers, except for the vertex->edge
// links, which are stored as List<int>'s. The arrays of POD data are stored in
// dvector's, so they grow in chunks, which is relatively efficient. The actual
// blocks are arrays, so they can be efficiently mem-copied into larger buffers
// if necessary.
//
// Reference counts for verts/tris/edges are stored as separate refcount_vector
// instances. 
//
// Vertices are stored as doubles, although this should be easily changed
// if necessary, as the internal data structure is not exposed
//
// Per-vertex Vertex Normals, Colors, and UVs are optional and stored as floats.
//
// For each vertex, vertex_edges[i] is the unordered list of connected edges. The
// elements of the list are indices into the edges list.
// This list is unsorted but can be traversed in-order (ie cw/ccw) at some additional cost. 
//
// Triangles are stored as 3 ints, with optionally a per-triangle integer group id.
//
// The edges of a triangle are similarly stored as 3 ints, in triangle_edes. If the 
// triangle is [v1,v2,v3], then the triangle edges [e1,e2,e3] are 
// e1=edge(v1,v2), e2=edge(v2,v3), e3=edge(v3,v1), where the e# are indexes into edges.
//
// Edges are stored as tuples of 4 ints. If the edge is between v1 and v2, with neighbour
// tris t1 and t2, then the edge is [min(v1,v2), max(v1,v2), t1, t2]. For a boundary
// edge, t2 is InvalidID. t1 is never InvalidID.
//
// Most of the class assumes that the mesh is manifold. Many functions will
// work if the topology is non-manifold, but behavior of operators like Split/Flip/Collapse
// edge is untested. 
//
// The function CheckValidity() does extensive sanity checking on the mesh data structure.
// Use this to test your code, both for mesh construction and editing!!
// 
//
// TODO:
//  - dvector w/ 'stride' option, so that we can guarantee that tuples are in single block.
//    The can have custom accessor that looks up entire tuple
class DMesh3
{
public:
	static constexpr double max_double = std::numeric_limits<double>::max();
	static constexpr int InvalidID = -1;
	static constexpr int NonManifoldID = -2;

	static const Vector3d InvalidVertex;
	static const Index3i InvalidTriangle;
	static const Index2i InvalidEdge;


protected:
	refcount_vector vertices_refcount;
    dvector<double> vertices;
	dvector<float> normals;
	dvector<float> colors;
	dvector<float> uv;

    // [TODO] this is optional if we only want to use this class as an iterable mesh-with-nbrs
    //   make it optional with a flag? (however find_edge depends on it...)
    small_list_set vertex_edges;

    refcount_vector triangles_refcount;
    dvector<int> triangles;
    dvector<int> triangle_edges;
	dvector<int> triangle_groups;

    refcount_vector edges_refcount;
    dvector<int> edges;

    int timestamp = 0;
    int shape_timestamp = 0;

    int max_group_id = 0;


    ///// <summary>
    ///// Support attaching arbitrary data to mesh. 
    ///// Note that metadata is currently **NOT** copied when copying a mesh.
    ///// </summary>
    //Dictionary<string, object> Metadata = nullptr;

public:
    DMesh3(bool bWantNormals = true, bool bWantColors = false, bool bWantUVs = false, bool bWantTriGroups = false)
    {
        vertices = dvector<double>();
		if ( bWantNormals)
			normals = dvector<float>();
		if ( bWantColors )
			colors = dvector<float>();
		if ( bWantUVs )
			uv = dvector<float>();

        vertex_edges = small_list_set();

        vertices_refcount = refcount_vector();

        triangles = dvector<int>();
        triangle_edges = dvector<int>();
        triangles_refcount = refcount_vector();
		if ( bWantTriGroups )
			triangle_groups = dvector<int>();
        max_group_id = 0;

        edges = dvector<int>();
        edges_refcount = refcount_vector();
    }
    DMesh3(MeshComponents flags) : 
        DMesh3( ((int)flags & (int)MeshComponents::VertexNormals) != 0,  ((int)flags & (int)MeshComponents::VertexColors) != 0,
                ((int)flags & (int)MeshComponents::VertexUVs) != 0,      ((int)flags & (int)MeshComponents::FaceGroups) != 0 )
    {
    }

    // normals/colors/uvs will only be copied if they exist
    DMesh3(const DMesh3 & copy, bool bCompact = false, bool bWantNormals = true, bool bWantColors = true, bool bWantUVs = true)
    {
        if (bCompact)
            CompactCopy(copy, bWantNormals, bWantColors, bWantUVs);
        else
            Copy(copy, bWantNormals, bWantColors, bWantUVs);
    }
    DMesh3(const DMesh3 & copy, bool bCompact, MeshComponents flags) : 
		DMesh3(copy, bCompact, ((int)flags & (int)MeshComponents::VertexNormals) != 0,  ((int)flags & (int)MeshComponents::VertexColors) != 0,
                ((int)flags & (int)MeshComponents::VertexUVs) != 0 )
    {
    }


    //DMesh3(IMesh copy, MeshHints hints, bool bWantNormals = true, bool bWantColors = true, bool bWantUVs = true)
    //{
    //    Copy(copy, hints, bWantNormals, bWantColors, bWantUVs);
    //}
    //DMesh3(IMesh copy, MeshHints hints, MeshComponents flags) : 
    //    this(copy, hints, (flags & MeshComponents.VertexNormals) != 0,  (flags & MeshComponents.VertexColors) != 0,
    //            (flags & MeshComponents.VertexUVs) != 0 )
    //{
    //}


public:
	//
	// iterators
	//   The functions vertices() / triangles() / edges() are provided so you can do:
	//      for ( int eid : edges() ) { ... }
	//   and other related begin() / end() idioms
	//
	typedef typename refcount_vector::index_enumerable vertex_iterator;
	vertex_iterator VertexIndices() const { return vertices_refcount.indices(); }

	typedef typename refcount_vector::index_enumerable triangle_iterator;
	triangle_iterator TriangleIndices() const { return triangles_refcount.indices(); }

	typedef typename refcount_vector::index_enumerable edge_iterator;
	edge_iterator EdgeIndices() const { return edges_refcount.indices(); }




protected:
	// internal

	void set_triangle(int tid, int v0, int v1, int v2)
	{
		int i = 3 * tid;
		triangles[i] = v0;
		triangles[i + 1] = v1;
		triangles[i + 2] = v2;
	}
	void set_triangle_edges(int tid, int e0, int e1, int e2)
	{
		int i = 3 * tid;
		triangle_edges[i] = e0;
		triangle_edges[i + 1] = e1;
		triangle_edges[i + 2] = e2;
	}

	int add_edge(int vA, int vB, int tA, int tB = InvalidID)
	{
		if (vB < vA) {
			int t = vB; vB = vA; vA = t;
		}
		int eid = edges_refcount.allocate();
		int i = 4 * eid;
		edges.insertAt(vA, i);
		edges.insertAt(vB, i + 1);
		edges.insertAt(tA, i + 2);
		edges.insertAt(tB, i + 3);

		vertex_edges.Insert(vA, eid);
		vertex_edges.Insert(vB, eid);
		return eid;
	}

	int replace_tri_vertex(int tID, int vOld, int vNew) {
		int i = 3 * tID;
		if (triangles[i] == vOld) { triangles[i] = vNew; return 0; }
		if (triangles[i + 1] == vOld) { triangles[i + 1] = vNew; return 1; }
		if (triangles[i + 2] == vOld) { triangles[i + 2] = vNew; return 2; }
		return -1;
	}

	int add_triangle_only(int a, int b, int c, int e0, int e1, int e2) {
		int tid = triangles_refcount.allocate();
		int i = 3 * tid;
		triangles.insertAt(c, i + 2);
		triangles.insertAt(b, i + 1);
		triangles.insertAt(a, i);
		triangle_edges.insertAt(e2, i + 2);
		triangle_edges.insertAt(e1, i + 1);
		triangle_edges.insertAt(e0, i + 0);
		return tid;
	}



	void allocate_edges_list(int vid)
	{
		if (vid < vertex_edges.Size())
			vertex_edges.Clear(vid);
		vertex_edges.AllocateAt(vid);
	}
	std::vector<int> vertex_edges_list(int vid) const
	{
		std::vector<int> list;
		for (int eid : vertex_edges.values(vid))
			list.push_back(vid);
		return list;
	}
	//List<int> vertex_vertices_list(int vid)
	//{
	//	List<int> vnbrs = List<int>();
	//	foreach(int eid in vertex_edges.ValueItr(vid))
	//		vnbrs.Add(edge_other_v(eid, vid));
	//	return vnbrs;
	//}


	void set_edge_vertices(int eID, int a, int b) {
		int i = 4 * eID;
		edges[i] = std::min(a, b);
		edges[i + 1] = std::max(a, b);
	}
	void set_edge_triangles(int eID, int t0, int t1) {
		int i = 4 * eID;
		edges[i + 2] = t0;
		edges[i + 3] = t1;
	}

	int replace_edge_vertex(int eID, int vOld, int vNew) {
		int i = 4 * eID;
		int a = edges[i], b = edges[i + 1];
		if (a == vOld) {
			edges[i] = std::min(b, vNew);
			edges[i + 1] = std::max(b, vNew);
			return 0;
		}
		else if (b == vOld) {
			edges[i] = std::min(a, vNew);
			edges[i + 1] = std::max(a, vNew);
			return 1;
		}
		else
			return -1;
	}


	int replace_edge_triangle(int eID, int tOld, int tNew) {
		int i = 4 * eID;
		int a = edges[i + 2], b = edges[i + 3];
		if (a == tOld) {
			if (tNew == InvalidID) {
				edges[i + 2] = b;
				edges[i + 3] = InvalidID;
			}
			else
				edges[i + 2] = tNew;
			return 0;
		}
		else if (b == tOld) {
			edges[i + 3] = tNew;
			return 1;
		}
		else
			return -1;
	}

	int replace_triangle_edge(int tID, int eOld, int eNew) {
		int i = 3 * tID;
		if (triangle_edges[i] == eOld) {
			triangle_edges[i] = eNew;
			return 0;
		}
		else if (triangle_edges[i + 1] == eOld) {
			triangle_edges[i + 1] = eNew;
			return 1;
		}
		else if (triangle_edges[i + 2] == eOld) {
			triangle_edges[i + 2] = eNew;
			return 2;
		}
		else
			return -1;
	}





public:
	// TODO make this work
    struct CompactInfo
    {
		std::map<int, int> MapV;
	};

    CompactInfo CompactCopy(const DMesh3 & copy, bool bNormals = true, bool bColors = true, bool bUVs = true)
    {
		// TODO can't do until CompactInfo works
        //if ( copy.IsCompact() ) {
        //    Copy(copy, bNormals, bColors, bUVs);
        //    CompactInfo ci = CompactInfo() { MapV = IdentityIndexMap() };
        //    return ci;
        //}

        vertices = dvector<double>();
        vertex_edges = small_list_set();
        vertices_refcount = refcount_vector();
        triangles = dvector<int>();
        triangle_edges = dvector<int>();
        triangles_refcount = refcount_vector();
        edges = dvector<int>();
        edges_refcount = refcount_vector();
        max_group_id = 0;

		normals = dvector<float>();
		colors = dvector<float>();
		uv = dvector<float>();
		triangle_groups = dvector<int>();

        // [TODO] if we ksome of these were dense we could copy directly...

        NewVertexInfo vinfo;
		std::vector<int> mapV; mapV.resize(copy.MaxVertexID());
		for ( int vid : VertexIndices() ) {
            copy.GetVertex(vid, vinfo, bNormals, bColors, bUVs);
            mapV[vid] = AppendVertex(vinfo);
        }

        // [TODO] would be much faster to explicitly copy triangle & edge data structures!!
		for ( int tid : TriangleIndices() ) {
            Index3i t = copy.GetTriangle(tid);
			t = Index3i(mapV[t.x()], mapV[t.y()], mapV[t.z()]);
            int g = (copy.HasTriangleGroups()) ? copy.GetTriangleGroup(tid) : InvalidID;
            AppendTriangle(t, g);
            max_group_id = std::max(max_group_id, g+1);
        }

        //return CompactInfo() {
        //    MapV = IndexMap(mapV, this.MaxVertexID)
        //};
		return CompactInfo();
    }


    void Copy(DMesh3 copy, bool bNormals = true, bool bColors = true, bool bUVs = true)
    {
        vertices = dvector<double>(copy.vertices);

        normals = (bNormals && copy.HasVertexNormals()) ? dvector<float>(copy.normals) : dvector<float>();
        colors = (bColors && copy.HasVertexColors()) ? dvector<float>(copy.colors) : dvector<float>();
        uv = (bUVs && copy.HasVertexUVs()) ? dvector<float>(copy.uv) : dvector<float>();

        vertices_refcount = refcount_vector(copy.vertices_refcount);

        vertex_edges = small_list_set(copy.vertex_edges);

        triangles = dvector<int>(copy.triangles);
        triangle_edges = dvector<int>(copy.triangle_edges);
        triangles_refcount = refcount_vector(copy.triangles_refcount);
        if (copy.HasTriangleGroups())
            triangle_groups = dvector<int>(copy.triangle_groups);
        max_group_id = copy.max_group_id;

        edges = dvector<int>(copy.edges);
        edges_refcount = refcount_vector(copy.edges_refcount);
    }


    /// <summary>
    /// Copy IMesh into this mesh. Currently always compacts.
    /// [TODO] if we get dense hint, we could be smarter w/ vertex map, etc
    /// </summary>
    //CompactInfo Copy(IMesh copy, MeshHints hints, bool bNormals = true, bool bColors = true, bool bUVs = true)
    //{
    //    vertices = dvector<double>();
    //    vertex_edges = small_list_set();
    //    vertices_refcount = refcount_vector();
    //    triangles = dvector<int>();
    //    triangle_edges = dvector<int>();
    //    triangles_refcount = refcount_vector();
    //    edges = dvector<int>();
    //    edges_refcount = refcount_vector();
    //    max_group_id = 0;

    //    normals = (bNormals && copy.HasVertexNormals) ? dvector<float>() : nullptr;
    //    colors = (bColors && copy.HasVertexColors) ? dvector<float>() : nullptr;
    //    uv = (bUVs && copy.HasVertexUVs) ? dvector<float>() : nullptr;
    //    triangle_groups = (copy.HasTriangleGroups) ? dvector<int>() : nullptr;


    //    // [TODO] if we ksome of these were dense we could copy directly...

    //    NewVertexInfo vinfo = NewVertexInfo();
    //    int[] mapV = int[copy.MaxVertexID];
    //    foreach (int vid in copy.VertexIndices()) {
    //        vinfo = copy.GetVertexAll(vid);
    //        mapV[vid] = AppendVertex(vinfo);
    //    }

    //    // [TODO] would be much faster to explicitly copy triangle & edge data structures!!

    //    foreach (int tid in copy.TriangleIndices()) {
    //        Index3i t = copy.GetTriangle(tid);
    //        t.a = mapV[t.a]; t.b = mapV[t.b]; t.c = mapV[t.c];
    //        int g = (copy.HasTriangleGroups) ? copy.GetTriangleGroup(tid) : InvalidID;
    //        AppendTriangle(t, g);
    //        max_group_id = std::max(max_group_id, g + 1);
    //    }

    //    return CompactInfo() {
    //        MapV = IndexMap(mapV, this.MaxVertexID)
    //    };
    //}



protected:
	void updateTimeStamp(bool bShapeChange) {
        timestamp++;
        if (bShapeChange)
            shape_timestamp++;
	}

public:
    /// <summary>
    /// Timestamp is incremented any time any change is made to the mesh
    /// </summary>
    int Timestamp() const {
        return timestamp;
    }

    /// <summary>
    /// ShapeTimestamp is incremented any time any vertex position is changed or the mesh topology is modified
    /// </summary>
    int ShapeTimestamp() const {
        return shape_timestamp;
    }


    // IMesh impl

    int VertexCount() const {
        return (int)vertices_refcount.count();
    }
    int TriangleCount() const {
        return (int)triangles_refcount.count();
    }
	int EdgeCount() const {
		return (int)edges_refcount.count();
	}

    // these values are (max_used+1), ie so an iteration should be < MaxTriangleID, not <=
	int MaxVertexID() const {
		return (int)vertices_refcount.max_index();
	}
	int MaxTriangleID() const {
		return (int)triangles_refcount.max_index();
	}
	int MaxEdgeID() const {
		return (int)edges_refcount.max_index();
	}
    int MaxGroupID() const {
        return max_group_id; 
    }

	bool HasVertexColors() const { return colors.size() == vertices.size(); }
	bool HasVertexNormals() const { return normals.size() == vertices.size(); }
	bool HasVertexUVs() const { return uv.size() / 2 == vertices.size() / 3; }
	bool HasTriangleGroups() const { return triangle_groups.size() == triangles.size() / 3; }

    int Components() const {
        int c = 0;
        if (HasVertexNormals()) c |= (int)MeshComponents::VertexNormals;
        if (HasVertexColors()) c |= (int)MeshComponents::VertexColors;
        if (HasVertexUVs()) c |= (int)MeshComponents::VertexUVs;
        if (HasTriangleGroups()) c |= (int)MeshComponents::FaceGroups;
        return c;
    }

    // info

    bool IsVertex(int vID) const {
        return vertices_refcount.isValid(vID);
    }
    bool IsTriangle(int tID) const {
        return triangles_refcount.isValid(tID);
    }
    bool IsEdge(int eID) const {
        return edges_refcount.isValid(eID);
    }

protected:
	void debug_check_is_vertex(int v) const {
		gDevAssert(IsVertex(v));
	}
	void debug_check_is_triangle(int t) const {
		gDevAssert(IsTriangle(t));
	}
	void debug_check_is_edge(int e) const {
		gDevAssert(IsEdge(e));
	}

public:

    // getters


    Vector3d GetVertex(int vID) const {
        debug_check_is_vertex(vID);
        int i = 3 * vID;
        return Vector3d(vertices[i], vertices[i + 1], vertices[i + 2]);
    }
    Vector3f GetVertexf(int vID) const {
        debug_check_is_vertex(vID);
        int i = 3 * vID;
        return Vector3f((float)vertices[i], (float)vertices[i + 1], (float)vertices[i + 2]);
    }

    void SetVertex(int vID, Vector3d vNewPos) {
		gDevAssert(IsFinite(vNewPos));
        debug_check_is_vertex(vID);

		int i = 3*vID;
		vertices[i] = vNewPos.x(); vertices[i+1] = vNewPos.y(); vertices[i+2] = vNewPos.z();
        updateTimeStamp(true);
	}

	Vector3f GetVertexNormal(int vID) const {
		if (HasVertexNormals() == false)
			return Vector3f::UnitY();

        debug_check_is_vertex(vID);
        int i = 3 * vID;
        return Vector3f(normals[i], normals[i + 1], normals[i + 2]);
	}

	void SetVertexNormal(int vID, Vector3f vNewNormal) {
		if ( HasVertexNormals() ) {
            debug_check_is_vertex(vID);
            int i = 3*vID;
			normals[i] = vNewNormal.x(); normals[i+1] = vNewNormal.y(); normals[i+2] = vNewNormal.z();
            updateTimeStamp(false);
		}
	}

    Vector3f GetVertexColor(int vID) const {
		if (HasVertexColors() == false)
			return Vector3f::Ones();
        debug_check_is_vertex(vID);
        int i = 3 * vID;
        return Vector3f(colors[i], colors[i + 1], colors[i + 2]);
	}

	void SetVertexColor(int vID, Vector3f vNewColor) {
		if ( HasVertexColors() ) {
            debug_check_is_vertex(vID);
            int i = 3*vID;
			colors[i] = vNewColor.x(); colors[i+1] = vNewColor.y(); colors[i+2] = vNewColor.z();
            updateTimeStamp(false);
		}
	}

	Vector2f GetVertexUV(int vID) const {
		if (HasVertexUVs() == false)
			return Vector2f::Zero();
        debug_check_is_vertex(vID);
        int i = 2 * vID;
        return Vector2f(uv[i], uv[i + 1]);
	}

	void SetVertexUV(int vID, Vector2f vNewUV) {
		if ( HasVertexUVs() ) {
            debug_check_is_vertex(vID);
            int i = 2*vID;
			uv[i] = vNewUV.x(); uv[i+1] = vNewUV.y();
            updateTimeStamp(false);
		}
	}

    bool GetVertex(int vID, NewVertexInfo & vinfo, bool bWantNormals, bool bWantColors, bool bWantUVs) const
    {
        if (vertices_refcount.isValid(vID) == false)
            return false;
        vinfo.v = Vector3d(vertices[3 * vID], vertices[3 * vID + 1], vertices[3 * vID + 2]);
        vinfo.bHaveN = vinfo.bHaveUV = vinfo.bHaveC = false;
        if (HasVertexNormals() && bWantNormals) {
            vinfo.bHaveN = true;
            vinfo.n = Vector3f(normals[3 * vID], normals[3 * vID + 1], normals[3 * vID + 2]);
        }
        if (HasVertexColors() && bWantColors) {
            vinfo.bHaveC = true;
            vinfo.c = Vector3f(colors[3 * vID], colors[3 * vID + 1], colors[3 * vID + 2]);
        }
        if (HasVertexUVs() && bWantUVs) {
            vinfo.bHaveUV = true;
            vinfo.uv = Vector2f(uv[2 * vID], uv[2 * vID + 1]);
        }
        return true;
    }

    int GetVtxEdgeCount(int vID) const {
        return vertices_refcount.isValid(vID) ? vertex_edges.Count(vID) : -1;
    }

    int GetMaxVtxEdgeCount() const {
        int max = 0;
		for ( int vid : VertexIndices() )
            max = std::max(max, vertex_edges.Count(vid));
        return max;
    }

	NewVertexInfo GetVertexAll(int i) const {
		NewVertexInfo vi = NewVertexInfo();
		vi.v = GetVertex(i);
		if ( HasVertexNormals() ) {
			vi.bHaveN = true;
			vi.n = GetVertexNormal(i);
		} else
			vi.bHaveN = false;
		if ( HasVertexColors()) {
			vi.bHaveC = true;
			vi.c = GetVertexColor(i);
		} else
			vi.bHaveC = false;
		if ( HasVertexUVs()) {
			vi.bHaveUV = true;
			vi.uv = GetVertexUV(i);
		} else
			vi.bHaveUV = false;
		return vi;
	}


    /// <summary>
    /// Compute a normal/tangent frame at vertex that is "stable" as long as
    /// the mesh topology doesn't change, meaning that one axis of the frame
    /// will be computed from projection of outgoing edge.
    /// Requires that vertex normals are available.
    /// by default, frame.Z is normal, and .X points along mesh edge
    /// if bFrameNormalY, then frame.Y is normal (X still points along mesh edge)
    /// </summary>
    Frame3d GetVertexFrame(int vID, bool bFrameNormalY = false) const
    {
        gDevAssert(HasVertexNormals());

        int vi = 3 * vID;
        Vector3d v = Vector3d(vertices[vi], vertices[vi + 1], vertices[vi + 2]);
        Vector3d normal = Vector3d(normals[vi], normals[vi + 1], normals[vi + 2]);
        int eid = vertex_edges.First(vID);
        int ovi = 3 * edge_other_v(eid, vID);
        Vector3d ov = Vector3d(vertices[ovi], vertices[ovi + 1], vertices[ovi + 2]);
        Vector3d edge = (ov - v);
        edge.normalize();

        Vector3d other = normal.cross(edge);
        edge = other.cross(normal);
        if (bFrameNormalY)
            return Frame3d(v, edge, normal, -other);
        else 
            return Frame3d(v, edge, other, normal);
    }




    Index3i GetTriangle(int tID) const {
        debug_check_is_triangle(tID);
        int i = 3 * tID;
        return Index3i(triangles[i], triangles[i + 1], triangles[i + 2]);
    }

    Index3i GetTriEdges(int tID) const {
        debug_check_is_triangle(tID);
        int i = 3 * tID;
        return Index3i(triangle_edges[i], triangle_edges[i + 1], triangle_edges[i + 2]);
    }

    int GetTriEdge(int tid, int j) const {
        debug_check_is_triangle(tid);
        return triangle_edges[3*tid+j];
    }


    Index3i GetTriNeighbourTris(int tID) const {
        if (triangles_refcount.isValid(tID)) {
            int tei = 3 * tID;
            Index3i nbr_t = Index3i::Zero();
            for (int j = 0; j < 3; ++j) {
                int ei = 4 * triangle_edges[tei + j];
                nbr_t[j] = (edges[ei + 2] == tID) ? edges[ei + 3] : edges[ei + 2];
            }
            return nbr_t;
        } else
            return InvalidTriangle;
    }

    //IEnumerable<int> TriTrianglesItr(int tID) {
    //    if (triangles_refcount.isValid(tID)) {
    //        int tei = 3 * tID;
    //        for (int j = 0; j < 3; ++j) {
    //            int ei = 4 * triangle_edges[tei + j];
    //            int nbr_t = (edges[ei + 2] == tID) ? edges[ei + 3] : edges[ei + 2];
    //            if (nbr_t != DMesh3.InvalidID)
    //                yield return nbr_t;
    //        }
    //    }
    //}



    int GetTriangleGroup(int tID) const { 
		return (HasTriangleGroups()) ? -1 
            : ( triangles_refcount.isValid(tID) ? triangle_groups[tID] : 0 );
	}

	void SetTriangleGroup(int tid, int group_id) {
		if (HasTriangleGroups()) {
            debug_check_is_triangle(tid);
            triangle_groups[tid] = group_id;
            max_group_id = std::max(max_group_id, group_id+1);
            updateTimeStamp(false);
		}
	}

    int AllocateTriangleGroup() {
        return ++max_group_id;
    }


    void GetTriVertices(int tID, Vector3d & v0, Vector3d & v1, Vector3d & v2) const {
        int ai = 3 * triangles[3 * tID];
        v0.x() = vertices[ai]; v0.y() = vertices[ai + 1]; v0.z() = vertices[ai + 2];
        int bi = 3 * triangles[3 * tID + 1];
        v1.x() = vertices[bi]; v1.y() = vertices[bi + 1]; v1.z() = vertices[bi + 2];
        int ci = 3 * triangles[3 * tID + 2];
        v2.x() = vertices[ci]; v2.y() = vertices[ci + 1]; v2.z() = vertices[ci + 2];
    }

    Vector3d GetTriVertex(int tid, int j) const {
        int a = triangles[3 * tid + j];
        return Vector3d(vertices[3 * a], vertices[3 * a + 1], vertices[3 * a + 2]);
    }

    Vector3d GetTriBaryPoint(int tID, double bary0, double bary1, double bary2) const {
        int ai = 3 * triangles[3 * tID], 
            bi = 3 * triangles[3 * tID + 1], 
            ci = 3 * triangles[3 * tID + 2];
        return Vector3d(
            (bary0*vertices[ai] + bary1*vertices[bi] + bary2*vertices[ci]),
            (bary0*vertices[ai + 1] + bary1*vertices[bi + 1] + bary2*vertices[ci + 1]),
            (bary0*vertices[ai + 2] + bary1*vertices[bi + 2] + bary2*vertices[ci + 2]));
    }

    Vector3d GetTriNormal(int tID) const
    {
		Vector3d v0, v1, v2;
        GetTriVertices(tID, v0, v1, v2);
        return g3::Normal(v0, v1, v2);
    }

    double GetTriArea(int tID) const
    {
		Vector3d v0, v1, v2;
        GetTriVertices(tID, v0, v1, v2);
        return Area(v0, v1, v2);
    }

	/// <summary>
	/// Compute triangle normal, area, and centroid all at once. Re-uses vertex
	/// lookups and computes normal & area simultaneously. *However* does not produce
	/// the same normal/area as separate calls, because of this.
	/// </summary>
	void GetTriInfo(int tID, Vector3d & normal, double & fArea, Vector3d & vCentroid) const
	{
		Vector3d v0, v1, v2;
		GetTriVertices(tID, v0, v1, v2);
		vCentroid = (1.0 / 3.0) * (v0 + v1 + v2);
		fArea = Area(v0, v1, v2);
		normal = Normal(v0, v1, v2);
		//normal = FastNormalArea(ref v0, ref v1, ref v2, out fArea);
	}


    /// <summary>
    /// interpolate vertex normals of triangle using barycentric coordinates
    /// </summary>
    Vector3d GetTriBaryNormal(int tID, double bary0, double bary1, double bary2) { 
        int ai = 3 * triangles[3 * tID], 
            bi = 3 * triangles[3 * tID + 1], 
            ci = 3 * triangles[3 * tID + 2];
        Vector3d n = Vector3d(
            (bary0*normals[ai] + bary1*normals[bi] + bary2*normals[ci]),
            (bary0*normals[ai + 1] + bary1*normals[bi + 1] + bary2*normals[ci + 1]),
            (bary0*normals[ai + 2] + bary1*normals[bi + 2] + bary2*normals[ci + 2]));
        n.normalize();
        return n;
    }

    /// <summary>
    /// efficiently compute centroid of triangle
    /// </summary>
    Vector3d GetTriCentroid(int tID) const
    {
        int ai = 3 * triangles[3 * tID], 
            bi = 3 * triangles[3 * tID + 1], 
            ci = 3 * triangles[3 * tID + 2];
        double f = (1.0 / 3.0);
        return Vector3d(
            (vertices[ai] + vertices[bi] + vertices[ci]) * f,
            (vertices[ai + 1] + vertices[bi + 1] + vertices[ci + 1]) * f,
            (vertices[ai + 2] + vertices[bi + 2] + vertices[ci + 2]) * f );
    }


    /// <summary>
    /// Compute interpolated vertex attributes at point of triangle
    /// </summary>
    void GetTriBaryPoint(int tID, double bary0, double bary1, double bary2, NewVertexInfo & vinfo)
    {
        vinfo = NewVertexInfo();
        int ai = 3 * triangles[3 * tID],
            bi = 3 * triangles[3 * tID + 1],
            ci = 3 * triangles[3 * tID + 2];
        vinfo.v = Vector3d(
            (bary0 * vertices[ai] + bary1 * vertices[bi] + bary2 * vertices[ci]),
            (bary0 * vertices[ai + 1] + bary1 * vertices[bi + 1] + bary2 * vertices[ci + 1]),
            (bary0 * vertices[ai + 2] + bary1 * vertices[bi + 2] + bary2 * vertices[ci + 2]));
        vinfo.bHaveN = HasVertexNormals();
        if (vinfo.bHaveN) {
            vinfo.n = Vector3f(
                (float)(bary0 * normals[ai] + bary1 * normals[bi] + bary2 * normals[ci]),
				(float)(bary0 * normals[ai + 1] + bary1 * normals[bi + 1] + bary2 * normals[ci + 1]),
				(float)(bary0 * normals[ai + 2] + bary1 * normals[bi + 2] + bary2 * normals[ci + 2]));
            vinfo.n.normalize();
        }
        vinfo.bHaveC = HasVertexColors();
        if (vinfo.bHaveC) {
            vinfo.c = Vector3f(
				(float)(bary0 * colors[ai] + bary1 * colors[bi] + bary2 * colors[ci]),
				(float)(bary0 * colors[ai + 1] + bary1 * colors[bi + 1] + bary2 * colors[ci + 1]),
				(float)(bary0 * colors[ai + 2] + bary1 * colors[bi + 2] + bary2 * colors[ci + 2]));
        }
        vinfo.bHaveUV = HasVertexUVs();
        if (vinfo.bHaveUV) {
            ai = 2 * triangles[3 * tID];
            bi = 2 * triangles[3 * tID + 1];
            ci = 2 * triangles[3 * tID + 2];
            vinfo.uv = Vector2f(
				(float)(bary0 * uv[ai] + bary1 * uv[bi] + bary2 * uv[ci]),
				(float)(bary0 * uv[ai + 1] + bary1 * uv[bi + 1] + bary2 * uv[ci + 1]));
        }
    }


    /// <summary>
    /// construct bounding box of triangle as efficiently as possible
    /// </summary>
    AxisAlignedBox3d GetTriBounds(int tID) const
    {
        int vi = 3 * triangles[3 * tID];
        double x = vertices[vi], y = vertices[vi + 1], z = vertices[vi + 2];
        double minx = x, maxx = x, miny = y, maxy = y, minz = z, maxz = z;
        for (int i = 1; i < 3; ++i) {
            vi = 3 * triangles[3 * tID + i];
            x = vertices[vi]; y = vertices[vi + 1]; z = vertices[vi + 2];
            if (x < minx) minx = x; else if (x > maxx) maxx = x;
            if (y < miny) miny = y; else if (y > maxy) maxy = y;
            if (z < minz) minz = z; else if (z > maxz) maxz = z;
        }
        return AxisAlignedBox3d(minx, miny, minz, maxx, maxy, maxz);
    }


    /// <summary>
    /// Construct stable frame at triangle centroid, where frame.Z is face normal,
    /// and frame.X is aligned with edge nEdge of triangle.
    /// </summary>
    Frame3d GetTriFrame(int tID, int nEdge = 0)
    {
        int ti = 3 * tID;
        int a = 3 * triangles[ti + (nEdge % 3)];
        int b = 3 * triangles[ti + ((nEdge+1) % 3)];
        int c = 3 * triangles[ti + ((nEdge+2) % 3)];
        Vector3d v1 = Vector3d(vertices[a], vertices[a + 1], vertices[a + 2]);
        Vector3d v2 = Vector3d(vertices[b], vertices[b + 1], vertices[b + 2]);
        Vector3d v3 = Vector3d(vertices[c], vertices[c + 1], vertices[c + 2]);

        Vector3d edge1 = v2 - v1;  edge1.normalize();
        Vector3d edge2 = v3 - v2;  edge2.normalize();
        Vector3d normal = edge1.cross(edge2); normal.normalize();

        Vector3d other = normal.cross(edge1);

        Vector3d center = (v1 + v2 + v3) / 3;
        return Frame3d(center, edge1, other, normal);
    }



    /// <summary>
    /// compute solid angle of oriented triangle tID relative to point p - see WindingNumber()
    /// </summary>
    double GetTriSolidAngle(int tID, const Vector3d & p) const
    {
        int ti = 3 * tID;
        int ta = 3 * triangles[ti];
        Vector3d a = Vector3d(vertices[ta] - p.x(), vertices[ta + 1] - p.y(), vertices[ta + 2] - p.z());
        int tb = 3 * triangles[ti + 1];
        Vector3d b = Vector3d(vertices[tb] - p.x(), vertices[tb + 1] - p.y(), vertices[tb + 2] - p.z());
        int tc = 3 * triangles[ti + 2];
        Vector3d c = Vector3d(vertices[tc] - p.x(), vertices[tc + 1] - p.y(), vertices[tc + 2] - p.z());
        // note: top and bottom are reversed here from formula in the paper? but it doesn't work otherwise...
        double la = a.norm(), lb = b.norm(), lc = c.norm();
        double bottom = (la * lb * lc) + a.dot(b) * lc + b.dot(c) * la + c.dot(a) * lb;
        double top = a.x() * (b.y() * c.z() - c.y() * b.z()) - a.y() * (b.x() * c.z() - c.x() * b.z()) + a.z() * (b.x() * c.y() - c.x() * b.y());
        return 2.0 * atan2(top, bottom);
    }



    /// <summary>
    /// compute internal angle at vertex i of triangle (where i is 0,1,2);
    /// TODO can be more efficient here, probably...
    /// </summary>
    double GetTriInternalAngleR(int tID, int i)
    {
        int ti = 3 * tID;
        int ta = 3 * triangles[ti];
        Vector3d a = Vector3d(vertices[ta], vertices[ta + 1], vertices[ta + 2]);
        int tb = 3 * triangles[ti + 1];
        Vector3d b = Vector3d(vertices[tb], vertices[tb + 1], vertices[tb + 2]);
        int tc = 3 * triangles[ti + 2];
        Vector3d c = Vector3d(vertices[tc], vertices[tc + 1], vertices[tc + 2]);
        if ( i == 0 )
            return VectorAngleR( (b-a).normalized(), (c-a).normalized());
        else if ( i == 1 )
            return VectorAngleR( (a-b).normalized(), (c-b).normalized());
        else
            return VectorAngleR( (a-c).normalized(), (b-c).normalized());
    }



    Index2i GetEdgeV(int eID) const {
        debug_check_is_edge(eID);
        int i = 4 * eID;
        return Index2i(edges[i], edges[i + 1]);
    }
    bool GetEdgeV(int eID, Vector3d & a, Vector3d & b) const {
        debug_check_is_edge(eID);
        int iv0 = 3 * edges[4 * eID];
        a.x() = vertices[iv0]; a.y() = vertices[iv0 + 1]; a.z() = vertices[iv0 + 2];
        int iv1 = 3 * edges[4 * eID + 1];
        b.x() = vertices[iv1]; b.y() = vertices[iv1 + 1]; b.z() = vertices[iv1 + 2];
        return true;
    }

    Index2i GetEdgeT(int eID) const {
        debug_check_is_edge(eID);
        int i = 4 * eID;
        return Index2i(edges[i + 2], edges[i + 3]);
    }

    /// <summary>
    /// return [v0,v1,t0,t1], or Index4i.Max if eid is invalid
    /// </summary>
    Index4i GetEdge(int eID) const
    {
        debug_check_is_edge(eID);
        int i = 4 * eID;
        return Index4i(edges[i], edges[i + 1], edges[i + 2], edges[i + 3]);
    }

	bool GetEdge(int eID, int & a, int & b, int & t0, int & t1) const {
        debug_check_is_edge(eID);
		int i = eID*4;
		a = edges[i]; b = edges[i+1]; t0 = edges[i+2]; t1 = edges[i+3];
		return true;
	}

    // return same indices as GetEdgeV, but oriented based on attached triangle
    Index2i GetOrientedBoundaryEdgeV(int eID) const
    {
        if ( edges_refcount.isValid(eID) ) {
            int ei = 4 * eID;
            if ( edges[ei+3] == InvalidID) {
                int a = edges[ei], b = edges[ei + 1];
                int ti = 3 * edges[ei + 2];
                Index3i tri = Index3i(triangles[ti], triangles[ti + 1], triangles[ti + 2]);
                int ai = find_edge_index_in_tri(a, b, tri);
                return Index2i(tri[ai], tri[(ai + 1) % 3]);
            }
        }
        gDevAssert(false);
        return InvalidEdge;
    }
			
    // average of 1 or 2 face normals
    Vector3d GetEdgeNormal(int eID) const
    {
        if (edges_refcount.isValid(eID)) {
            int ei = 4 * eID;
            Vector3d n = GetTriNormal(edges[ei + 2]);
            if (edges[ei + 3] != InvalidID) {
                n += GetTriNormal(edges[ei + 3]);
                n.normalize();
            }
            return n;
        }
        gDevAssert(false);
        return Vector3d::Zero();
    }

	Vector3d GetEdgePoint(int eID, double t) const
	{
		if (edges_refcount.isValid(eID)) {
			int ei = 4 * eID;
			int iv0 = 3 * edges[ei];
			int iv1 = 3 * edges[ei + 1];
			double mt = 1.0 - t;
			return Vector3d(
				mt*vertices[iv0] + t*vertices[iv1],
				mt*vertices[iv0 + 1] + t*vertices[iv1 + 1],
				mt*vertices[iv0 + 2] + t*vertices[iv1 + 2]);
		}
        gDevAssert(false);
		return Vector3d::Zero();
	}


    // mesh-building


    /// <summary>
    /// Append vertex at position, returns vid
    /// </summary>
    int AppendVertex(const Vector3d & v) {
		NewVertexInfo vinfo(v);
        return AppendVertex(vinfo);
    }

    /// <summary>
    /// Append vertex at position and other fields, returns vid
    /// </summary>
    int AppendVertex(const NewVertexInfo & info)
    {
        int vid = vertices_refcount.allocate();
		int i = 3*vid;
        vertices.insertAt(info.v[2], i + 2);
        vertices.insertAt(info.v[1], i + 1);
        vertices.insertAt(info.v[0], i);

		if ( HasVertexNormals() ) {
			Vector3f n = (info.bHaveN) ? info.n : Vector3f::UnitY();
			normals.insertAt(n[2], i + 2);
			normals.insertAt(n[1], i + 1);
			normals.insertAt(n[0], i);
		}

		if ( HasVertexColors() ) {
			Vector3f c = (info.bHaveC) ? info.c : Vector3f::Ones();
			colors.insertAt(c[2], i + 2);
			colors.insertAt(c[1], i + 1);
			colors.insertAt(c[0], i);
		}

		if ( HasVertexUVs() ) {
			Vector2f u = (info.bHaveUV) ? info.uv : Vector2f::Zero();
			int j = 2*vid;
			uv.insertAt(u[1], j + 1);
			uv.insertAt(u[0], j);
		}

        allocate_edges_list(vid);

        updateTimeStamp(true);
        return vid;
    }


    /// <summary>
    /// copy vertex fromVID from existing source mesh, returns vid
    /// </summary>
    int AppendVertex(DMesh3 from, int fromVID)
    {
        int bi = 3 * fromVID;

        int vid = vertices_refcount.allocate();
		int i = 3*vid;
        vertices.insertAt(from.vertices[bi+2], i + 2);
        vertices.insertAt(from.vertices[bi+1], i + 1);
        vertices.insertAt(from.vertices[bi], i);
		if ( HasVertexNormals() ) {
            if ( from.HasVertexNormals() ) {
                normals.insertAt(from.normals[bi + 2], i + 2);
                normals.insertAt(from.normals[bi + 1], i + 1);
                normals.insertAt(from.normals[bi], i);
            } else {
                normals.insertAt(0, i + 2);
                normals.insertAt(1, i + 1);       // y-up
                normals.insertAt(0, i);
            }
		}

		if ( HasVertexColors() ) {
            if (from.HasVertexColors()) {
                colors.insertAt(from.colors[bi + 2], i + 2);
                colors.insertAt(from.colors[bi + 1], i + 1);
                colors.insertAt(from.colors[bi], i);
            } else {
                colors.insertAt(1, i + 2);
                colors.insertAt(1, i + 1);       // white
                colors.insertAt(1, i);
            }
		}

		if ( HasVertexUVs() ) {
			int j = 2*vid;
            if (from.HasVertexUVs()) {
                int bj = 2 * fromVID;
                uv.insertAt(from.uv[bj + 1], j + 1);
                uv.insertAt(from.uv[bj], j);
            } else {
                uv.insertAt(0, j + 1);
                uv.insertAt(0, j);
            }
		}

        allocate_edges_list(vid);

        updateTimeStamp(true);
        return vid;
    }



    /// <summary>
    /// insert vertex at given index, assuming it is unused
    /// If bUnsafe, we use fast id allocation that does not update free list.
    /// You should only be using this between BeginUnsafeVerticesInsert() / EndUnsafeVerticesInsert() calls
    /// </summary>
    MeshResult InsertVertex(int vid, const NewVertexInfo & info, bool bUnsafe = false)
    {
        if (vertices_refcount.isValid(vid))
            return MeshResult::Failed_VertexAlreadyExists;

        bool bOK = (bUnsafe) ? vertices_refcount.allocate_at_unsafe(vid) :
                                vertices_refcount.allocate_at(vid);
        if (bOK == false)
            return MeshResult::Failed_CannotAllocateVertex;

        int i = 3 * vid;
        vertices.insertAt(info.v[2], i + 2);
        vertices.insertAt(info.v[1], i + 1);
        vertices.insertAt(info.v[0], i);

        if (HasVertexNormals()) {
            Vector3f n = (info.bHaveN) ? info.n : Vector3f::UnitY();
            normals.insertAt(n[2], i + 2);
            normals.insertAt(n[1], i + 1);
            normals.insertAt(n[0], i);
        }

        if (HasVertexColors()) {
            Vector3f c = (info.bHaveC) ? info.c : Vector3f::Ones();
            colors.insertAt(c[2], i + 2);
            colors.insertAt(c[1], i + 1);
            colors.insertAt(c[0], i);
        }

        if (HasVertexUVs()) {
            Vector2f u = (info.bHaveUV) ? info.uv : Vector2f::Zero();
            int j = 2 * vid;
            uv.insertAt(u[1], j + 1);
            uv.insertAt(u[0], j);
        }

        allocate_edges_list(vid);

        updateTimeStamp(true);
        return MeshResult::Ok;
    }



    virtual void BeginUnsafeVerticesInsert() {
        // do nothing...
    }
    virtual void EndUnsafeVerticesInsert() {
        vertices_refcount.rebuild_free_list();
    }



    int AppendTriangle(int v0, int v1, int v2, int gid = -1) {
        return AppendTriangle(Index3i(v0, v1, v2), gid);
    }
    int AppendTriangle(Index3i tv, int gid = -1) {
        if (IsVertex(tv[0]) == false || IsVertex(tv[1]) == false || IsVertex(tv[2]) == false) {
            gDevAssert(false);
            return InvalidID;
        }
        if (tv[0] == tv[1] || tv[0] == tv[2] || tv[1] == tv[2]) {
            gDevAssert(false);
            return InvalidID;
        }

        // look up edges. if any already have two triangles, this would 
        // create non-manifold geometry and so we do not allow it
        int e0 = find_edge(tv[0], tv[1]);
        int e1 = find_edge(tv[1], tv[2]);
        int e2 = find_edge(tv[2], tv[0]);
        if ((e0 != InvalidID && IsBoundaryEdge(e0) == false)
                || (e1 != InvalidID && IsBoundaryEdge(e1) == false)
                || (e2 != InvalidID && IsBoundaryEdge(e2) == false)) {
            return NonManifoldID;
        }

        // now safe to insert triangle
        int tid = triangles_refcount.allocate();
		int i = 3*tid;
        triangles.insertAt(tv[2], i + 2);
        triangles.insertAt(tv[1], i + 1);
        triangles.insertAt(tv[0], i);
        if ( HasTriangleGroups() ) {
            triangle_groups.insertAt(gid, tid);
            max_group_id = std::max(max_group_id, gid+1);
        }

        // increment ref counts and update/create edges
        vertices_refcount.increment(tv[0]);
        vertices_refcount.increment(tv[1]);
        vertices_refcount.increment(tv[2]);

        add_tri_edge(tid, tv[0], tv[1], 0, e0);
        add_tri_edge(tid, tv[1], tv[2], 1, e1);
        add_tri_edge(tid, tv[2], tv[0], 2, e2);

        updateTimeStamp(true);
        return tid;
    }
    // helper fn for above, just makes code cleaner
    void add_tri_edge(int tid, int v0, int v1, int j, int eid)
    {
        if (eid != InvalidID) {
            edges[4 * eid + 3] = tid;
            triangle_edges.insertAt(eid, 3 * tid + j);
        } else
            triangle_edges.insertAt(add_edge(v0, v1, tid), 3 * tid + j);
    }





    /// <summary>
    /// Insert triangle at given index, assuming it is unused.
    /// If bUnsafe, we use fast id allocation that does not update free list.
    /// You should only be using this between BeginUnsafeTrianglesInsert() / EndUnsafeTrianglesInsert() calls
    /// </summary>
    MeshResult InsertTriangle(int tid, Index3i tv, int gid = -1, bool bUnsafe = false)
    {
        if (triangles_refcount.isValid(tid))
            return MeshResult::Failed_TriangleAlreadyExists;

        if (IsVertex(tv[0]) == false || IsVertex(tv[1]) == false || IsVertex(tv[2]) == false) {
            gDevAssert(false);
            return MeshResult::Failed_NotAVertex;
        }
        if (tv[0] == tv[1] || tv[0] == tv[2] || tv[1] == tv[2]) {
            gDevAssert(false);
            return MeshResult::Failed_InvalidNeighbourhood;
        }

        // look up edges. if any already have two triangles, this would 
        // create non-manifold geometry and so we do not allow it
        int e0 = find_edge(tv[0], tv[1]);
        int e1 = find_edge(tv[1], tv[2]);
        int e2 = find_edge(tv[2], tv[0]);
        if ((e0 != InvalidID && IsBoundaryEdge(e0) == false)
                || (e1 != InvalidID && IsBoundaryEdge(e1) == false)
                || (e2 != InvalidID && IsBoundaryEdge(e2) == false)) {
            return MeshResult::Failed_WouldCreateNonmanifoldEdge;
        }

        bool bOK = (bUnsafe) ? triangles_refcount.allocate_at_unsafe(tid) :
                                triangles_refcount.allocate_at(tid);
        if (bOK == false)
            return MeshResult::Failed_CannotAllocateTriangle;

        // now safe to insert triangle
        int i = 3 * tid;
        triangles.insertAt(tv[2], i + 2);
        triangles.insertAt(tv[1], i + 1);
        triangles.insertAt(tv[0], i);
        if (HasTriangleGroups()) {
            triangle_groups.insertAt(gid, tid);
            max_group_id = std::max(max_group_id, gid + 1);
        }

        // increment ref counts and update/create edges
        vertices_refcount.increment(tv[0]);
        vertices_refcount.increment(tv[1]);
        vertices_refcount.increment(tv[2]);

        add_tri_edge(tid, tv[0], tv[1], 0, e0);
        add_tri_edge(tid, tv[1], tv[2], 1, e1);
        add_tri_edge(tid, tv[2], tv[0], 2, e2);

        updateTimeStamp(true);
        return MeshResult::Ok;
    }


    virtual void BeginUnsafeTrianglesInsert() {
        // do nothing...
    }
    virtual void EndUnsafeTrianglesInsert() {
        triangles_refcount.rebuild_free_list();
    }





    void EnableVertexNormals(Vector3f initial_normal)
    {
        if (HasVertexNormals())
            return;
        normals = dvector<float>();
        int NV = MaxVertexID();
        normals.resize(3*NV);
        for (int i = 0; i < NV; ++i) {
            int vi = 3 * i;
            normals[vi] = initial_normal.x();
            normals[vi + 1] = initial_normal.y();
            normals[vi + 2] = initial_normal.z();
        }
    }
    void DiscardVertexNormals() {
        normals = dvector<float>();
    }

    void EnableVertexColors(Vector3f initial_color)
    {
        if (HasVertexColors())
            return;
        colors = dvector<float>();
        int NV = MaxVertexID();
        colors.resize(3*NV);
        for (int i = 0; i < NV; ++i) {
            int vi = 3 * i;
            colors[vi] = initial_color.x();
            colors[vi + 1] = initial_color.y();
            colors[vi + 2] = initial_color.z();
        }
    }
    void DiscardVertexColors() {
        colors = dvector<float>();
    }

    void EnableVertexUVs(Vector2f initial_uv)
    {
        if (HasVertexUVs())
            return;
        uv = dvector<float>();
        int NV = MaxVertexID();
        uv.resize(2*NV);
        for (int i = 0; i < NV; ++i) {
            int vi = 2 * i;
            uv[vi] = initial_uv.x();
            uv[vi + 1] = initial_uv.y();
        }
    }
    void DiscardVertexUVs() {
        uv = dvector<float>();
    }

    void EnableTriangleGroups(int initial_group = 0)
    {
        if (HasTriangleGroups())
            return;
        triangle_groups = dvector<int>();
        int NT = MaxTriangleID();
        triangle_groups.resize(NT);
        for (int i = 0; i < NT; ++i)
            triangle_groups[i] = initial_group;
        max_group_id = 0;
    }
    void DiscardTriangleGroups() {
        triangle_groups = dvector<int>();
        max_group_id = 0;
    }






    // iterators

	/// <summary> Enumerate boundary edges </summary>
	refcount_vector::filtered_enumerable BoundaryEdgeIndices() {
		return edges_refcount.filtered_indices( [=](int eid) {
			return edges[4 * eid + 3] == InvalidID;
		});
	}

	template <typename T> using value_iteration = refcount_vector::mapped_enumerable<T>;

    /// <summary> Enumerate vertices </summary>
	value_iteration<Vector3d> Vertices() {
		return vertices_refcount.mapped_indices<Vector3d>( [=](int vid) {
			int i = 3 * vid;
			return Vector3d(vertices[i], vertices[i + 1], vertices[i + 2]);
		});
    }

    /// <summary> Enumerate triangles </summary>
	value_iteration<Index3i> Triangles() {
		return triangles_refcount.mapped_indices<Index3i>( [=](int tid) {
			int i = 3 * tid;
			return Index3i(triangles[i], triangles[i + 1], triangles[i + 2]);
		});
    }

	/// <summary> Enumerate edges. return value is [v0,v1,t0,t1], where t1 will be InvalidID if this is a boundary edge </summary>
	value_iteration<Index4i> Edges() {
		return edges_refcount.mapped_indices<Index4i>([=](int eid) {
			int i = 4 * eid;
			return Index4i(edges[i], edges[i + 1], edges[i + 2], edges[i + 3]);
		});
	}


    // queries

    /// <summary>
    /// Find edgeid for edge [a,b]
    /// </summary>
    int FindEdge(int vA, int vB) const {
        debug_check_is_vertex(vA);
        debug_check_is_vertex(vB);
        return find_edge(vA, vB);
    }

    /// <summary>
    /// Find edgeid for edge [a,b] from triangle that contains the edge.
    /// This is faster than FindEdge() because it is constant-time
    /// </summary>
    int FindEdgeFromTri(int vA, int vB, int tID) const  {
        return find_edge_from_tri(vA, vB, tID);
    }

	/// <summary>
    /// If edge has vertices [a,b], and is connected two triangles [a,b,c] and [a,b,d],
    /// this returns [c,d], or [c,InvalidID] for a boundary edge
    /// </summary>
    Index2i GetEdgeOpposingV(int eID) const
    {
        // [TODO] there was a comment here saying this does more work than necessary??
		// ** it is important that verts returned maintain [c,d] order!!
		int i = 4*eID;
        int a = edges[i], b = edges[i + 1];
        int t0 = edges[i + 2], t1 = edges[i + 3];
		int c = find_tri_other_vtx(a, b, triangles, t0);
        if (t1 != InvalidID) {
			int d = find_tri_other_vtx(a, b, triangles, t1);
            return Index2i(c, d);
        } else
            return Index2i(c, InvalidID);
    }


    /// <summary>
    /// Find triangle made up of any permutation of vertices [a,b,c]
    /// </summary>
    int FindTriangle(int a, int b, int c) const
    {
        int eid = find_edge(a, b);
        if (eid == InvalidID)
            return InvalidID;
        int ei = 4 * eid;

        // triangles attached to edge [a,b] must contain verts a and b...
        int ti = 3 * edges[ei + 2];
        if (triangles[ti] == c || triangles[ti + 1] == c || triangles[ti + 2] == c )
            return edges[ei + 2];
        if (edges[ei + 3] != InvalidID) {
            ti = 3 * edges[ei + 3];
            if (triangles[ti] == c || triangles[ti + 1] == c || triangles[ti + 2] == c )
                return edges[ei + 3];
        }

        return InvalidID;
    }



    /// <summary>
    /// Enumerate "other" vertices of edges connected to vertex (ie vertex one-ring)
    /// </summary>
	small_list_set::value_enumerable VtxVerticesItr(int vID) const {
		gDevAssert(vertices_refcount.isValid(vID));
		return vertex_edges.values(vID, [vID, this](int eid) { return edge_other_v(eid, vID); } );
	}


    /// <summary>
    /// Enumerate edge ids connected to vertex (ie edge one-ring)
    /// </summary>
	small_list_set::value_enumerable VtxEdgesItr(int vID) const {
		gDevAssert(vertices_refcount.isValid(vID));
		return vertex_edges.values(vID);
    }


    /// <summary>
    /// Returns count of boundary edges at vertex, and 
    /// the first two boundary edges if found. 
    /// If return is > 2, call VtxAllBoundaryEdges
    /// </summary>
    int VtxBoundaryEdges(int vID, int & e0, int & e1) const
    {
        if ( vertices_refcount.isValid(vID) ) {
            int count = 0;
			for (int eid : vertex_edges.values(vID) ) {
                int ei = 4 * eid;
                if ( edges[ei+3] == InvalidID ) {
                    if (count == 0)
                        e0 = eid;
                    else if (count == 1)
                        e1 = eid;
                    count++;
                }
            }
            return count;
        }
        gDevAssert(false);
        return -1;
    }

    /// <summary>
    /// Find edge ids of boundary edges connected to vertex.
    /// e needs to be large enough (ie call VtxBoundaryEdges, or as large as max one-ring)
    /// returns count, ie number of elements of e that were filled
    /// </summary>
    int VtxAllBoundaryEdges(int vID, int * e) const
    {
        if (vertices_refcount.isValid(vID)) {
            int count = 0;
			for (int eid : vertex_edges.values(vID)) {
                int ei = 4 * eid;
                if ( edges[ei+3] == InvalidID ) 
                    e[count++] = eid;
            }
            return count;
        }
        gDevAssert(false);
        return -1;
    }


    /// <summary>
    /// Get triangle one-ring at vertex. 
    /// bUseOrientation is more efficient but returns incorrect result if vertex is a bowtie
    /// </summary>
    MeshResult GetVtxTriangles(int vID, std::vector<int> & vTriangles, bool bUseOrientation) const
    {
        if (!IsVertex(vID))
            return MeshResult::Failed_NotAVertex;

        if (bUseOrientation) {
			for (int eid : vertex_edges.values(vID)) {
                int vOther = edge_other_v(eid, vID);
				int i = 4*eid;
                int et0 = edges[i + 2];
                if (tri_has_sequential_v(et0, vID, vOther))
                    vTriangles.push_back(et0);
                int et1 = edges[i + 3];
                if (et1 != InvalidID && tri_has_sequential_v(et1, vID, vOther))
                    vTriangles.push_back(et1);
            }
        } else {
            // brute-force method
			for (int eid : vertex_edges.values(vID)) {
				int i = 4*eid;					
                int t0 = edges[i + 2];
                if ( std::find(vTriangles.begin(), vTriangles.end(), t0) == vTriangles.end() )
                    vTriangles.push_back(t0);
                int t1 = edges[i + 3];
                if (t1 != InvalidID && std::find(vTriangles.begin(), vTriangles.end(), t1) == vTriangles.end())
                    vTriangles.push_back(t1);
            }
        }
        return MeshResult::Ok;
    }


    /// <summary>
    /// return # of triangles attached to vID, or -1 if invalid vertex
    /// if bBruteForce = true, explicitly checks, which creates a list and is expensive
    /// default is false, uses orientation, no memory allocation
    /// </summary>
    int GetVtxTriangleCount(int vID, bool bBruteForce = false) const
    {
        if ( bBruteForce ) {
            std::vector<int> vTriangles;
            if (GetVtxTriangles(vID, vTriangles, false) != MeshResult::Ok)
                return -1;
            return (int)vTriangles.size();
        }

        if (!IsVertex(vID))
            return -1;
        int N = 0;
		for (int eid : vertex_edges.values(vID)) {
            int vOther = edge_other_v(eid, vID);
			int i = 4*eid;
            int et0 = edges[i + 2];
            if (tri_has_sequential_v(et0, vID, vOther))
                N++;
            int et1 = edges[i + 3];
            if (et1 != InvalidID && tri_has_sequential_v(et1, vID, vOther))
                N++;
        }
        return N;
    }

	
	using vtx_triangles_enumerable = expand_enumerable<int,int, small_list_set::value_iterator>;

    /// <summary>
    /// iterate over triangle IDs of vertex one-ring
    /// </summary>
	vtx_triangles_enumerable VtxTrianglesItr(int vID) {
		gDevAssert(vertices_refcount.isValid(vID));

		std::function<int(int, int&)> expand_f = [=](int eid, int & k) {
			int vOther = edge_other_v(eid, vID);
			int i = 4 * eid;
			if (k == -1) {
				int et0 = edges[i + 2];
				if (tri_has_sequential_v(et0, vID, vOther)) {
					k = 0;
					return et0;
				}
			}
			if (k != 1) {
				int et1 = edges[i + 3];
				if (et1 != InvalidID && tri_has_sequential_v(et1, vID, vOther)) {
					k = 1;
					return et1;
				}
			}
			k = -1;
			return -1;
		};

		return vtx_triangles_enumerable(vertex_edges.values(vID), expand_f);
	}


    /// <summary>
    ///  from edge and vert, returns other vert, two opposing verts, and two triangles
    /// </summary>
    void GetVtxNbrhood(int eID, int vID, int & vOther, int & oppV1, int & oppV2, int & t1, int & t2) const
	{
		int i = 4*eID;
		vOther = (edges[i] == vID) ? edges[i+1] : edges[i];
		t1 = edges[i + 2];
		oppV1 = find_tri_other_vtx(vID, vOther, triangles, t1);
		t2 = edges[i + 3];
		if ( t2 != InvalidID )
			oppV2 = find_tri_other_vtx(vID, vOther, triangles, t2);
		else
			t2 = InvalidID;
	}


    /// <summary>
    /// Fastest possible one-ring centroid. This is used inside many other algorithms
    /// so it helps to have it be maximally efficient
    /// </summary>
    void VtxOneRingCentroid(int vID, Vector3d & centroid) const
    {
        centroid = Vector3d::Zero();
        if (vertices_refcount.isValid(vID)) {
            int n = 0;
			for (int eid : vertex_edges.values(vID)) {
                int other_idx = 3 * edge_other_v(eid, vID);
                centroid.x() += vertices[other_idx];
                centroid.y() += vertices[other_idx + 1];
                centroid.z() += vertices[other_idx + 2];
                n++;
            }
            if (n > 0) {
				centroid *= 1.0 / n;
            }
        }
    }



    bool tri_has_v(int tID, int vID) const {
		int i = 3*tID;
        return triangles[i] == vID 
            || triangles[i + 1] == vID
            || triangles[i + 2] == vID;
    }

    bool tri_is_boundary(int tID) const {
		int i = 3*tID;
        return IsBoundaryEdge(triangle_edges[i])
            || IsBoundaryEdge(triangle_edges[i + 1])
            || IsBoundaryEdge(triangle_edges[i + 2]);
    }

    bool tri_has_neighbour_t(int tCheck, int tNbr) const {
		int i = 3*tCheck;
        return edge_has_t(triangle_edges[i], tNbr)
            || edge_has_t(triangle_edges[i + 1], tNbr)
            || edge_has_t(triangle_edges[i + 2], tNbr);
    }

    bool tri_has_sequential_v(int tID, int vA, int vB) const
    {
		int i = 3*tID;
        int v0 = triangles[i], v1 = triangles[i + 1], v2 = triangles[i + 2];
        if (v0 == vA && v1 == vB) return true;
        if (v1 == vA && v2 == vB) return true;
        if (v2 == vA && v0 == vB) return true;
        return false;
    }

	//! returns edge ID
	int find_tri_neighbour_edge(int tID, int vA, int vB) const
	{
		int i = 3*tID;
		int tv0 = triangles[i], tv1 = triangles[i+1];
		if ( same_pair_unordered(tv0, tv1, vA, vB) ) return triangle_edges[3*tID];
		int tv2 = triangles[i+2];
		if ( same_pair_unordered(tv1, tv2, vA, vB) ) return triangle_edges[3*tID+1];
		if ( same_pair_unordered(tv2, tv0, vA, vB) ) return triangle_edges[3*tID+2];
		return InvalidID;	
	}

	// returns 0/1/2
	int find_tri_neighbour_index(int tID, int vA, int vB) const
	{
		int i = 3*tID;
		int tv0 = triangles[i], tv1 = triangles[i+1];
		if ( same_pair_unordered(tv0, tv1, vA, vB) ) return 0;
		int tv2 = triangles[i+2];
		if ( same_pair_unordered(tv1, tv2, vA, vB) ) return 1;
		if ( same_pair_unordered(tv2, tv0, vA, vB) ) return 2;
		return InvalidID;	
	}


    bool IsBoundaryEdge(int eid) const {
        return edges[4 * eid + 3] == InvalidID;
    }

    bool edge_has_v(int eid, int vid) const {
		int i = 4*eid;
        return (edges[i] == vid) || (edges[i + 1] == vid);
    }
    bool edge_has_t(int eid, int tid) const {
		int i = 4*eid;
        return (edges[i + 2] == tid) || (edges[i + 3] == tid);
    }
    int edge_other_v(int eID, int vID) const
    {
		int i = 4*eID;
        int ev0 = edges[i], ev1 = edges[i + 1];
        return (ev0 == vID) ? ev1 : ((ev1 == vID) ? ev0 : InvalidID);
    }
    int edge_other_t(int eID, int tid) const {
		int i = 4*eID;
        int et0 = edges[i + 2], et1 = edges[i + 3];
        return (et0 == tid) ? et1 : ((et1 == tid) ? et0 : InvalidID);
    }


    bool IsBoundaryVertex(int vID) const {
		for (int eid : vertex_edges.values(vID)) {
            if (edges[4 * eid + 3] == InvalidID)
                return true;
        }
        return false;
    }


    /// <summary>
    /// Returns true if any edge of triangle is a boundary edge
    /// </summary>
    bool IsBoundaryTriangle(int tID) const
    {
        debug_check_is_triangle(tID);
        int i = 3 * tID;
        return IsBoundaryEdge(triangle_edges[i]) || IsBoundaryEdge(triangle_edges[i + 1]) || IsBoundaryEdge(triangle_edges[i + 2]);
    }



    int find_edge(int vA, int vB) const
    {
        // [RMS] edge vertices must be sorted (min,max),
        //   that means we only need one index-check in inner loop.
        //   commented out code is robust to incorrect ordering, but slower.
        int vO = std::max(vA, vB);
        int vI = std::min(vA, vB);
		for (int eid : vertex_edges.values(vI)) {
            if (edges[4 * eid + 1] == vO)
                //if (edge_has_v(eid, vO))
                return eid;
        }
        return InvalidID;

        // this is slower, likely because it creates func<> every time. can we do w/o that?
        //return vertex_edges.Find(vI, (eid) => { return edges[4 * eid + 1] == vO; }, InvalidID);
    }

    int find_edge_from_tri(int vA, int vB, int tID) const
    {
        int i = 3 * tID;
        int t0 = triangles[i], t1 = triangles[i + 1];
        if (same_pair_unordered(vA, vB, t0, t1))
            return triangle_edges[i];
        int t2 = triangles[i + 2];
        if (same_pair_unordered(vA, vB, t1, t2))
            return triangle_edges[i+1];
        if (same_pair_unordered(vA, vB, t2, t0))
            return triangle_edges[i+2];
        return InvalidID;
    }





    // queries


    /// <summary>
    /// Returns true if the two triangles connected to edge have different group IDs
    /// </summary>
    bool IsGroupBoundaryEdge(int eID) const
    {
		gDevAssert(IsEdge(eID));
		gDevAssert(HasTriangleGroups());

        int et1 = edges[4 * eID + 3];
        if (et1 == InvalidID)
            return false;
        int g1 = triangle_groups[et1];
        int et0 = edges[4 * eID + 2];
        int g0 = triangle_groups[et0];
        return g1 != g0;
    }


    /// <summary>
    /// returns true if vertex has more than one tri group in its tri nbrhood
    /// </summary>
    bool IsGroupBoundaryVertex(int vID) const
    {
		gDevAssert(IsVertex(vID));
		gDevAssert(HasTriangleGroups());

		int group_id = InvalidGroupID;
		for (int eID : vertex_edges.values(vID)) {
            int et0 = edges[4 * eID + 2];
            int g0 = triangle_groups[et0];
            if (group_id != g0) {
                if (group_id == InvalidGroupID)
                    group_id = g0;
                else
                    return true;        // saw multiple group IDs
            }
            int et1 = edges[4 * eID + 3];
            if (et1 != InvalidID) {
                int g1 = triangle_groups[et1];
                if (group_id != g1)
                    return true;        // saw multiple group IDs
            }
        }
        return false;
    }



    /// <summary>
    /// returns true if more than two group boundary edges meet at vertex (ie 3+ groups meet at this vertex)
    /// </summary>
    bool IsGroupJunctionVertex(int vID) const
    {
		gDevAssert(IsVertex(vID));
		gDevAssert(HasTriangleGroups());

		Index2i groups(InvalidGroupID, InvalidGroupID);
		for (int eID : vertex_edges.values(vID)) {
            Index2i et = Index2i(edges[4 * eID + 2], edges[4 * eID + 3]);
            for (int k = 0; k < 2; ++k) {
                if (et[k] == InvalidID)
                    continue;
                int g0 = triangle_groups[et[k]];
                if (g0 != groups[0] && g0 != groups[1]) {
                    if (groups[0] != InvalidGroupID && groups[1] != InvalidGroupID)
                        return true;
                    if (groups[0] == InvalidGroupID)
                        groups[0] = g0;
                    else
                        groups[1] = g0;
                }
            }
        }
        return false;
    }


    /// <summary>
    /// returns up to 4 group IDs at vertex. Returns false if > 4 encountered
    /// </summary>
    bool GetVertexGroups(int vID, Index4i & groups) const
    {
		gDevAssert(IsVertex(vID));
		gDevAssert(HasTriangleGroups());

		groups = Index4i(InvalidGroupID, InvalidGroupID, InvalidGroupID, InvalidGroupID);
        int ng = 0;

		for (int eID : vertex_edges.values(vID)) {
            int et0 = edges[4 * eID + 2];
            int g0 = triangle_groups[et0];
            if ( Contains(groups, g0) == false )
                groups[ng++] = g0;
            if (ng == 4)
                return false;
            int et1 = edges[4 * eID + 3];
            if ( et1 != InvalidID ) {
                int g1 = triangle_groups[et1];
                if ( Contains(groups, g1) == false)
                    groups[ng++] = g1;
                if (ng == 4)
                    return false;
            }
        }
        return true;
    }



    /// <summary>
    /// returns all group IDs at vertex
    /// </summary>
    bool GetAllVertexGroups(int vID, std::vector<int> & groups) const
    {
		gDevAssert(IsVertex(vID));
		gDevAssert(HasTriangleGroups());

		for (int eID : vertex_edges.values(vID)) {
            int et0 = edges[4 * eID + 2];
            int g0 = triangle_groups[et0];
            if ( Contains(groups, g0) == false)
                groups.push_back(g0);
            int et1 = edges[4 * eID + 3];
            if ( et1 != InvalidID ) {
                int g1 = triangle_groups[et1];
                if (Contains(groups, g1) == false)
                    groups.push_back(g1);
            }
        }
        return true;
    }
	std::vector<int> GetAllVertexGroups(int vID) const {
		std::vector<int> result;
        GetAllVertexGroups(vID, result);
        return result;
    }



    /// <summary>
    /// returns true if vID is a "bowtie" vertex, ie multiple disjoint triangle sets in one-ring
    /// </summary>
    bool IsBowtieVertex(int vID) const
    {
		gDevAssert(vertices_refcount.isValid(vID));

		int nEdges = vertex_edges.Count(vID);
		if (nEdges == 0)
			return false;

        // find a boundary edge to start at
        int start_eid = -1;
        bool start_at_boundary = false;
		for (int eid : vertex_edges.values(vID)) {
            if (edges[4 * eid + 3] == InvalidID) {
                start_at_boundary = true;
                start_eid = eid;
                break;
            }
        }
        // if no boundary edge, start at arbitrary edge
        if (start_eid == -1)
            start_eid = vertex_edges.First(vID);
        // initial triangle
        int start_tid = edges[4 * start_eid + 2];

        int prev_tid = start_tid;
        int prev_eid = start_eid;

        // walk forward to next edge. if we hit start edge or boundary edge,
        // we are done the walk. count number of edges as we go.
        int count = 1;
        while (true) {
            int i = 3 * prev_tid;
            Index3i tv = Index3i(triangles[i], triangles[i+1], triangles[i+2]);
            Index3i te = Index3i(triangle_edges[i], triangle_edges[i+1], triangle_edges[i+2]);
            int vert_idx = find_tri_index(vID, tv);
            int e1 = te[vert_idx], e2 = te[(vert_idx+2) % 3];
            int next_eid = (e1 == prev_eid) ? e2 : e1;
            if (next_eid == start_eid)
                break;
            Index2i next_eid_tris = GetEdgeT(next_eid);
            int next_tid = (next_eid_tris[0] == prev_tid) ? next_eid_tris[1] : next_eid_tris[0];
            if (next_tid == InvalidID) {
                break;
            }
            prev_eid = next_eid;
            prev_tid = next_tid;
            count++;
        }

        // if we did not see all edges at vertex, we have a bowtie
        int target_count = (start_at_boundary) ? nEdges - 1 : nEdges;
        bool is_bowtie = (target_count != count);
        return is_bowtie;
    }


    /// <summary>
    /// Computes bounding box of all vertices.
    /// </summary>
    AxisAlignedBox3d GetBounds() const
    {
        double x = 0, y = 0, z = 0;
        for ( int vi : VertexIndices() ) {
            x = vertices[3*vi]; y = vertices[3*vi + 1]; z = vertices[3*vi + 2];
            break;
        }
        double minx = x, maxx = x, miny = y, maxy = y, minz = z, maxz = z;
		for (int vi : VertexIndices()) {
            x = vertices[3*vi]; y = vertices[3*vi + 1]; z = vertices[3*vi + 2];
            if (x < minx) minx = x; else if (x > maxx) maxx = x;
            if (y < miny) miny = y; else if (y > maxy) maxy = y;
            if (z < minz) minz = z; else if (z > maxz) maxz = z;
        }
        return AxisAlignedBox3d(minx, miny, minz, maxx, maxy, maxz);
    }

    AxisAlignedBox3d cached_bounds;
    int cached_bounds_timestamp = -1;

    /// <summary>
    /// cached bounding box, lazily re-computed on access if mesh has changed
    /// </summary>
    AxisAlignedBox3d CachedBounds()
    {
        if (cached_bounds_timestamp != Timestamp()) {
            cached_bounds = GetBounds();
            cached_bounds_timestamp = Timestamp();
        }
        return cached_bounds;
    }


    bool cached_is_closed = false;
    int cached_is_closed_timestamp = -1;

    bool IsClosed() const {
        if (TriangleCount() == 0)
            return false;
        // [RMS] under possibly-mistaken belief that foreach() has some overhead...
        if (MaxEdgeID() / EdgeCount() > 5) {
            for (int eid : EdgeIndices() )
                if (IsBoundaryEdge(eid))
                    return false;
        } else {
            int N = MaxEdgeID();
            for (int i = 0; i < N; ++i)
                if (edges_refcount.isValid(i) && IsBoundaryEdge(i))
                    return false;
        }
        return true;            
    }

    bool CachedIsClosed() {
        if (cached_is_closed_timestamp != Timestamp()) {
            cached_is_closed = IsClosed();
            cached_is_closed_timestamp = Timestamp();
        }
        return cached_is_closed;
    }




    /// <summary> returns true if vertices, edges, and triangles are all "dense" (Count == MaxID) </summary>
    bool IsCompact() const {
        return vertices_refcount.is_dense() && edges_refcount.is_dense() && triangles_refcount.is_dense();
    }

    /// <summary> Returns true if vertex count == max vertex id </summary>
    bool IsCompactV() const {
        return vertices_refcount.is_dense();
    }

    /// <summary> returns true if triangle count == max triangle id </summary>
    bool IsCompactT() const {
        return triangles_refcount.is_dense();
    }

    /// <summary> returns measure of compactness in range [0,1], where 1 is fully compacted </summary>
    double CompactMetric() const {
        return ((double)VertexCount() / (double)MaxVertexID() + (double)TriangleCount() / (double)MaxTriangleID()) * 0.5;
    }



    /// <summary>
    /// Compute mesh winding number, from Jacobson et al, Robust Inside-Outside Segmentation using Generalized Winding Numbers
    /// http://igl.ethz.ch/projects/winding-number/
    /// returns ~0 for points outside a closed, consistently oriented mesh, and a positive or negative integer
    /// for points inside, with value > 1 depending on how many "times" the point inside the mesh (like in 2D polygon winding)
    /// </summary>
    double WindingNumber(Vector3d v) const
    {
        double sum = 0;
        for ( int tid : TriangleIndices() )
            sum += GetTriSolidAngle(tid, v);
        return sum / (4.0 * Math<double>::PI);
    }




    // Metadata support
	// [RMS] disabled for now

    //bool HasMetadata {
    //    get { return Metadata != nullptr && Metadata.Keys.Count > 0; }
    //}
    //void AttachMetadata(string key, object o)
    //{
    //    if (Metadata == nullptr)
    //        Metadata = Dictionary<string, object>();
    //    Metadata.Add(key, o);
    //}
    //object FindMetadata(string key)
    //{
    //    if (Metadata == nullptr)
    //        return nullptr;
    //    object o = nullptr;
    //    bool bFound = Metadata.TryGetValue(key, out o);
    //    return (bFound) ? o : nullptr;
    //}
    //bool RemoveMetadata(string key)
    //{
    //    if (Metadata == nullptr)
    //        return false;
    //    return Metadata.Remove(key);
    //}
    //void ClearMetadata()
    //{
    //    if (Metadata != nullptr) {
    //        Metadata.Clear();
    //        Metadata = nullptr;
    //    }
    //}







    // direct access to internal dvectors - dangerous!!

    const dvector<double> & VerticesBuffer() {
        return vertices;
    }
    const refcount_vector & VerticesRefCounts() {
        return vertices_refcount;
    }
    const dvector<float> & NormalsBuffer() {
        return normals;
    }
    const dvector<float> & ColorsBuffer() {
        return colors;
    }
    const dvector<float> & UVBuffer() {
        return uv;
    }

    const dvector<int> & TrianglesBuffer() {
        return triangles;
    }
    const refcount_vector & TrianglesRefCounts() {
        return triangles_refcount;
    }
    const dvector<int> & GroupsBuffer() {
        return triangle_groups;
    }

    const dvector<int> & EdgesBuffer() {
        return edges;
    }
    const refcount_vector & EdgesRefCounts() {
        return edges_refcount;
    }
    const small_list_set & VertexEdges() {
        return vertex_edges;
    }



    /// <summary>
    /// Rebuild mesh topology.
    /// assumes that we have initialized vertices, triangles, and edges buffers,
    /// and edges refcounts. Rebuilds vertex and tri refcounts, triangle edges, vertex edges.
    /// </summary>
    void RebuildFromEdgeRefcounts()
    {
        int MaxVID = (int)vertices.length() / 3;
        int MaxTID = (int)triangles.length() / 3;

        triangle_edges.resize(triangles.length());
        triangles_refcount.RawRefCounts().resize(MaxTID);

        vertex_edges.Resize(MaxVID);
        vertices_refcount.RawRefCounts().resize(MaxVID);

        int MaxEID = (int)edges.length() / 4;
        for ( int eid = 0; eid < MaxEID; ++eid ) {
            if (edges_refcount.isValid(eid) == false)
                continue;
            int va = edges[4 * eid];
            int vb = edges[4 * eid + 1];
            int t0 = edges[4 * eid + 2];
            int t1 = edges[4 * eid + 3];

            // set vertex and tri refcounts to 1
            // find edges [a,b] in each triangle and set its tri-edge to this edge

            if (vertices_refcount.isValidUnsafe(va) == false) {
                allocate_edges_list(va);
                vertices_refcount.set_Unsafe(va, 1);
            }
            if (vertices_refcount.isValidUnsafe(vb) == false) {
                allocate_edges_list(vb);
                vertices_refcount.set_Unsafe(vb, 1);
            }
            triangles_refcount.set_Unsafe(t0, 1);
            Index3i tri0 = GetTriangle(t0);
            int idx0 = find_edge_index_in_tri(va, vb, tri0);
            triangle_edges[3 * t0 + idx0] = eid;

            if (t1 != InvalidID) {
                triangles_refcount.set_Unsafe(t1, 1);
                Index3i tri1 = GetTriangle(t1);
                int idx1 = find_edge_index_in_tri(va, vb, tri1);
                triangle_edges[3 * t1 + idx1] = eid;
            }

            // add this edge to both vertices
            vertex_edges.Insert(va, eid);
            vertex_edges.Insert(vb, eid);
        }

        // iterate over triangles and increment vtx refcount for each tri
        bool has_groups = HasTriangleGroups();
        max_group_id = 0;
        for ( int tid = 0; tid < MaxTID; ++tid ) {
            if (triangles_refcount.isValid(tid) == false)
                continue;
            int a = triangles[3 * tid], b = triangles[3 * tid + 1], c = triangles[3 * tid + 2];
            vertices_refcount.increment(a);
            vertices_refcount.increment(b);
            vertices_refcount.increment(c);

            if (has_groups)
                max_group_id = std::max(max_group_id, triangle_groups[tid]);
        }
        max_group_id++;

        vertices_refcount.rebuild_free_list();
        triangles_refcount.rebuild_free_list();
        edges_refcount.rebuild_free_list();

        updateTimeStamp(true);
    }



    /// <summary>
    /// Compact mesh in-place, by moving vertices around and rewriting indices.
    /// Should be faster if the amount of compacting is not too significant, and
    /// is useful in some places.
    /// [TODO] vertex_edges is not compacted. does not affect indices, but does keep memory.
    /// 
    /// If bComputeCompactInfo=false, the returned CompactInfo is not initialized
    /// </summary>
	// [RMS] disabled for now
    CompactInfo CompactInPlace(bool bComputeCompactInfo = false)
    {
        //IndexMap mapV = (bComputeCompactInfo) ? IndexMap(MaxVertexID, VertexCount) : nullptr;
		std::map<int, int> mapV;
        CompactInfo ci = CompactInfo();
        ci.MapV = mapV;

        // find first free vertex, and last used vertex
        int iLastV = MaxVertexID() - 1, iCurV = 0;
        while (vertices_refcount.isValidUnsafe(iLastV) == false)
            iLastV--;
        while (vertices_refcount.isValidUnsafe(iCurV))
            iCurV++;

        dvector<short> vref = vertices_refcount.RawRefCounts();

        while (iCurV < iLastV) {
            int kc = iCurV * 3, kl = iLastV * 3;
            vertices[kc] = vertices[kl];  vertices[kc+1] = vertices[kl+1];  vertices[kc+2] = vertices[kl+2];
            if ( HasVertexNormals() ) {
                normals[kc] = normals[kl];  normals[kc+1] = normals[kl+1];  normals[kc+2] = normals[kl+2];
            }
            if (HasVertexColors()) {
                colors[kc] = colors[kl];  colors[kc+1] = colors[kl+1];  colors[kc+2] = colors[kl+2];
            }
            if (HasVertexUVs()) {
                int ukc = iCurV * 2, ukl = iLastV * 2;
                uv[ukc] = uv[ukl]; uv[ukc+1] = uv[ukl+1];
            }

            for ( int eid : vertex_edges.values(iLastV) ) {
                // replace vertex in edges
                replace_edge_vertex(eid, iLastV, iCurV);

                // replace vertex in triangles
                int t0 = edges[4*eid + 2];
                replace_tri_vertex(t0, iLastV, iCurV);
                int t1 = edges[4*eid + 3];
                if ( t1 != InvalidID )
                    replace_tri_vertex(t1, iLastV, iCurV);
            }

            // shift vertex refcount to position
            vref[iCurV] = vref[iLastV];
            vref[iLastV] = refcount_vector::invalid;

            // move edge list
            vertex_edges.Move(iLastV, iCurV);

            if (bComputeCompactInfo)
                mapV[iLastV] = iCurV;

            // move cur forward one, last back one, and  then search for next valid
            iLastV--; iCurV++;
            while (vertices_refcount.isValidUnsafe(iLastV) == false)
                iLastV--;
            while (vertices_refcount.isValidUnsafe(iCurV) && iCurV < iLastV)
                iCurV++;
        }

        // trim vertices data structures
        vertices_refcount.trim(VertexCount());
        vertices.resize(VertexCount() * 3);
        if (HasVertexNormals())
            normals.resize(VertexCount() * 3);
        if (HasVertexColors())
            colors.resize(VertexCount() * 3);
        if (HasVertexUVs())
            uv.resize(VertexCount() * 2);

        // [TODO] vertex_edges!!!

        /** shift triangles **/

        // find first free triangle, and last valid triangle
        int iLastT = MaxTriangleID() - 1, iCurT = 0;
        while (triangles_refcount.isValidUnsafe(iLastT) == false)
            iLastT--;
        while (triangles_refcount.isValidUnsafe(iCurT))
            iCurT++;

        dvector<short> tref = triangles_refcount.RawRefCounts();

        while (iCurT < iLastT) {
            int kc = iCurT * 3, kl = iLastT * 3;

            // shift triangle
            for (int j = 0; j < 3; ++j) {
                triangles[kc + j] = triangles[kl + j];
                triangle_edges[kc + j] = triangle_edges[kl + j];
            }
            if (HasTriangleGroups())
                triangle_groups[iCurT] = triangle_groups[iLastT];

            // update edges
            for ( int j = 0; j < 3; ++j ) {
                int eid = triangle_edges[kc + j];
                replace_edge_triangle(eid, iLastT, iCurT);
            }

            // shift triangle refcount to position
            tref[iCurT] = tref[iLastT];
            tref[iLastT] = refcount_vector::invalid;

            // move cur forward one, last back one, and  then search for next valid
            iLastT--; iCurT++;
            while (triangles_refcount.isValidUnsafe(iLastT) == false)
                iLastT--;
            while (triangles_refcount.isValidUnsafe(iCurT) && iCurT < iLastT)
                iCurT++;
        }

        // trim triangles data structures
        triangles_refcount.trim(TriangleCount());
        triangles.resize(TriangleCount() * 3);
        triangle_edges.resize(TriangleCount() * 3);
        if (HasTriangleGroups())
            triangle_groups.resize(TriangleCount());

        /** shift edges **/

        // find first free edge, and last used edge
        int iLastE = MaxEdgeID() - 1, iCurE = 0;
        while (edges_refcount.isValidUnsafe(iLastE) == false)
            iLastE--;
        while (edges_refcount.isValidUnsafe(iCurE))
            iCurE++;

        dvector<short> eref = edges_refcount.RawRefCounts();

        while (iCurE < iLastE) {
            int kc = iCurE * 4, kl = iLastE * 4;

            // shift edge
            for (int j = 0; j < 4; ++j) {
                edges[kc + j] = edges[kl + j];
            }

            // replace edge in vertex edges lists
            int v0 = edges[kc], v1 = edges[kc + 1];
            vertex_edges.Replace(v0, [iLastE](int eid) { return eid == iLastE; }, iCurE);
            vertex_edges.Replace(v1, [iLastE](int eid) { return eid == iLastE; }, iCurE);

            // replace edge in triangles
            replace_triangle_edge(edges[kc + 2], iLastE, iCurE);
            if (edges[kc + 3] != InvalidID)
                replace_triangle_edge(edges[kc + 3], iLastE, iCurE);

            // shift triangle refcount to position
            eref[iCurE] = eref[iLastE];
            eref[iLastE] = refcount_vector::invalid;

            // move cur forward one, last back one, and  then search for next valid
            iLastE--; iCurE++;
            while (edges_refcount.isValidUnsafe(iLastE) == false)
                iLastE--;
            while (edges_refcount.isValidUnsafe(iCurE) && iCurE < iLastE)
                iCurE++;
        }

        // trim edge data structures
        edges_refcount.trim(EdgeCount());
        edges.resize(EdgeCount() * 4);

        return ci;
    }











    // edits

    MeshResult ReverseTriOrientation(int tID) {
        if (!IsTriangle(tID))
            return MeshResult::Failed_NotATriangle;
        internal_reverse_tri_orientation(tID);
        updateTimeStamp(true);
        return MeshResult::Ok;
    }
    void internal_reverse_tri_orientation(int tID) {
        Index3i t = GetTriangle(tID);
        set_triangle(tID, t[1], t[0], t[2]);
        Index3i te = GetTriEdges(tID);
        set_triangle_edges(tID, te[0], te[2], te[1]);
    }

	void ReverseOrientation(bool bFlipNormals = true) {
		for ( int tid : TriangleIndices() ) {
			internal_reverse_tri_orientation(tid);
		}
		if ( bFlipNormals && HasVertexNormals() ) {
			for ( int vid : VertexIndices() ) {
				int i = 3*vid;
				normals[i] = -normals[i];
				normals[i+1] = -normals[i+1];
				normals[i+2] = -normals[i+2];
			}
		}
        updateTimeStamp(true);
	}




    /// <summary>
    /// Remove vertex vID, and all connected triangles if bRemoveAllTriangles = true
    /// (if false, them throws exception if there are still any triangles!)
    /// if bPreserveManifold, checks that we will not create a bowtie vertex first
    /// </summary>
    MeshResult RemoveVertex(int vID, bool bRemoveAllTriangles = true, bool bPreserveManifold = false)
    {
        if (vertices_refcount.isValid(vID) == false)
            return MeshResult::Failed_NotAVertex;

        if ( bRemoveAllTriangles ) {

            // if any one-ring vtx is a boundary vtx and one of its outer-ring edges is an
            // interior edge then we will create a bowtie if we remove that triangle
            if ( bPreserveManifold ) {
                for ( int tid : VtxTrianglesItr(vID) ) {
                    Index3i tri = GetTriangle(tid);
                    int j = find_tri_index(vID, tri);
                    int oa = tri[(j + 1) % 3], ob = tri[(j + 2) % 3];
                    int eid = find_edge(oa,ob);
                    if (IsBoundaryEdge(eid))
                        continue;
                    if (IsBoundaryVertex(oa) || IsBoundaryVertex(ob))
                        return MeshResult::Failed_WouldCreateBowtie;
                }
            }

            std::vector<int> tris;
            GetVtxTriangles(vID, tris, true);
            for (int tID : tris) {
                MeshResult result = RemoveTriangle(tID, false, bPreserveManifold);
                if (result != MeshResult::Ok)
                    return result;
            }
        }

		gDevAssert(vertices_refcount.refCount(vID) == 1);
        //if ( vertices_refcount.refCount(vID) != 1)
        //    throw std::exception("DMesh3.RemoveVertex: vertex is still referenced");

        vertices_refcount.decrement(vID);
        gDevAssert(vertices_refcount.isValid(vID) == false);
        vertex_edges.Clear(vID);

        updateTimeStamp(true);
        return MeshResult::Ok;
    }



    /// <summary>
    /// Remove a tID from the mesh. Also removes any unreferenced edges after tri is removed.
    /// If bRemoveIsolatedVertices is false, then if you remove all tris from a vert, that vert is also removed.
    /// If bPreserveManifold, we check that you will not create a bowtie vertex (and return false).
    ///   If this check is not done, you have to make sure you don't create a bowtie, because other
    ///   code assumes we don't have bowties, and will not handle it properly
    /// </summary>
    MeshResult RemoveTriangle(int tID, bool bRemoveIsolatedVertices = true, bool bPreserveManifold = false)
    {
        if ( ! triangles_refcount.isValid(tID) ) {
            gDevAssert(false);
            return MeshResult::Failed_NotATriangle;
        }

        Index3i tv = GetTriangle(tID);
        Index3i te = GetTriEdges(tID);

        // if any tri vtx is a boundary vtx connected to two interior edges, then
        // we cannot remove this triangle because it would create a bowtie vertex!
        // (that vtx already has 2 boundary edges, and we would add two more)
        if (bPreserveManifold) {
            for (int j = 0; j < 3; ++j) {
                if (IsBoundaryVertex(tv[j])) {
                    if (IsBoundaryEdge(te[j]) == false && IsBoundaryEdge(te[(j + 2) % 3]) == false)
                        return MeshResult::Failed_WouldCreateBowtie;
                }
            }
        }

        // Remove triangle from its edges. if edge has no triangles left,
        // then it is removed.
        for (int j = 0; j < 3; ++j) {
            int eid = te[j];
            replace_edge_triangle(eid, tID, InvalidID);
            if (edges[4 * eid + 2] == InvalidID) {
                int a = edges[4 * eid];
                vertex_edges.Remove(a, eid);

                int b = edges[4 * eid + 1];
                vertex_edges.Remove(b, eid);

                edges_refcount.decrement(eid);
            }
        }

        // free this triangle
		triangles_refcount.decrement( tID );
		gDevAssert( triangles_refcount.isValid( tID ) == false );

        // Decrement vertex refcounts. If any hit 1 and we got remove-isolated flag,
        // we need to remove that vertex
        for (int j = 0; j < 3; ++j) {
            int vid = tv[j];
            vertices_refcount.decrement(vid);
            if ( bRemoveIsolatedVertices && vertices_refcount.refCount(vid) == 1) {
                vertices_refcount.decrement(vid);
                gDevAssert(vertices_refcount.isValid(vid) == false);
                vertex_edges.Clear(vid);
            }
        }

        updateTimeStamp(true);
        return MeshResult::Ok;
    }






    virtual MeshResult SetTriangle(int tID, Index3i newv, bool bRemoveIsolatedVertices = true)
    {
        Index3i tv = GetTriangle(tID);
        Index3i te = GetTriEdges(tID);
        if (tv[0] == newv[0] && tv[1] == newv[1])
            te[0] = -1;
        if (tv[1] == newv[1] && tv[2] == newv[2])
            te[1] = -1;
        if (tv[2] == newv[2] && tv[0] == newv[0])
            te[2] = -1;

        if (!triangles_refcount.isValid(tID)) {
            gDevAssert(false);
            return MeshResult::Failed_NotATriangle;
        }
        if (IsVertex(newv[0]) == false || IsVertex(newv[1]) == false || IsVertex(newv[2]) == false) {
            gDevAssert(false);
            return MeshResult::Failed_NotAVertex;
        }
        if (newv[0] == newv[1] || newv[0] == newv[2] || newv[1] == newv[2]) {
            gDevAssert(false);
            return MeshResult::Failed_BrokenTopology;
        }
        // look up edges. if any already have two triangles, this would 
        // create non-manifold geometry and so we do not allow it
        int e0 = find_edge(newv[0], newv[1]);
        int e1 = find_edge(newv[1], newv[2]);
        int e2 = find_edge(newv[2], newv[0]);
        if ((te[0] != -1 && e0 != InvalidID && IsBoundaryEdge(e0) == false)
                || (te[1] != -1 && e1 != InvalidID && IsBoundaryEdge(e1) == false)
                || (te[2] != -1 && e2 != InvalidID && IsBoundaryEdge(e2) == false)) {
            return MeshResult::Failed_BrokenTopology;
        }


        // [TODO] check that we are not going to create invalid stuff...

        // Remove triangle from its edges. if edge has no triangles left, then it is removed.
        for (int j = 0; j < 3; ++j) {
            int eid = te[j];
            if (eid == -1)      // we don't need to modify this edge
                continue;
            replace_edge_triangle(eid, tID, InvalidID);
            if (edges[4 * eid + 2] == InvalidID) {
                int a = edges[4 * eid];
                vertex_edges.Remove(a, eid);

                int b = edges[4 * eid + 1];
                vertex_edges.Remove(b, eid);

                edges_refcount.decrement(eid);
            }
        }

        // Decrement vertex refcounts. If any hit 1 and we got remove-isolated flag,
        // we need to remove that vertex
        for (int j = 0; j < 3; ++j) {
            int vid = tv[j];
            if (vid == newv[j])     // we don't need to modify this vertex
                continue;
            vertices_refcount.decrement(vid);
            if (bRemoveIsolatedVertices && vertices_refcount.refCount(vid) == 1) {
                vertices_refcount.decrement(vid);
                gDevAssert(vertices_refcount.isValid(vid) == false);
                vertex_edges.Clear(vid);
            }
        }


        // ok now re-insert with vertices
        int i = 3 * tID;
        for (int j = 0; j < 3; ++j) {
            if (newv[j] != tv[j]) {
                triangles[i + j] = newv[j];
                vertices_refcount.increment(newv[j]);
            }
        }

        if ( te[0] != -1 )
            add_tri_edge(tID, newv[0], newv[1], 0, e0);
        if (te[1] != -1)
            add_tri_edge(tID, newv[1], newv[2], 1, e1);
        if (te[2] != -1)
            add_tri_edge(tID, newv[2], newv[0], 2, e2);

        updateTimeStamp(true);
        return MeshResult::Ok;
    }







    struct EdgeSplitInfo {
		bool bIsBoundary;
		int vNew;
		int eNewBN;      // edge [vNew,vB] (original was AB)
		int eNewCN;      // edge [vNew,vC] (C is "first" other vtx in ring)
		int eNewDN;		// edge [vNew,vD] (D is "second" other, which doesn't exist on bdry)
        int eNewT2;
        int eNewT3;
	};
	MeshResult SplitEdge(int vA, int vB, EdgeSplitInfo & split)
	{
		int eid = find_edge(vA, vB);
		if ( eid == InvalidID ) {
			split = EdgeSplitInfo();
			return MeshResult::Failed_NotAnEdge;
		}
		return SplitEdge(eid, split);
	}
    /// <summary>
    /// Split edge eab. 
    /// split_t defines position along edge, and is assumed to be based on order of vertices returned by GetEdgeV()
    /// </summary>
	MeshResult SplitEdge(int eab, EdgeSplitInfo & split, double split_t = 0.5)
	{
		split = EdgeSplitInfo();
		if (! IsEdge(eab) )
			return MeshResult::Failed_NotAnEdge;

		// look up primary edge & triangle
		int eab_i = 4*eab;
		int a = edges[eab_i], b = edges[eab_i + 1];
		int t0 = edges[eab_i + 2];
        if (t0 == InvalidID)
            return MeshResult::Failed_BrokenTopology;
		Index3i T0tv = GetTriangle(t0);
		int c = orient_tri_edge_and_find_other_vtx(a, b, T0tv);
        if (vertices_refcount.rawRefCount(c) > 32764)
            return MeshResult::Failed_HitValenceLimit;
        if (a != edges[eab_i])
            split_t = 1.0 - split_t;    // if we flipped a/b order we need to reverse t

        // quite a bit of code is duplicated between boundary and non-boundary case, but it
        //  is too hard to follow later if we factor it out...
        if ( IsBoundaryEdge(eab) ) {

            // create vertex
            Vector3d vNew = Lerp(GetVertex(a), GetVertex(b), split_t);
            int f = AppendVertex(vNew);
            if (HasVertexNormals())
                SetVertexNormal(f, Lerp(GetVertexNormal(a), GetVertexNormal(b), (float)split_t).normalized());
            if (HasVertexColors())
                SetVertexColor(f, Lerp(GetVertexColor(a), GetVertexColor(b), (float)split_t));
            if (HasVertexUVs())
                SetVertexUV(f, Lerp(GetVertexUV(a), GetVertexUV(b), (float)split_t));

            // look up edge bc, which needs to be modified
            Index3i T0te = GetTriEdges(t0);
			int ebc = T0te[ find_edge_index_in_tri(b, c, T0tv) ];

			// rewrite existing triangle
			replace_tri_vertex(t0, b, f);

			// add second triangle
			int t2 = add_triangle_only(f,b,c, InvalidID, InvalidID, InvalidID);
			if ( HasTriangleGroups() )
				triangle_groups.insertAt(triangle_groups[t0], t2);

			// rewrite edge bc, create edge af
			replace_edge_triangle(ebc, t0, t2);
			int eaf = eab; 
			replace_edge_vertex(eaf, b, f);
            vertex_edges.Remove(b, eab);
            vertex_edges.Insert(f, eaf);

			// create edges fb and fc 
			int efb = add_edge(f, b, t2);
			int efc = add_edge(f, c, t0, t2);

			// update triangle edge-nbrs
			replace_triangle_edge(t0, ebc, efc);
			set_triangle_edges(t2, efb, ebc, efc);

			// update vertex refcounts
			vertices_refcount.increment(c);
			vertices_refcount.increment(f, 2);

			split.bIsBoundary = true;
			split.vNew = f;
            split.eNewBN = efb;
			split.eNewCN = efc;
			split.eNewDN = InvalidID;
            split.eNewT2 = t2;
            split.eNewT3 = InvalidID;

			updateTimeStamp(true);
			return MeshResult::Ok;

		} else {		// interior triangle branch
				
			// look up other triangle
			int t1 = edges[eab_i + 3];
			Index3i T1tv = GetTriangle(t1);
			int d = find_tri_other_vtx( a, b, T1tv );
            if (vertices_refcount.rawRefCount(d) > 32764) 
                return MeshResult::Failed_HitValenceLimit;

            // create vertex
            Vector3d vNew = Lerp(GetVertex(a), GetVertex(b), split_t);
            int f = AppendVertex(vNew);
            if (HasVertexNormals())
                SetVertexNormal(f, Lerp(GetVertexNormal(a), GetVertexNormal(b), (float)split_t).normalized());
            if (HasVertexColors())
                SetVertexColor(f, Lerp(GetVertexColor(a), GetVertexColor(b), (float)split_t));
            if (HasVertexUVs())
                SetVertexUV(f, Lerp(GetVertexUV(a), GetVertexUV(b), (float)split_t));

            // look up edges that we are going to need to update
            // [TODO OPT] could use ordering to reduce # of compares here
            Index3i T0te = GetTriEdges(t0);
			int ebc = T0te[find_edge_index_in_tri( b, c, T0tv )];
			Index3i T1te = GetTriEdges(t1);
			int edb = T1te[find_edge_index_in_tri( d, b, T1tv )];

			// rewrite existing triangles
			replace_tri_vertex(t0, b, f);
			replace_tri_vertex(t1, b, f);

			// add two triangles to close holes we just created
			int t2 = add_triangle_only(f,b,c, InvalidID, InvalidID, InvalidID);
			int t3 = add_triangle_only(f, d, b, InvalidID, InvalidID, InvalidID);
			if ( HasTriangleGroups() ) {
				triangle_groups.insertAt(triangle_groups[t0], t2);
				triangle_groups.insertAt(triangle_groups[t1], t3);
			}

			// update the edges we found above, to point to triangles
			replace_edge_triangle(ebc, t0, t2);
			replace_edge_triangle(edb, t1, t3);

			// edge eab became eaf
			int eaf = eab; //Edge * eAF = eAB;
			replace_edge_vertex(eaf, b, f);

            // update a/b/f vertex-edges
            vertex_edges.Remove(b, eab);
            vertex_edges.Insert(f, eaf);

			// create edges connected to f  (also updates vertex-edges)
			int efb = add_edge( f, b, t2, t3 );
			int efc = add_edge( f, c, t0, t2 );
			int edf = add_edge( d, f, t1, t3 );

			// update triangle edge-nbrs
			replace_triangle_edge(t0, ebc, efc);
			replace_triangle_edge(t1, edb, edf);
			set_triangle_edges(t2, efb, ebc, efc);
			set_triangle_edges(t3, edf, edb, efb);

			// update vertex refcounts
			vertices_refcount.increment( c );
			vertices_refcount.increment( d );
			vertices_refcount.increment( f, 4 );

			split.bIsBoundary = false;
			split.vNew = f;
            split.eNewBN = efb;
			split.eNewCN = efc;
			split.eNewDN = edf;
            split.eNewT2 = t2;
            split.eNewT3 = t3;

            updateTimeStamp(true);
			return MeshResult::Ok;
		}

	}






	struct EdgeFlipInfo {
		int eID;
		int v0,v1;
		int ov0,ov1;
		int t0,t1;
	};
	MeshResult FlipEdge(int vA, int vB, EdgeFlipInfo & flip) {
		int eid = find_edge(vA, vB);
		if ( eid == InvalidID ) {
			flip = EdgeFlipInfo();
			return MeshResult::Failed_NotAnEdge;
		}
		return FlipEdge(eid, flip);
	}
	MeshResult FlipEdge(int eab, EdgeFlipInfo & flip) 
	{
		flip = EdgeFlipInfo();
		if (! IsEdge(eab) )
			return MeshResult::Failed_NotAnEdge;
		if ( IsBoundaryEdge(eab) )
			return MeshResult::Failed_IsBoundaryEdge;

		// find oriented edge [a,b], tris t0,t1, and other verts c in t0, d in t1
		int eab_i = 4*eab;
		int a = edges[eab_i], b = edges[eab_i + 1];
		int t0 = edges[eab_i + 2], t1 = edges[eab_i + 3];
		Index3i T0tv = GetTriangle(t0), T1tv = GetTriangle(t1);
		int c = orient_tri_edge_and_find_other_vtx( a, b, T0tv );
		int d = find_tri_other_vtx(a, b, T1tv);
		if ( c == InvalidID || d == InvalidID ) {
			return MeshResult::Failed_BrokenTopology;
		}

		int flipped = find_edge(c,d);
		if ( flipped != InvalidID )
			return MeshResult::Failed_FlippedEdgeExists;

		// find edges bc, ca, ad, db
		int ebc = find_tri_neighbour_edge(t0, b, c);
		int eca = find_tri_neighbour_edge(t0, c, a);
		int ead = find_tri_neighbour_edge(t1,a,d);
		int edb = find_tri_neighbour_edge(t1,d,b);

		// update triangles
		set_triangle(t0, c, d, b);
		set_triangle(t1, d, c, a);

		// update edge AB, which becomes flipped edge CD
		set_edge_vertices(eab, c, d);
		set_edge_triangles(eab, t0,t1);
		int ecd = eab;

		// update the two other edges whose triangle nbrs have changed
		if ( replace_edge_triangle(eca, t0,t1) == -1 )
			throw std::exception("DMesh3.FlipEdge: first replace_edge_triangle failed");
		if ( replace_edge_triangle(edb, t1, t0) == -1 )
			throw std::exception("DMesh3.FlipEdge: second replace_edge_triangle failed");

		// update triangle nbr lists (these are edges)
		set_triangle_edges(t0, ecd, edb, ebc);
		set_triangle_edges(t1, ecd, eca, ead);

		// remove old eab from verts a and b, and decrement ref counts
		if ( vertex_edges.Remove(a, eab) == false ) 
			throw std::exception("DMesh3.FlipEdge: first edge list remove failed");
		if ( vertex_edges.Remove(b, eab) == false ) 
			throw std::exception("DMesh3.FlipEdge: second edge list remove failed");
		vertices_refcount.decrement(a);
		vertices_refcount.decrement(b);
		if ( IsVertex(a) == false || IsVertex(b) == false )
			throw std::exception("DMesh3.FlipEdge: either a or b is not a vertex?");

        // add edge ecd to verts c and d, and increment ref counts
        vertex_edges.Insert(c, ecd);
        vertex_edges.Insert(d, ecd);
		vertices_refcount.increment(c);
		vertices_refcount.increment(d);

		// success! collect up results
		flip.eID = eab;
		flip.v0 = a; flip.v1 = b;
		flip.ov0 = c; flip.ov1 = d;
		flip.t0 = t0; flip.t1 = t1;

		updateTimeStamp(true);
		return MeshResult::Ok;
	}



	void debug_fail(const std::string & s) {
	#if DEBUG
		System.Console.WriteLine("DMesh3.CollapseEdge: check failed: " + s);
		gDevAssert(false);
		//throw Exception("DMesh3.CollapseEdge: check failed: " + s);
	#endif
	}


	void check_tri(int t) {
		Index3i tv = GetTriangle(t);
		if ( tv[0] == tv[1] || tv[0] == tv[2] || tv[1] == tv[2] )
			gDevAssert(false);
	}
	void check_edge(int e) {
		Index2i tv = GetEdgeT(e);
		if ( tv[0] == -1 )
			gDevAssert(false);
	}


	struct EdgeCollapseInfo {
		int vKept;
		int vRemoved;
		bool bIsBoundary;

        int eCollapsed;              // edge we collapsed
        int tRemoved0, tRemoved1;    // tris we removed (second may be invalid)
        int eRemoved0, eRemoved1;    // edges we removed (second may be invalid)
        int eKept0, eKept1;          // edges we kept (second may be invalid)
	};
	MeshResult CollapseEdge(int vKeep, int vRemove, EdgeCollapseInfo & collapse)
	{
		collapse = EdgeCollapseInfo();
				
		if ( IsVertex(vKeep) == false || IsVertex(vRemove) == false )
			return MeshResult::Failed_NotAnEdge;

		int b = vKeep;		// renaming for sanity. We remove a and keep b
		int a = vRemove;

		int eab = find_edge( a, b );
		if (eab == InvalidID)
			return MeshResult::Failed_NotAnEdge;

		int t0 = edges[4*eab+2];
        if (t0 == InvalidID)
            return MeshResult::Failed_BrokenTopology;
		Index3i T0tv = GetTriangle(t0);
		int c = find_tri_other_vtx(a, b, T0tv);

		// look up opposing triangle/vtx if we are not in boundary case
		bool bIsBoundaryEdge = false;
		int d = InvalidID;
		int t1 = edges[4*eab+3];
		if (t1 != InvalidID) {
			Index3i T1tv = GetTriangle(t1);
			d = find_tri_other_vtx( a, b, T1tv );
			if (c == d)
				return MeshResult::Failed_FoundDuplicateTriangle;
		} else {
			bIsBoundaryEdge = true;
		}

		// We cannot collapse if edge lists of a and b share vertices other
		//  than c and d  (because then we will make a triangle [x b b].
		//  Unfortunately I cannot see a way to do this more efficiently than brute-force search
		//  [TODO] if we had tri iterator for a, couldn't we check each tri for b  (skipping t0 and t1) ?
        int edges_a_count = vertex_edges.Count(a); 
        int eac = InvalidID, ead = InvalidID, ebc = InvalidID, ebd = InvalidID;
        for ( int eid_a : vertex_edges.values(a) ) { 
			int vax =  edge_other_v(eid_a, a);
            if ( vax == c ) {
                eac = eid_a;
                continue;
            }
            if ( vax == d ) {
                ead = eid_a;
                continue;
            }
			if ( vax == b )
				continue;
            for (int eid_b : vertex_edges.values(b)) {
				if ( edge_other_v(eid_b, b) == vax )
					return MeshResult::Failed_InvalidNeighbourhood;
			}
		}

        // [RMS] I am not sure this tetrahedron case will detect bowtie vertices.
        // But the single-triangle case does

		// We cannot collapse if we have a tetrahedron. In this case a has 3 nbr edges,
		//  and edge cd exists. But that is not conclusive, we also have to check that
		//  cd is an internal edge, and that each of its tris contain a or b
		if (edges_a_count == 3 && bIsBoundaryEdge == false) {
			int edc = find_edge( d, c );
			int edc_i = 4*edc;
			if (edc != InvalidID && edges[edc_i+3] != InvalidID ) {
				int edc_t0 = edges[edc_i+2];
				int edc_t1 = edges[edc_i+3];

				if ( (tri_has_v(edc_t0,a) && tri_has_v(edc_t1, b)) 
					|| (tri_has_v(edc_t0, b) && tri_has_v(edc_t1, a)) )
				return MeshResult::Failed_CollapseTetrahedron;
			}

		} else if (bIsBoundaryEdge == true && IsBoundaryEdge(eac) ) {
            // Cannot collapse edge if we are down to a single triangle
            ebc = find_edge_from_tri(b, c, t0);
            if ( IsBoundaryEdge(ebc) )
				return MeshResult::Failed_CollapseTriangle;
		}

		// [RMS] this was added from C++ version...seems like maybe I never hit
		//   this case? Conceivably could be porting bug but looking at the
		//   C++ code I cannot see how we could possibly have caught this case...
		//
		// cannot collapse an edge where both vertices are boundary vertices
		// because that would create a bowtie
		//
		// NOTE: potentially scanning all edges here...couldn't we
		//  pick up eac/bc/ad/bd as we go? somehow?
		if ( bIsBoundaryEdge == false && IsBoundaryVertex(a) && IsBoundaryVertex(b) )
			return MeshResult::Failed_InvalidNeighbourhood;


		// 1) remove edge ab from vtx b
		// 2) find edges ad and ac, and tris tad, tac across those edges  (will use later)
		// 3) for other edges, replace a with b, and add that edge to b
		// 4) replace a with b in all triangles connected to a
		int tad = InvalidID, tac = InvalidID;
        for ( int eid : vertex_edges.values(a)) { 
			int o = edge_other_v(eid, a);
			if (o == b) {
				if ( vertex_edges.Remove(b,eid) != true )
					debug_fail("remove case o == b");
			} else if (o == c) {
				if ( vertex_edges.Remove(c, eid) != true )
					debug_fail("remove case o == c");
				tac = edge_other_t(eid, t0);
			} else if (o == d) {
				if (vertex_edges.Remove(d, eid) != true )
					debug_fail("remove case o == c, step 1");
				tad = edge_other_t(eid, t1);
			} else {
				if ( replace_edge_vertex(eid, a, b) == -1 )
					debug_fail("remove case else");
                vertex_edges.Insert(b, eid);
			}

			// [TODO] perhaps we can already have unique tri list because of the manifold-nbrhood check we need to do...
			for (int j = 0; j < 2; ++j) {
				int t_j = edges[4*eid + 2 + j];
				if (t_j != InvalidID && t_j != t0 && t_j != t1) {
					if ( tri_has_v(t_j, a) ) {
						if ( replace_tri_vertex(t_j, a, b) == -1 )
							debug_fail("remove last check");
						vertices_refcount.increment(b);
						vertices_refcount.decrement(a);
					}
				}
			}
		}

		if (bIsBoundaryEdge == false) {

            // remove all edges from vtx a, then remove vtx a
            vertex_edges.Clear(a);
			gDevAssert( vertices_refcount.refCount(a) == 3 );		// in t0,t1, and initial ref
			vertices_refcount.decrement( a, 3 );
			gDevAssert( vertices_refcount.isValid( a ) == false );

			// remove triangles T0 and T1, and update b/c/d refcounts
			triangles_refcount.decrement( t0 );
			triangles_refcount.decrement( t1 );
			vertices_refcount.decrement( c );
			vertices_refcount.decrement( d );
			vertices_refcount.decrement( b, 2 );
			gDevAssert( triangles_refcount.isValid( t0 ) == false );
			gDevAssert( triangles_refcount.isValid( t1 ) == false );

			// remove edges ead, eab, eac
			edges_refcount.decrement( ead );
			edges_refcount.decrement( eab );
			edges_refcount.decrement( eac );
			gDevAssert( edges_refcount.isValid( ead ) == false );
			gDevAssert( edges_refcount.isValid( eab ) == false );
			gDevAssert( edges_refcount.isValid( eac ) == false );

			// replace t0 and t1 in edges ebd and ebc that we kept
			ebd = find_edge_from_tri( b, d, t1 );
            if ( ebc == InvalidID )   // we may have already looked this up
				ebc = find_edge_from_tri( b, c, t0 );

			if( replace_edge_triangle(ebd, t1, tad ) == -1 )
				debug_fail("isboundary=false branch, ebd replace triangle");

			if ( replace_edge_triangle(ebc, t0, tac ) == -1 )
				debug_fail("isboundary=false branch, ebc replace triangle");

			// update tri-edge-nbrs in tad and tac
			if (tad != InvalidID) {
				if ( replace_triangle_edge(tad, ead, ebd ) == -1 )
					debug_fail("isboundary=false branch, ebd replace triangle");
			}
			if (tac != InvalidID) {
				if ( replace_triangle_edge(tac, eac, ebc ) == -1 )
					debug_fail("isboundary=false branch, ebd replace triangle");
			}

		} else {

            //  this is basically same code as above, just not referencing t1/d

            // remove all edges from vtx a, then remove vtx a
            vertex_edges.Clear(a);
            gDevAssert( vertices_refcount.refCount( a ) == 2 );		// in t0 and initial ref
			vertices_refcount.decrement( a, 2 );
			gDevAssert( vertices_refcount.isValid( a ) == false );

			// remove triangle T0 and update b/c refcounts
			triangles_refcount.decrement( t0 );
			vertices_refcount.decrement( c );
			vertices_refcount.decrement( b );
			gDevAssert( triangles_refcount.isValid( t0 ) == false );

			// remove edges eab and eac
			edges_refcount.decrement( eab );
			edges_refcount.decrement( eac );
			gDevAssert( edges_refcount.isValid( eab ) == false );
			gDevAssert( edges_refcount.isValid( eac ) == false );

			// replace t0 in edge ebc that we kept
			ebc = find_edge_from_tri( b, c, t0 );
			if ( replace_edge_triangle(ebc, t0, tac ) == -1 )
				debug_fail("isboundary=false branch, ebc replace triangle");

			// update tri-edge-nbrs in tac
			if (tac != InvalidID) {
				if ( replace_triangle_edge(tac, eac, ebc ) == -1 )
					debug_fail("isboundary=true branch, ebd replace triangle");
			}
		}

		collapse.vKept = vKeep;
		collapse.vRemoved = vRemove;
		collapse.bIsBoundary = bIsBoundaryEdge;
        collapse.eCollapsed = eab;
        collapse.tRemoved0 = t0; collapse.tRemoved1 = t1;
        collapse.eRemoved0 = eac; collapse.eRemoved1 = ead;
        collapse.eKept0 = ebc; collapse.eKept1 = ebd;

		updateTimeStamp(true);
		return MeshResult::Ok;
	}





	struct MergeEdgesInfo
	{
		int eKept;
		int eRemoved;

		Vector2i vKept;
		Vector2i vRemoved;           // either may be InvalidID if it was same as vKept

		Vector2i eRemovedExtra;      // bonus collapsed edge, or InvalidID
		Vector2i eKeptExtra;			// edge paired w/ eRemovedExtra
	};
	MeshResult MergeEdges(int eKeep, int eDiscard, MergeEdgesInfo & merge_info) 
	{
		merge_info = MergeEdgesInfo();
		if (IsEdge(eKeep) == false || IsEdge(eDiscard) == false)
			return MeshResult::Failed_NotAnEdge;

		Index4i edgeinfo_keep = GetEdge(eKeep);
		Index4i edgeinfo_discard = GetEdge(eDiscard);
		if (edgeinfo_keep[3] != InvalidID || edgeinfo_discard[3] != InvalidID)
			return MeshResult::Failed_NotABoundaryEdge;

		int a = edgeinfo_keep[0], b = edgeinfo_keep[1];
		int tab = edgeinfo_keep[2];
		int eab = eKeep;
		int c = edgeinfo_discard[0], d = edgeinfo_discard[1];
		int tcd = edgeinfo_discard[2];
		int ecd = eDiscard;

		// Need to correctly orient a,b and c,d and then check that 
		// we will not join triangles with incompatible winding order
		// I can't see how to do this purely topologically. 
		// So relying on closest-pairs testing.
		orient_tri_edge(a, b, GetTriangle(tab));
		//int tcd_otherv = orient_tri_edge_and_find_other_vtx(ref c, ref d, GetTriangle(tcd));
		orient_tri_edge(c, d, GetTriangle(tcd));
		int x = c; c = d; d = x;   // joinable bdry edges have opposing orientations, so flip to get ac and b/d correspondences
		Vector3d Va = GetVertex(a), Vb = GetVertex(b), Vc = GetVertex(c), Vd = GetVertex(d);
		if ( ((Va-Vc).squaredNorm() + (Vb-Vd).squaredNorm()) >
			 ((Va-Vd).squaredNorm() + (Vb-Vc).squaredNorm()) )
			return MeshResult::Failed_SameOrientation;

		// alternative that detects normal flip of triangle tcd. This is a more 
		// robust geometric test, but fails if tri is degenerate...also more expensive
		//Vector3d otherv = GetVertex(tcd_otherv);
		//Vector3d Ncd = MathUtil.FastNormalDirection(GetVertex(c), GetVertex(d), otherv);
		//Vector3d Nab = MathUtil.FastNormalDirection(GetVertex(a), GetVertex(b), otherv);
		//if (Ncd.Dot(Nab) < 0)
		//return MeshResult::Failed_SameOrientation;

		merge_info.eKept = eab;
		merge_info.eRemoved = ecd;

        // if a/c or b/d are connected by an existing edge, we can't merge
        if (a != c && find_edge(a,c) != InvalidID )
            return MeshResult::Failed_InvalidNeighbourhood;
        if (b != d && find_edge(b, d) != InvalidID)
            return MeshResult::Failed_InvalidNeighbourhood;

        // if vertices at either end already share a common neighbour vertex, and we 
        // do the merge, that would create duplicate edges. This is something like the
        // 'link condition' in edge collapses. 
        // Note that we have to catch cases where both edges to the shared vertex are
        // boundary edges, in that case we will also merge this edge later on
        if ( a != c ) {
            int ea = 0, ec = 0, other_v = (b == d) ? b : -1;
            for ( int cnbr : VtxVerticesItr(c) ) {
                if (cnbr != other_v && (ea = find_edge(a, cnbr)) != InvalidID) {
                    ec = find_edge(c, cnbr);
                    if (IsBoundaryEdge(ea) == false || IsBoundaryEdge(ec) == false)
                        return MeshResult::Failed_InvalidNeighbourhood;
                }
            }
        }
        if ( b != d ) {
            int eb = 0, ed = 0, other_v = (a == c) ? a : -1;
            for ( int dnbr : VtxVerticesItr(d)) {
                if (dnbr != other_v && (eb = find_edge(b, dnbr)) != InvalidID) {
                    ed = find_edge(d, dnbr);
                    if (IsBoundaryEdge(eb) == false || IsBoundaryEdge(ed) == false)
                        return MeshResult::Failed_InvalidNeighbourhood;
                }
            }
        }
            

        // [TODO] this acts on each interior tri twice. could avoid using vtx-tri iterator?
        if (a != c) {
			// replace c w/ a in edges and tris connected to c, and move edges to a
            for ( int eid : vertex_edges.values(c)) { 
				if (eid == eDiscard)
					continue;
				replace_edge_vertex(eid, c, a);
				short rc = 0;
				if (replace_tri_vertex(edges[4 * eid + 2], c, a) >= 0)
					rc++;
				if (edges[4 * eid + 3] != InvalidID) {
					if (replace_tri_vertex(edges[4 * eid + 3], c, a) >= 0)
						rc++;
				}
                vertex_edges.Insert(a, eid);
				if (rc > 0) {
					vertices_refcount.increment(a, rc);
					vertices_refcount.decrement(c, rc);
				}
			}
            vertex_edges.Clear(c);
			vertices_refcount.decrement(c);
			merge_info.vRemoved[0] = c;
		} else {
            vertex_edges.Remove(a, ecd);
			merge_info.vRemoved[0] = InvalidID;
		}
		merge_info.vKept[0] = a;

		if (d != b) {
			// replace d w/ b in edges and tris connected to d, and move edges to b
            for (int eid : vertex_edges.values(d)) { 
				if (eid == eDiscard)
					continue;
				replace_edge_vertex(eid, d, b);
				short rc = 0;
				if (replace_tri_vertex(edges[4 * eid + 2], d, b) >= 0)
					rc++;
				if (edges[4 * eid + 3] != InvalidID) {
					if (replace_tri_vertex(edges[4 * eid + 3], d, b) >= 0)
						rc++;
				}
                vertex_edges.Insert(b, eid);
				if (rc > 0) {
					vertices_refcount.increment(b, rc);
					vertices_refcount.decrement(d, rc);
				}

			}
            vertex_edges.Clear(d);
			vertices_refcount.decrement(d);
			merge_info.vRemoved[1] = d;
		} else {
            vertex_edges.Remove(b, ecd);
			merge_info.vRemoved[1] = InvalidID;
		}
		merge_info.vKept[1] = b;

		// replace edge cd with edge ab in triangle tcd
		replace_triangle_edge(tcd, ecd, eab);
		edges_refcount.decrement(ecd);

		// update edge-tri adjacency
		set_edge_triangles(eab, tab, tcd);

		// Once we merge ab to cd, there may be additional edges (now) connected
		// to either a or b that are connected to the same vertex on their 'other' side.
		// So we now have two boundary edges connecting the same two vertices - disaster!
		// We need to find and merge these edges. 
		// Q: I don't think it is possible to have multiple such edge-pairs at a or b
		//    But I am not certain...is a bit tricky to handle because we modify edges_v...
		merge_info.eRemovedExtra = Vector2i(InvalidID, InvalidID);
		merge_info.eKeptExtra = merge_info.eRemovedExtra;
		for (int vi = 0; vi < 2; ++vi) {
			int v1 = a, v2 = c;   // vertices of merged edge
			if ( vi == 1 ) {
				v1 = b; v2 = d;
			}
			if (v1 == v2)
				continue;
            std::vector<int> edges_v = vertex_edges_list(v1);
            int Nedges = (int)edges_v.size();
			bool found = false;
			// in this loop, we compare 'other' vert_1 and vert_2 of edges around v1.
			// problem case is when vert_1 == vert_2  (ie two edges w/ same other vtx).
            //restart_merge_loop:
			for (int i = 0; i < Nedges && found == false; ++i) {
				int edge_1 = edges_v[i];
				if ( IsBoundaryEdge(edge_1) == false)
					continue;
				int vert_1 = edge_other_v(edge_1, v1);
				for (int j = i + 1; j < Nedges; ++j) {
					int edge_2 = edges_v[j];
					int vert_2 = edge_other_v(edge_2, v1);
					if (vert_1 == vert_2 && IsBoundaryEdge(edge_2)) { // if ! boundary here, we are in deep trouble...
						// replace edge_2 w/ edge_1 in tri, update edge and vtx-edge-nbr lists
						int tri_1 = edges[4 * edge_1 + 2];
						int tri_2 = edges[4 * edge_2 + 2];
						replace_triangle_edge(tri_2, edge_2, edge_1);
						set_edge_triangles(edge_1, tri_1, tri_2);
                        vertex_edges.Remove(v1, edge_2);
                        vertex_edges.Remove(vert_1, edge_2);
						edges_refcount.decrement(edge_2);
						merge_info.eRemovedExtra[vi] = edge_2;
						merge_info.eKeptExtra[vi] = edge_1;

                        //edges_v = vertex_edges_list(v1);      // this code allows us to continue checking, ie in case we had
                        //Nedges = edges_v.Count;               // multiple such edges. but I don't think it's possible.
                        //goto restart_merge_loop;
                        found = true;			  // exit outer i loop
                        break;					  // exit inner j loop
                    }
				}
			}
		}

        updateTimeStamp(true);
        return MeshResult::Ok;
	}







    struct PokeTriangleInfo
    {
        int new_vid;
        int new_t1, new_t2;
        Index3i new_edges;
	};
    virtual MeshResult PokeTriangle(int tid, PokeTriangleInfo & result)
    {
        return PokeTriangle(tid, Vector3d::Ones() / 3.0, result);
    }
    virtual MeshResult PokeTriangle(int tid, const Vector3d & baryCoordinates, PokeTriangleInfo & result)
    {
        result = PokeTriangleInfo();

        if (!IsTriangle(tid))
            return MeshResult::Failed_NotATriangle;

        Index3i tv = GetTriangle(tid);
        Index3i te = GetTriEdges(tid);

        // create vertex with interpolated vertex attribs
        NewVertexInfo vinfo;
        GetTriBaryPoint(tid, baryCoordinates[0], baryCoordinates[1], baryCoordinates[2], vinfo);
        int center = AppendVertex(vinfo);

        // add in edges to center vtx, do not connect to triangles yet
        int eaC = add_edge(tv[0], center, -1, -1);
        int ebC = add_edge(tv[1], center, -1, -1);
        int ecC = add_edge(tv[2], center, -1, -1);
        vertices_refcount.increment(tv[0]);
        vertices_refcount.increment(tv[1]);
        vertices_refcount.increment(tv[2]);
        vertices_refcount.increment(center, 3);

        // old triangle becomes tri along first edge
        set_triangle(tid, tv[0], tv[1], center);
        set_triangle_edges(tid, te[0], ebC, eaC);

        // add two triangles
        int t1 = add_triangle_only(tv[1], tv[2], center, te[1], ecC, ebC );
        int t2 = add_triangle_only(tv[2], tv[0], center, te[2], eaC, ecC);

        // second and third edges of original tri have neighbours
        replace_edge_triangle(te[1], tid, t1);
        replace_edge_triangle(te[2], tid, t2);

        // set the triangles for the edges we created above
        set_edge_triangles(eaC, tid, t2);
        set_edge_triangles(ebC, tid, t1);
        set_edge_triangles(ecC, t1, t2);

        // transfer groups
        if ( HasTriangleGroups() ) {
            int g = triangle_groups[tid];
            triangle_groups.insertAt(g, t1);
            triangle_groups.insertAt(g, t2);
        }

        result.new_vid = center;
        result.new_t1 = t1;
        result.new_t2 = t2;
        result.new_edges = Index3i(eaC, ebC, ecC);

        updateTimeStamp(true);
        return MeshResult::Ok;
    }



	std::string MeshInfoString()
	{
		std::ostringstream b;
		b << "Vertices  " << " count " << VertexCount() << " max " << MaxVertexID() << " " << vertices_refcount.UsageStats() << std::endl;
		b << "Triangles " << " count " << TriangleCount() << " max " << MaxTriangleID() << " " << triangles_refcount.UsageStats() << std::endl;
		b << "Edges     " << " count " << EdgeCount() << " max " << MaxEdgeID() << " " << edges_refcount.UsageStats() << std::endl;
		b << "Normals " << HasVertexNormals() << "  Colors " << HasVertexColors() << "  UVs " << HasVertexUVs() << "  Groups " << HasTriangleGroups() << std::endl;
		b << "Closed " << CachedIsClosed() << " Compact " << IsCompact() << " timestamp " << timestamp << " shape_timestamp " << shape_timestamp << "  MaxGroupID " << MaxGroupID() << std::endl;
		b << "VertexEdges " << vertex_edges.MemoryUsage() << std::endl;
		return b.str();
	}

    
	/// <summary>
	/// Check if this m2 is the same as this mesh. By default only checks
	/// vertices and triangles, turn on other parameters w/ flags
	/// </summary>
	bool IsSameMesh(const DMesh3 & m2, bool bCheckConnectivity, bool bCheckEdgeIDs = false,
		bool bCheckNormals = false, bool bCheckColors = false, bool bCheckUVs = false,
		bool bCheckGroups = false,
		float Epsilon = Wml::Mathf::EPSILON)
	{
		if (VertexCount() != m2.VertexCount())
			return false;
		if (TriangleCount() != m2.TriangleCount())
			return false;
		for (int vid : VertexIndices()) {
			if (m2.IsVertex(vid) == false || EpsilonEqual(GetVertex(vid), m2.GetVertex(vid), Epsilon) == false)
				return false;
		}
		for (int tid : TriangleIndices()) {
			if (m2.IsTriangle(tid) == false || (GetTriangle(tid) != m2.GetTriangle(tid)) )
				return false;
		}
		if (bCheckConnectivity) {
			for (int eid : EdgeIndices()) {
				Index4i e = GetEdge(eid);
				int other_eid = m2.FindEdge(e[0], e[1]);
				if (other_eid == InvalidID)
					return false;
				Index4i oe = m2.GetEdge(other_eid);
				if (std::min(e[2], e[3]) != std::min(oe[2], oe[3]) || std::max(e[2], e[3]) != std::max(oe[2], oe[3]))
					return false;
			}
		}
		if (bCheckEdgeIDs) {
			if (EdgeCount() != m2.EdgeCount())
				return false;
			for (int eid : EdgeIndices()) {
				if (m2.IsEdge(eid) == false || GetEdge(eid) != m2.GetEdge(eid) )
					return false;
			}
		}
		if (bCheckNormals) {
			if (HasVertexNormals() != m2.HasVertexNormals())
				return false;
			if (HasVertexNormals()) {
				for (int vid : VertexIndices()) {
					if ( EpsilonEqual(GetVertexNormal(vid), m2.GetVertexNormal(vid), Epsilon) == false)
						return false;
				}
			}
		}
		if (bCheckColors) {
			if (HasVertexColors() != m2.HasVertexColors())
				return false;
			if (HasVertexColors()) {
				for (int vid : VertexIndices()) {
					if ( EpsilonEqual( GetVertexColor(vid), m2.GetVertexColor(vid), Epsilon) == false)
						return false;
				}
			}
		}
		if (bCheckUVs) {
			if (HasVertexUVs() != m2.HasVertexUVs())
				return false;
			if (HasVertexUVs()) {
				for (int vid : VertexIndices()) {
					if ( EpsilonEqual( GetVertexUV(vid), m2.GetVertexUV(vid), Epsilon) == false)
						return false;
				}
			}
		}
		if (bCheckGroups) {
			if (HasTriangleGroups() != m2.HasTriangleGroups())
				return false;
			if (HasTriangleGroups()) {
				for (int tid : TriangleIndices()) {
					if (GetTriangleGroup(tid) != m2.GetTriangleGroup(tid))
						return false;
				}
			}
		}
		return true;
	}



	
	enum class FailMode { DevAssert, Throw, ReturnOnly };


	/// <summary>
	// This function checks that the mesh is well-formed, ie all internal data
	// structures are consistent
	/// </summary>
	bool CheckValidity(bool bAllowNonManifoldVertices = false, FailMode eFailMode = FailMode::Throw) const {

		std::vector<int> triToVtxRefs;
		triToVtxRefs.resize(MaxVertexID());
		//int[] triToVtxRefs = new int[this.MaxVertexID];

		bool is_ok = true;
		std::function<void(bool)> CheckOrFailF = [&](bool b) { 
			is_ok = is_ok && b; 
		};
		if (eFailMode == FailMode::DevAssert) {
			CheckOrFailF = [&](bool b) {
				gDevAssert(b);
				is_ok = is_ok && b;
			};
		} else if (eFailMode == FailMode::Throw) {
			CheckOrFailF = [&](bool b) {
				if (b == false)
					throw std::exception("DMesh3.CheckValidity: check failed!");
				is_ok = is_ok && b;
			};
		}

		//if (normals != null)
		//	CheckOrFailF(normals.size == vertices.size);
		//if (colors != null)
		//	CheckOrFailF(colors.size == vertices.size);
		//if (uv != null)
		//	CheckOrFailF(uv.size / 2 == vertices.size / 3);
		//if (triangle_groups != null)
		//	CheckOrFailF(triangle_groups.size == triangles.size / 3);

		for (int tID : TriangleIndices()) {

			CheckOrFailF(IsTriangle(tID));
			CheckOrFailF(triangles_refcount.refCount(tID) == 1);

			// vertices must exist
			Index3i tv = GetTriangle(tID);
			for (int j = 0; j < 3; ++j) {
				CheckOrFailF(IsVertex(tv[j]));
				triToVtxRefs[tv[j]] += 1;
			}

			// edges must exist and reference this tri
			Index3i e;
			for (int j = 0; j < 3; ++j) {
				int a = tv[j], b = tv[(j + 1) % 3];
				e[j] = FindEdge(a, b);
				CheckOrFailF(e[j] != InvalidID);
				CheckOrFailF(edge_has_t(e[j], tID));
				CheckOrFailF(e[j] == FindEdgeFromTri(a, b, tID));
			}
			CheckOrFailF(e[0] != e[1] && e[0] != e[2] && e[1] != e[2]);

			// tri nbrs must exist and reference this tri, or same edge must be boundary edge
			Index3i te = GetTriEdges(tID);
			for (int j = 0; j < 3; ++j) {
				int eid = te[j];
				CheckOrFailF(IsEdge(eid));
				int tOther = edge_other_t(eid, tID);
				if (tOther == InvalidID) {
					CheckOrFailF(IsBoundaryTriangle(tID));
					continue;
				}

				CheckOrFailF(tri_has_neighbour_t(tOther, tID) == true);

				// edge must have same two verts as tri for same index
				int a = tv[j], b = tv[(j + 1) % 3];
				Index2i ev = GetEdgeV(te[j]);
				CheckOrFailF(same_pair_unordered(a, b, ev[0], ev[1]));

				// also check that nbr edge has opposite orientation
				Index3i othertv = GetTriangle(tOther);
				int found = find_tri_ordered_edge(b, a, othertv);
				CheckOrFailF(found != InvalidID);
			}
		}


		// edge verts/tris must exist
		for (int eID : EdgeIndices()) {
			CheckOrFailF(IsEdge(eID));
			CheckOrFailF(edges_refcount.refCount(eID) == 1);
			Index2i ev = GetEdgeV(eID);
			Index2i et = GetEdgeT(eID);
			CheckOrFailF(IsVertex(ev[0]));
			CheckOrFailF(IsVertex(ev[1]));
			CheckOrFailF(et[0] != InvalidID);
			CheckOrFailF(ev[0] < ev[1]);
			CheckOrFailF(IsTriangle(et[0]));
			if (et[1] != InvalidID) {
				CheckOrFailF(IsTriangle(et[1]));
			}
		}

		// verify compact check
		bool is_compact = vertices_refcount.is_dense();
		if (is_compact) {
			for (int vid = 0; vid < vertices.length() / 3; ++vid) {
				CheckOrFailF(vertices_refcount.isValid(vid));
			}
		}

		// vertex edges must exist and reference this vert
		for (int vID : VertexIndices()) {
			CheckOrFailF(IsVertex(vID));

			Vector3d v = GetVertex(vID);
			CheckOrFailF( _isnan(v.squaredNorm()) == false);
			CheckOrFailF( isfinite(v.squaredNorm()) );

			for (int edgeid : vertex_edges.values(vID)) {
				CheckOrFailF(IsEdge(edgeid));
				CheckOrFailF(edge_has_v(edgeid, vID));

				int otherV = edge_other_v(edgeid, vID);
				int e2 = find_edge(vID, otherV);
				CheckOrFailF(e2 != InvalidID);
				CheckOrFailF(e2 == edgeid);
				e2 = find_edge(otherV, vID);
				CheckOrFailF(e2 != InvalidID);
				CheckOrFailF(e2 == edgeid);
			}

			for (int nbr_vid : VtxVerticesItr(vID)) {
				CheckOrFailF(IsVertex(nbr_vid));
				int edge = find_edge(vID, nbr_vid);
				CheckOrFailF(IsEdge(edge));
			}

			std::vector<int> vTris, vTris2;
			GetVtxTriangles(vID, vTris, false);
			GetVtxTriangles(vID, vTris2, true);
			CheckOrFailF(vTris.size() == vTris2.size());
			//System.Console.WriteLine(string.Format("{0} {1} {2}", vID, vTris.Count, GetVtxEdges(vID).Count));
			if (bAllowNonManifoldVertices)
				CheckOrFailF(vTris.size() <= GetVtxEdgeCount(vID));
			else
				CheckOrFailF(vTris.size() == GetVtxEdgeCount(vID) || vTris.size() == GetVtxEdgeCount(vID) - 1);
			CheckOrFailF(vertices_refcount.refCount(vID) == vTris.size() + 1);
			CheckOrFailF(triToVtxRefs[vID] == vTris.size());
			for (int tID : vTris) {
				CheckOrFailF(tri_has_v(tID, vID));
			}

			// check that edges around vert only references tris above, and reference all of them!
			std::vector<int> vRemoveTris(vTris);
			for (int edgeid : vertex_edges.values(vID)) {
				Index2i edget = GetEdgeT(edgeid);
				CheckOrFailF( Contains(vTris,edget[0]) );
				if (edget[1] != InvalidID)
					CheckOrFailF( Contains(vTris, edget[1]) );
				Remove(vRemoveTris, edget[0]);
				if (edget[1] != InvalidID)
					Remove(vRemoveTris, edget[1]);
			}
			CheckOrFailF(vRemoveTris.size() == 0);
		}

		return is_ok;
	}







};   // end DMesh3


}