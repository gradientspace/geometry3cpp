#pragma once

#include <g3types.h>
#include <GeometryInterfaces.h>
#include <rcvector.h>
#include <object_pool.h>
#include <fixed_list.h>
#include <XDMesh3.h>


namespace g3 {


enum MeshConfig 
{
	VertexNormals       = 0b00001,
	VertexColors        = 0b00010,
	VertexUVs           = 0b00100,
	TriangleColors      = 0b01000,

	VertexNormalsColors = 0b00011
};
typedef unsigned int MeshConfigFlags;


struct EdgeFlipInfo
{
    EdgeID eID;
    Vector2i v;
    Vector2i ov;
    Vector2i t;
};

struct EdgeSplitInfo
{
	VertexID vNew;
	bool bIsBoundary;
};

struct EdgeCollapseInfo
{
	VertexID vKept;
	VertexID vRemoved;
	bool bIsBoundary;
};

template<typename Real, template<typename> class VectorType = dvector >
class OldDMesh3 : public IDynamicMesh<Real>
{
public:
	typedef OldDMesh3<Real, VectorType> DMesh3T;
	static Vector3<Real> InvalidVertex;
	static Vector3i InvalidTriangle;
	static Vector2i InvalidEdge;

	virtual ~OldDMesh3();

	OldDMesh3( MeshConfigFlags flags = 0 );
	OldDMesh3(const DMesh3T & copy);
	OldDMesh3(const DMesh3T && move);

	const DMesh3T & operator=( const DMesh3T & copy );
	const DMesh3T & operator=( const DMesh3T && moved );

	//
	// configuration
	//

	bool GetEnableVertexNormals() const		{ return ( m_flags & VertexNormals) != 0; }
	bool GetEnableVertexColors() const		{ return ( m_flags & VertexColors) != 0; }
	bool GetEnableVertexUVs() const			{ return ( m_flags & VertexUVs) != 0; }
	bool GetEnableTriangleColors() const	{ return ( m_flags & TriangleColors) != 0; }

	void SetEnableVertexNormals(bool bEnable);
	void SetEnableVertexColors(bool bEnable);
	void SetEnableVertexUVs(bool bEnable);
	void SetEnableTriangleColors(bool bEnable);


	virtual unsigned int GetVertexCount() const override;
	virtual unsigned int GetTriangleCount() const override;
	virtual unsigned int GetEdgeCount() const;

	virtual VertexID GetMaxVertexID() const;
	virtual TriangleID GetMaxTriangleID() const;
	virtual EdgeID GetMaxEdgeID() const;

	virtual bool IsVertex( VertexID vID ) const;
	virtual bool IsTriangle( TriangleID tID ) const;
	virtual bool IsEdge( EdgeID eID ) const;


	//
	// construction
	//

	virtual VertexID AppendVertex( const Vector3<Real> & v ) override;
	virtual TriangleID AppendTriangle( const Vector3i & t, GroupID gID = 0 ) override;

	MeshResult RemoveTriangle(TriangleID tID, bool bDeleteUnrefVertices = true);
    MeshResult ReverseTriOrientation(TriangleID tID);

	MeshResult FlipEdge( VertexID vA, VertexID vB, EdgeFlipInfo & flip );
	MeshResult FlipEdge( EdgeID eID, EdgeFlipInfo & flip );

	MeshResult SplitEdge( VertexID vA, VertexID vB, EdgeSplitInfo & split );
	MeshResult SplitEdge( EdgeID eID, EdgeSplitInfo & split );

	MeshResult CollapseEdge( VertexID vKeep, VertexID vRemove, EdgeCollapseInfo & collapse );


	virtual const Vector3<Real> & GetVertex( VertexID vID ) const;
	virtual const Vector3i & GetTriangle( TriangleID tID ) const;
	virtual const Vector2i & GetEdgeV( EdgeID eID ) const;
	virtual Vector2i GetEdgeOpposingV( EdgeID eID ) const;
	virtual const Vector2i & GetEdgeT( EdgeID eID ) const;

	virtual void SetVertex( VertexID vID, const Vector3<Real> & v );

	EdgeID FindEdge(VertexID vA, VertexID vB) const;
	MeshResult GetVtxEdges(VertexID vID, std::vector<int> & vEdges) const;
    MeshResult GetVtxTriangles(VertexID vID, std::vector<TriangleID> & vTriangles, bool bUseOrientation = true) const;
    Vector3i GetTriTriangles(TriangleID tID) const;
    inline Vector3i GetTriEdges(TriangleID tID) const;
    
    MeshResult GetVtxOrderedTriangles(VertexID vID, std::vector<TriangleID> & vTriangles, bool bIsInteriorHint = false) const;
    
protected:
	unsigned int m_flags;

	struct Vertex {
		Vector3<Real> v;

		unsigned short bits;
		short ref;
		inline Vector3<Real> printable() const { return v; }
		inline int get_refcount() const { return (int)ref; }
		inline void set_refcount(int i) { ref = (short)i; }
	};
	rcvector<Vertex, VectorType> m_vVertices;
	VectorType<Vector3f> m_vNormals;
	VectorType<Color4b> m_vColors;
	VectorType<Vector2f> m_vUVs;


	struct Triangle {
		Vector3i tv;			// triangle indices
		Vector3i te;			// triangle edges

		struct {
			unsigned int gid : 24;	// max 16 million groups
			unsigned int bits : 6;	// can trade bits here for group id
			int ref : 2;			// range of this field is [-2,2]
		} data;

		inline Vector3i printable() const { return tv; }
		inline int get_refcount() const { return data.ref; }
		inline void set_refcount( int i ) { data.ref = (i > 1) ? 1 : ((i < -1) ? -1 : i); }

        inline bool hasV(VertexID vID) const { return tv[0] == vID || tv[1] == vID || tv[2] == vID; }
        bool hasSequentialV(VertexID vA, VertexID vB) const;
		int replaceV(VertexID vOld, VertexID vNew);
		inline bool hasNeighbourE(EdgeID eID) const { return te[0] == eID || te[1] == eID || te[2] == eID; }
        int findNeighbourIndex(VertexID vA, VertexID vB) const;
        EdgeID findNeighbourEdge(VertexID vA, VertexID vB) const;
        int replaceE(EdgeID eOld, EdgeID eNew);
    };
	rcvector<Triangle, VectorType> m_vTriangles;
	VectorType<Color4b> m_vTriColors;

	struct Edge {
		Vector2i v;		// vertices of edge
		Vector2i t;		// triangles on either side of edge (second element may be InvalidID if boundary edge)

		unsigned short bits;
		short ref;
		inline Vector4i printable() const { return Vector4i(v[0],v[1],t[0],t[1]); }
		inline int get_refcount() const { return (int)ref; }
		inline void set_refcount(int i) { ref = (short)i; }

		inline bool hasV(VertexID vID) const { return v[0] == vID || v[1] == vID; }
        inline VertexID otherV(VertexID vID) const { return (v[0] == vID) ? v[1] : ((v[1] == vID) ? v[0] : InvalidID); }
		inline bool hasT(TriangleID tID) const { return t[0] == tID || t[1] == tID; }
        inline TriangleID otherT(TriangleID tID) const { return (t[0] == tID) ? t[1] : ((t[1] == tID) ? t[0] : InvalidID); }
		inline bool isBoundary() const { return t[1] == InvalidID; }
        inline int replaceTriangle(TriangleID tOld, TriangleID tNew);   // returns index, or -1
		inline int replaceVertex(VertexID vOld, VertexID vNew);   // returns index, or -1
	};
	rcvector<Edge, VectorType> m_vEdges;

	// idea here is that we have a list for every /used/ vertex anyway, might
	// as well just store them in another vector!
	// trade-off is that the basic list will still exist for 'free' vertex IDs.
	// this could be significant...but hopefully we make it up on all the pointers!
	typedef fixed_index_list<7> edge_list;
	object_pool<edge_list::node> edge_list_pool;
	VectorType<edge_list> m_vVertexEdges;
    
	EdgeID add_edge(VertexID vA, VertexID vB, TriangleID tA, TriangleID tB = InvalidID);
    void remove_edge(EdgeID eID);
	Edge * find_edge(VertexID vA, VertexID vB, EdgeID & eID);
	const Edge * find_edge(VertexID vA, VertexID vB, EdgeID & eID) const;
    EdgeID next_edge(TriangleID tID, VertexID vCenter, VertexID vOther) const;
    
    bool isBoundaryT(TriangleID tID) const;
    bool hasNeighbourT(TriangleID tCheck, TriangleID tNbr) const;
    
public:

	//
	// iterators
	//   The functions vertices() / triangles() / edges() are provided so you can do:
	//      for ( eid : edges() ) { ... }
	//   and other related begin() / end() idioms
	//
	typedef typename rcvector<Vertex, VectorType>::index_iterator vertex_iterator;
	typedef typename rcvector<Vertex, VectorType>::index_wrapper vertex_iterator_wrap;
	vertex_iterator begin_vertices();
	vertex_iterator end_vertices();
	vertex_iterator_wrap vertices();

	typedef typename rcvector<Triangle, VectorType>::index_iterator triangle_iterator;
	typedef typename rcvector<Triangle, VectorType>::index_wrapper triangle_iterator_wrap;
	triangle_iterator begin_triangles();
	triangle_iterator end_triangles();
	triangle_iterator_wrap triangles();

	typedef typename rcvector<Edge, VectorType>::index_iterator edge_iterator;
	typedef typename rcvector<Edge, VectorType>::index_wrapper edge_iterator_wrap;
	edge_iterator begin_edges();
	edge_iterator end_edges();
	edge_iterator_wrap edges();

	//
	// IDynamicMesh interface
	//
    virtual bool CheckValidity(bool bAssert = true) override;
};
	
	
typedef OldDMesh3<float> OldDMesh3f;
typedef OldDMesh3<double> OldDMesh3d;



template<typename Real, template<typename> class VectorType>
unsigned int OldDMesh3<Real, VectorType>::GetVertexCount( ) const
{
	return m_vVertices.count();
}
template<typename Real, template<typename> class VectorType>
unsigned int OldDMesh3<Real, VectorType>::GetTriangleCount( ) const
{
	return m_vTriangles.count();
}
template<typename Real, template<typename> class VectorType>
unsigned int OldDMesh3<Real, VectorType>::GetEdgeCount( ) const
{
	return m_vEdges.count();
}


template<typename Real, template<typename> class VectorType>
VertexID OldDMesh3<Real, VectorType>::GetMaxVertexID( ) const
{
	return m_vVertices.max_index();
}
template<typename Real, template<typename> class VectorType>
TriangleID OldDMesh3<Real, VectorType>::GetMaxTriangleID( ) const
{
	return m_vTriangles.max_index();
}
template<typename Real, template<typename> class VectorType>
EdgeID OldDMesh3<Real, VectorType>::GetMaxEdgeID( ) const
{
	return m_vEdges.max_index();
}


template<typename Real, template<typename> class VectorType>
bool OldDMesh3<Real, VectorType>::IsVertex( VertexID vID ) const
{
	return m_vVertices.isValid(vID);
}
template<typename Real, template<typename> class VectorType>
bool OldDMesh3<Real, VectorType>::IsTriangle( TriangleID tID ) const
{
	return m_vTriangles.isValid(tID);
}
template<typename Real, template<typename> class VectorType>
bool OldDMesh3<Real, VectorType>::IsEdge( EdgeID eID ) const
{
	return m_vEdges.isValid(eID);
}


template<typename Real, template<typename> class VectorType>
const Vector3<Real> & OldDMesh3<Real, VectorType>::GetVertex( VertexID vID ) const
{
	return m_vVertices.isValid(vID) ? 
		m_vVertices[vID].v : InvalidVertex;
}
template<typename Real, template<typename> class VectorType>
void OldDMesh3<Real, VectorType>::SetVertex( VertexID vID, const Vector3<Real> & v )
{
	if ( m_vVertices.isValid(vID) )
		m_vVertices[vID].v = v;
}



template<typename Real, template<typename> class VectorType>
const Vector3i & OldDMesh3<Real, VectorType>::GetTriangle( TriangleID tID ) const
{
	return m_vTriangles.isValid(tID) ? 
		m_vTriangles[tID].tv : InvalidTriangle;
}

template<typename Real, template<typename> class VectorType>
const Vector2i & OldDMesh3<Real, VectorType>::GetEdgeV( EdgeID eID ) const
{
	return m_vEdges.isValid(eID) ? 
		m_vEdges[eID].v : InvalidEdge;
}

template<typename Real, template<typename> class VectorType>
const Vector2i & OldDMesh3<Real, VectorType>::GetEdgeT( EdgeID eID ) const
{
	return m_vEdges.isValid(eID) ? 
		m_vEdges[eID].t : InvalidEdge;
}


template<typename Real, template<typename> class VectorType>
EdgeID OldDMesh3<Real, VectorType>::FindEdge( VertexID vA, VertexID vB ) const
{
	EdgeID eID;
	const Edge * e = find_edge(vA, vB, eID );
	return (e == nullptr) ? InvalidID : eID;
}

template<typename Real, template<typename> class VectorType>
MeshResult OldDMesh3<Real, VectorType>::GetVtxEdges(VertexID vID,  std::vector<int> & v ) const
{
	if (! IsVertex(vID) )
		return MeshResult::Failed_NotAVertex;
	m_vVertexEdges[vID].get(v);
	return MeshResult::Ok;
}

template<typename Real, template<typename> class VectorType>
Vector3i OldDMesh3<Real, VectorType>::GetTriEdges(TriangleID tID) const
{
    if (! IsTriangle(tID) )
        return Vector3i(InvalidID, InvalidID, InvalidID);
    return m_vTriangles[tID].te;
}


template<typename Real, template<typename> class VectorType>
typename OldDMesh3<Real, VectorType>::vertex_iterator OldDMesh3<Real, VectorType>::begin_vertices()
{
	return m_vVertices.begin_indices();
}
template<typename Real, template<typename> class VectorType>
typename OldDMesh3<Real, VectorType>::vertex_iterator OldDMesh3<Real, VectorType>::end_vertices()
{
	return m_vVertices.end_indices();
}
template<typename Real, template<typename> class VectorType>
typename OldDMesh3<Real, VectorType>::vertex_iterator_wrap OldDMesh3<Real, VectorType>::vertices()
{
	return vertex_iterator_wrap(&m_vVertices);
}


template<typename Real, template<typename> class VectorType>
typename OldDMesh3<Real, VectorType>::triangle_iterator OldDMesh3<Real, VectorType>::begin_triangles()
{
	return m_vTriangles.begin_indices();
}
template<typename Real, template<typename> class VectorType>
typename OldDMesh3<Real, VectorType>::triangle_iterator OldDMesh3<Real, VectorType>::end_triangles()
{
	return m_vTriangles.end_indices();
}
template<typename Real, template<typename> class VectorType>
typename OldDMesh3<Real, VectorType>::triangle_iterator_wrap OldDMesh3<Real, VectorType>::triangles()
{
	return triangle_iterator_wrap(&m_vTriangles);
}


template<typename Real, template<typename> class VectorType>
typename OldDMesh3<Real, VectorType>::edge_iterator OldDMesh3<Real, VectorType>::begin_edges()
{
	return m_vEdges.begin_indices();
}
template<typename Real, template<typename> class VectorType>
typename OldDMesh3<Real, VectorType>::edge_iterator OldDMesh3<Real, VectorType>::end_edges()
{
	return m_vEdges.end_indices();
}
template<typename Real, template<typename> class VectorType>
typename OldDMesh3<Real, VectorType>::edge_iterator_wrap OldDMesh3<Real, VectorType>::edges()
{
	return edge_iterator_wrap(&m_vEdges);
}


	
}