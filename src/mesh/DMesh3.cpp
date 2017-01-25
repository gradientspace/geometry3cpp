#include <geometry3PCH.h>
#include <DMesh3.h>

#include <g3Colors.h>
#include <g3Debug.h>
#include <index_util.h>

using namespace g3;

#define DMESH_TEMPLATE template<typename Real, template<typename> class VectorType>
#define DMESH_TYPE DMesh3<Real,VectorType>


DMESH_TEMPLATE Vector3<Real> DMESH_TYPE::InvalidVertex = Vector3<Real>( std::numeric_limits<Real>::max(), 0, 0  );
DMESH_TEMPLATE Vector3i DMESH_TYPE::InvalidTriangle = Vector3i(InvalidID, InvalidID, InvalidID);
DMESH_TEMPLATE Vector2i DMESH_TYPE::InvalidEdge = Vector2i(InvalidID, InvalidID);


DMESH_TEMPLATE
DMESH_TYPE::~DMesh3()
{
}

DMESH_TEMPLATE
DMESH_TYPE::DMesh3(MeshConfigFlags flags)
{
	m_flags = flags;
}

DMESH_TEMPLATE
DMESH_TYPE::DMesh3( const DMesh3T & copy )
{
	*this = copy;
}


DMESH_TEMPLATE
DMESH_TYPE::DMesh3( const DMesh3T && moved )
{
	*this = std::move(moved);
}


DMESH_TEMPLATE
const DMESH_TYPE & DMESH_TYPE::operator=( const DMesh3T & copy )
{
	m_flags = copy.m_flags;

	m_vVertices = copy.m_vVertices;
	m_vNormals = copy.m_vNormals;
	m_vColors = copy.m_vColors;
	m_vUVs = copy.m_vUVs;
	m_vTriColors = copy.m_vTriColors;

	m_vTriangles = copy.m_vTriangles;
	m_vEdges = copy.m_vEdges;

	// manually copy the VertexEdges because there is no overall data structure to operator=
	edge_list_pool.free_pool();
	m_vVertexEdges.clear();
	m_vVertexEdges.resize( copy.m_vVertexEdges.size() );
	auto end = m_vVertices.end();
	for ( auto cur = m_vVertices.begin(); cur != end; cur++) {
		unsigned int i = cur.index();
		m_vVertexEdges[i].copy( copy.m_vVertexEdges[i], &edge_list_pool );
	}

	return *this;
}

DMESH_TEMPLATE 
const DMESH_TYPE & DMESH_TYPE::operator=( const DMesh3T && moved )
{
	m_flags = moved.m_flags;

	m_vVertices = std::move(moved.m_vVertices);
	m_vNormals = std::move(moved.m_vNormals);
	m_vColors = std::move(moved.m_vColors);
	m_vUVs = std::move(moved.m_vUVs);
	m_vTriColors = std::move(moved.m_vTriColors);

	m_vTriangles = std::move(moved.m_vTriangles);
	m_vEdges = std::move(moved.m_vEdges);

	edge_list_pool = std::move(moved.edge_list_pool);
	m_vVertexEdges = std::move(moved.m_vVertexEdges);

	return *this;
}





DMESH_TEMPLATE 
void DMESH_TYPE::SetEnableVertexNormals( bool bEnable )
{
	if (bEnable && (m_flags & VertexNormals) == 0 ) {
		m_flags |= VertexNormals;
		m_vNormals.resize( m_vVertices.max_index() );
	} else if (bEnable == false && (m_flags & VertexNormals) != 0) {
		m_flags ^= VertexNormals;
		m_vNormals.clear(true);
	}
}
DMESH_TEMPLATE 
void DMESH_TYPE::SetEnableVertexColors( bool bEnable )
{
	if (bEnable && (m_flags & VertexColors) == 0 ) {
		m_flags |= VertexColors;
		m_vColors.resize( m_vVertices.max_index() );
	} else if (bEnable == false && (m_flags & VertexColors) != 0) {
		m_flags ^= VertexColors;
		m_vColors.clear(true);
	}	
}
DMESH_TEMPLATE 
void DMESH_TYPE::SetEnableVertexUVs( bool bEnable )
{
	if (bEnable && (m_flags & VertexUVs) == 0 ) {
		m_flags |= VertexUVs;
		m_vUVs.resize( m_vVertices.max_index() );
	} else if (bEnable == false && (m_flags & VertexUVs) != 0) {
		m_flags ^= VertexUVs;
		m_vUVs.clear(true);
	}	
}
DMESH_TEMPLATE 
void DMESH_TYPE::SetEnableTriangleColors( bool bEnable )
{
	if (bEnable && (m_flags & TriangleColors) == 0 ) {
		m_flags |= TriangleColors;
		m_vTriColors.resize( m_vTriangles.max_index() );
	} else if (bEnable == false && (m_flags & TriangleColors) != 0) {
		m_flags ^= TriangleColors;
		m_vTriColors.clear(true);
	}	
}



DMESH_TEMPLATE 
VertexID DMESH_TYPE::AppendVertex( const Vector3<Real> & v )
{
	VertexID vID = m_vVertices.insert( { v, 0, 0 } );
	if ( GetEnableVertexNormals() )
		m_vNormals.insert( vID, Wml::Vector3f::UNIT_X );
	if ( GetEnableVertexColors() )
		m_vColors.insert( vID, Colors::White );
	if ( GetEnableVertexUVs() )
		m_vUVs.insert( vID, Wml::Vector2f::ZERO );

	if ( m_vVertexEdges.size() < vID )
		m_vVertexEdges.resize(vID);

	return vID;
}



DMESH_TEMPLATE 
TriangleID DMESH_TYPE::AppendTriangle( const Vector3i & tv, GroupID gID )
{
	if (IsVertex( tv[0] ) == false || IsVertex( tv[1] ) == false || IsVertex( tv[2] ) == false) {
		gDevAssert( false );
		return InvalidID;
	}
	if (tv[0] == tv[1] || tv[0] == tv[2] || tv[1] == tv[2]) {
		gDevAssert( false );
		return InvalidID;
	}

	// look up edges. if any already have two triangles, this would create non-manifold geometry
	// and so we do not allow it
	Edge * e[3];
	EdgeID eID[3];
	for (int j = 0; j < 3; ++j) {
		e[j] = find_edge( tv[j], tv[(j+1)%3], eID[j] );
		if (e[j] != nullptr && e[j]->isBoundary() == false)
			return NonManifoldID;
	}

	// now safe to insert triangle
	TriangleID tID = m_vTriangles.insert( { tv, Vector3i(InvalidID, InvalidID, InvalidID), {gID,0,0} } );
    Triangle & newt = m_vTriangles[tID];

	// increment ref counts and update/create edges
	for (int j = 0; j < 3; ++j) {
		m_vVertices.increment(tv[j]);

        if ( e[j] != nullptr ) {
			e[j]->t[1] = tID;
            newt.te[j] = eID[j];
        } else {
			newt.te[j] = add_edge(tv[j], tv[(j+1)%3], tID);
        }
	}

    ITimeStamped::updateTimeStamp();
    
	return tID;
}




DMESH_TEMPLATE
MeshResult DMESH_TYPE::RemoveTriangle(TriangleID tID, bool bDeleteUnrefVertices)
{
    if ( ! IsTriangle(tID) )
        return MeshResult::Failed_NotATriangle;
    
    // decrement tri refcount, which should free it
    Vector3i tv = m_vTriangles[tID].tv;
	Vector3i nbrs = m_vTriangles[tID].te;
    m_vTriangles.decrement(tID);
    gDevAssert( m_vTriangles.isValid(tID) == false );
    
    // look up edges
    Edge * e[3];
    for (int j = 0; j < 3; ++j)
        e[j] = &m_vEdges[nbrs[j]];
    
    // update edges
    for ( int j = 0; j < 3; ++j ) {
        if ( e[j]->t[1] == tID ) {
            // edge has two tris and we are second
            e[j]->t[1] = InvalidID;
        } else if ( e[j]->t[0] == tID && e[j]->t[1] != InvalidID ) {
            // two tris, we are first
            e[j]->t[0] = e[j]->t[1];
			e[j]->t[1] = InvalidID;
		} else {
            // edge has one tri and we are it - have to remove
            e[j]->t[0] = InvalidID;
        }
    }
    
    // remove edges that no longer have triangles
    for ( int j = 0; j < 3; ++j ) {
        if ( e[j]->t[0] == InvalidID )
            remove_edge(nbrs[j]);
    }
    
    // decrement vertex ref counts
    for ( int j = 0; j < 3; ++j ) {
        m_vVertices.decrement(tv[j]);
        
        if ( bDeleteUnrefVertices && m_vVertexEdges[tv[j]].empty() ) {
            m_vVertices.decrement(tv[j]);
            gDevAssert( IsVertex(tv[j]) == false );
        }
    }

    ITimeStamped::updateTimeStamp();
    
    return MeshResult::Ok;
}


DMESH_TEMPLATE
MeshResult DMESH_TYPE::ReverseTriOrientation(TriangleID tID)
{
    if ( ! IsTriangle(tID) )
        return MeshResult::Failed_NotATriangle;
    Triangle & t = m_vTriangles[tID];
    t.tv = Vector3i(t.tv[1], t.tv[0], t.tv[2]);
    t.te = Vector3i(t.te[0], t.te[2], t.te[1]);
    return MeshResult::Ok;
}



DMESH_TEMPLATE
MeshResult DMESH_TYPE::FlipEdge( VertexID vA, VertexID vB, EdgeFlipInfo & flip )
{
	EdgeID eab;
	Edge * eAB = find_edge(vA,vB, eab);
	if ( eAB == nullptr )
		return MeshResult::Failed_NotAnEdge;
	return FlipEdge(eab, flip);
}
DMESH_TEMPLATE
MeshResult DMESH_TYPE::FlipEdge( EdgeID eab, EdgeFlipInfo & flip )
{
	if (! IsEdge(eab) )
		return MeshResult::Failed_NotAnEdge;
	Edge * eAB = & m_vEdges[eab];
	if ( eAB->isBoundary() )
		return MeshResult::Failed_IsBoundaryEdge;

	// find oriented edge [a,b], tris t0,t1, and other verts c in t0, d in t1
	VertexID a = eAB->v[0], b = eAB->v[1];
	TriangleID t0 = eAB->t[0], t1 = eAB->t[1];
	Triangle & T0 = m_vTriangles[t0];
	Triangle & T1 = m_vTriangles[t1];
	VertexID c = orient_tri_edge_and_find_other_vtx( a, b, T0.tv );
	VertexID d = find_tri_other_vtx(a, b, T1.tv);
	if ( c == InvalidID || d == InvalidID ) {
		gDevAssert(false);
		return MeshResult::Failed_BrokenTopology;
	}

	EdgeID flipped;
	if ( find_edge(c,d, flipped) != nullptr )
		return MeshResult::Failed_FlippedEdgeExists;

	// find edges bc, ca, ad, db
    EdgeID ebc = T0.findNeighbourEdge(b,c);
    EdgeID eca = T0.findNeighbourEdge(c,a);
    Edge * eCA = &m_vEdges[eca];
    EdgeID ead = T1.findNeighbourEdge(a,d);
    EdgeID edb = T1.findNeighbourEdge(d,b);
    Edge * eDB = &m_vEdges[edb];

	// update triangles
	T0.tv = {c,d,b};
	T1.tv = {d,c,a};

	// update edge AB, which becomes flipped edge CD
	eAB->v = { std::min(c,d), std::max(c,d) };
	eAB->t = { t0, t1 };
	EdgeID ecd = eab;

	// update the two other edges whose triangle nbrs have changed
	if ( eCA->replaceTriangle(t0, t1) == -1 ) gDevAssert(false);
	if ( eDB->replaceTriangle(t1, t0) == -1 ) gDevAssert(false);

    // update triangle nbr lists (these are edges)
    T0.te = {ecd, edb, ebc};
    T1.te = {ecd, eca, ead};
    
	// remove old eab from verts a and b, and decrement ref counts
	if ( m_vVertexEdges[a].remove(eab, &edge_list_pool) == false ) gDevAssert(false);
	if ( m_vVertexEdges[b].remove(eab, &edge_list_pool) == false ) gDevAssert(false);
	m_vVertices.decrement(a);
	m_vVertices.decrement(b);
	gDevAssert( IsVertex(a) && IsVertex(b) );

	// add new edge ecd to verts c and d, and increment ref counts
	m_vVertexEdges[c].add(ecd, &edge_list_pool);
	m_vVertexEdges[d].add(ecd, &edge_list_pool);
	m_vVertices.increment(c);
	m_vVertices.increment(d);

	// success! collect up results
	flip.eID = eab;
	flip.v = {a,b};
	flip.ov = {c,d};
	flip.t = {t0,t1};
    
    ITimeStamped::updateTimeStamp();

	return MeshResult::Ok;
}



DMESH_TEMPLATE
MeshResult DMESH_TYPE::SplitEdge( VertexID vA, VertexID vB, EdgeSplitInfo & split )
{
	EdgeID eab;
	Edge * eAB = find_edge(vA,vB, eab);
	if ( eAB == nullptr )
		return MeshResult::Failed_NotAnEdge;
	return SplitEdge(eab, split);
}
template<typename Real, template<typename> class VectorType>
MeshResult DMesh3<Real, VectorType>::SplitEdge( EdgeID eab, EdgeSplitInfo & split )
{
	if (! IsEdge(eab) )
		return MeshResult::Failed_NotAnEdge;
	Edge * eAB = & m_vEdges[eab];

	// look up primary edge & triangle
	TriangleID t0 = eAB->t[0];
	Triangle & T0 = m_vTriangles[t0];
	VertexID a = eAB->v[0], b = eAB->v[1];
	VertexID c = orient_tri_edge_and_find_other_vtx(a, b, T0.tv);

	// create new vertex
	Vector3<Real> vNew = (Real)0.5 * ( m_vVertices[a].v + m_vVertices[b].v );
	VertexID f = AppendVertex( vNew );

	// quite a bit of code is duplicated between boundary and non-boundary case, but it
	//  is too hard to follow later if we factor it out...
	if ( eAB->isBoundary() ) {

		// look up edge bc, which needs to be modified
		EdgeID ebc = T0.te[ find_edge_index_in_tri(b, c, T0.tv) ];
		Edge & eBC = m_vEdges[ebc];

		// rewrite existing triangle
		T0.replaceV(b, f);

		// add new second triangle
		TriangleID t2 = m_vTriangles.insert( { Vector3i(f, b, c), Vector3i(InvalidID, InvalidID, InvalidID), {T0.data.gid,0,0} } );
		Triangle & T2 = m_vTriangles[t2];

		// rewrite edge bc, create edge af
		eBC.replaceTriangle(t0, t2);
		VertexID eaf = eab; Edge * eAF = eAB;
		eAF->replaceVertex(b, f);
		m_vVertexEdges[b].remove(eab, &edge_list_pool);
		m_vVertexEdges[f].add(eaf, &edge_list_pool);

		// create new edges fb and fc 
		EdgeID efb = add_edge(f, b, t2);
		EdgeID efc = add_edge(f, c, t0, t2);

		// update triangle edge-nbrs
		T0.replaceE(ebc, efc);
		T2.te = Vector3i(efb, ebc, efc);

		// update vertex refcounts
		m_vVertices.increment(c);
		m_vVertices.increment(f, 2);

		split.bIsBoundary = true;
		split.vNew = f;

		ITimeStamped::updateTimeStamp();
		return MeshResult::Ok;

	} else {		// interior triangle branch

		// look up other triangle
		TriangleID t1 = eAB->t[1];
		Triangle & T1 = m_vTriangles[t1];
		VertexID d = find_tri_other_vtx( a, b, T1.tv );

		// look up edges that we are going to need to update
		// [TODO OPT] could use ordering to reduce # of compares here
		EdgeID ebc = T0.te[find_edge_index_in_tri( b, c, T0.tv )];
		Edge & eBC = m_vEdges[ebc];
		EdgeID edb = T1.te[find_edge_index_in_tri( d, b, T1.tv )];
		Edge & eDB = m_vEdges[edb];

		// rewrite existing triangles
		T0.replaceV( b, f );
		T1.replaceV( b, f );

		// add two new triangles to close holes we just created
		TriangleID t2 = m_vTriangles.insert( { Vector3i( f, b, c ), Vector3i( InvalidID, InvalidID, InvalidID ), {T0.data.gid,0,0} } );
		Triangle & T2 = m_vTriangles[t2];
		TriangleID t3 = m_vTriangles.insert( { Vector3i( f, d, b ), Vector3i( InvalidID, InvalidID, InvalidID ), {T1.data.gid,0,0} } );
		Triangle & T3 = m_vTriangles[t3];

		// update the edges we found above, to point to new triangles
		eBC.replaceTriangle( t0, t2 );
		eDB.replaceTriangle( t1, t3 );

		// edge eab became eaf
		VertexID eaf = eab; Edge * eAF = eAB;
		eAF->replaceVertex( b, f );

		// update a/b/f vertex-edges
		m_vVertexEdges[b].remove( eab, &edge_list_pool );
		m_vVertexEdges[f].add( eaf, &edge_list_pool );

		// create new edges connected to f  (also updates vertex-edges)
		EdgeID efb = add_edge( f, b, t2, t3 );
		EdgeID efc = add_edge( f, c, t0, t2 );
		EdgeID edf = add_edge( d, f, t1, t3 );

		// update triangle edge-nbrs
		T0.replaceE( ebc, efc );
		T1.replaceE( edb, edf );
		T2.te = Vector3i( efb, ebc, efc );
		T3.te = Vector3i( edf, edb, efb );

		// update vertex refcounts
		m_vVertices.increment( c );
		m_vVertices.increment( d );
		m_vVertices.increment( f, 4 );

		split.bIsBoundary = false;
		split.vNew = f;

		ITimeStamped::updateTimeStamp();
		return MeshResult::Ok;
	}
}

#ifdef _DEBUG
#define DMESH3_CHECK_EQ(x, y) gDevAssert((x) == (y))
#define DMESH3_CHECK_NEQ(x, y) gDevAssert((x) != (y))
#else
#define DMESH3_CHECK_EQ(x, y) x
#define DMESH3_CHECK_NEQ(x, y) x
#endif

template<typename Real, template<typename> class VectorType>
MeshResult DMesh3<Real, VectorType>::CollapseEdge( VertexID vKeep, VertexID vRemove, EdgeCollapseInfo & collapse )
{
	if ( IsVertex(vKeep) == false || IsVertex(vRemove) == false )
		return MeshResult::Failed_NotAnEdge;

	VertexID b = vKeep;		// renaming for sanity. We remove a and keep b
	VertexID a = vRemove;
	edge_list & edges_b = m_vVertexEdges[b];

	EdgeID eab;
	Edge * eAB = find_edge( a, b, eab );
	if (eAB == nullptr)
		return MeshResult::Failed_NotAnEdge;

	TriangleID t0 = eAB->t[0];
	Triangle & T0 = m_vTriangles[t0];
	VertexID c = find_tri_other_vtx(a, b, T0.tv);

	// look up opposing triangle/vtx if we are not in boundary case
	bool bIsBoundary = false;
	VertexID d = InvalidID;
	TriangleID t1 = eAB->t[1];
	if (t1 != InvalidID) {
		Triangle & T1 = m_vTriangles[t1];
		d = find_tri_other_vtx( a, b, T1.tv );
		if (c == d)
			return MeshResult::Failed_FoundDuplicateTriangle;
	} else {
		bIsBoundary = true;
	}

	// We cannot collapse if edge lists of a and b share vertices other
	//  than c and d  (because then we will make a triangle [x b b].
	//  Unfortunately I cannot see a way to do this more efficiently than brute-force search
	//  [TODO] if we had tri iterator for a, couldn't we check each tri for b  (skipping t0 and t1) ?
	edge_list & edges_a = m_vVertexEdges[a];
	int edges_a_count = 0;
	for (auto eid_a : edges_a) {
		VertexID vax = m_vEdges[eid_a].otherV(a);
		edges_a_count++;
		if ( vax == b || vax == c || vax == d )
			continue;
		for (auto eid_b : edges_b) {
			if ( m_vEdges[eid_b].otherV(b) == vax )
				return MeshResult::Failed_InvalidNeighbourhood;
		}
	}

	// We cannot collapse if we have a tetrahedron. In this case a has 3 nbr edges,
	//  and edge cd exists. But that is not conclusive, we also have to check that
	//  cd is an internal edge, and that each of its tris contain a or b
	if (edges_a_count == 3 && bIsBoundary == false) {
		EdgeID edc;
		Edge * eDC = find_edge( d, c, edc );
		if (eDC != nullptr && eDC->t[1] != InvalidID ) {
			Triangle & X0 = m_vTriangles[eDC->t[0]];
			Triangle & X1 = m_vTriangles[eDC->t[1]];
			if ( (X0.hasV(a) && X1.hasV(b)) || (X0.hasV(b) && X1.hasV(a)) )
				return MeshResult::Failed_CollapseTetrahedron;
		}

	} else if (edges_a_count == 2 && bIsBoundary == true) {
		// cannot collapse edge if we are down to a single triangle
		if ( edges_b.size() == 2 && m_vVertexEdges[c].size() == 2 )
			return MeshResult::Failed_CollapseTriangle;
	}


	// 1) remove edge ab from vtx b
	// 2) find edges ad and ac, and tris tad, tac across those edges  (will use later)
	// 3) for other edges, replace a with b, and add that edge to b
	// 4) replace a with b in all triangles connected to a
	EdgeID ead = InvalidID, eac = InvalidID;
	TriangleID tad = InvalidID, tac = InvalidID;
	for (auto eid : edges_a) {
		Edge & e = m_vEdges[eid];
		VertexID o = e.otherV(a);
		if (o == b) {
			DMESH3_CHECK_EQ( edges_b.remove(eid, &edge_list_pool), true );
		} else if (o == c) {
			eac = eid;
			DMESH3_CHECK_EQ( m_vVertexEdges[c].remove(eid, &edge_list_pool), true );
			tac = e.otherT(t0);
		} else if (o == d) {
			ead = eid;
			DMESH3_CHECK_EQ( m_vVertexEdges[d].remove(eid, &edge_list_pool), true );
			DMESH3_CHECK_NEQ( tad = e.otherT(t1), InvalidID );
		} else {
			DMESH3_CHECK_NEQ( e.replaceVertex(a, b), -1 );
			edges_b.add(eid, &edge_list_pool);
		}

		// [TODO] perhaps we can already have unique tri list because of the manifold-nbrhood check we need to do...
		for (int j = 0; j < 2; ++j) {
			if (e.t[j] != InvalidID && e.t[j] != t0 && e.t[j] != t1) {
				Triangle & T = m_vTriangles[e.t[j]];
				if (T.hasV( a )) {
					DMESH3_CHECK_NEQ( T.replaceV( a, b ), -1 );
					m_vVertices.increment(b);
					m_vVertices.decrement(a);
				}
			}
		}
	}

	if (bIsBoundary == false) {

		// remove all edges from vtx a, then remove vtx a
		edges_a.clear( &edge_list_pool );
		gDevAssert( m_vVertices.refCount( a ) == 3 );		// in t0,t1, and initial ref
		m_vVertices.decrement( a, 3 );
		gDevAssert( m_vVertices.isValid( a ) == false );

		// remove triangles T0 and T1, and update b/c/d refcounts
		m_vTriangles.decrement( t0 );
		m_vTriangles.decrement( t1 );
		m_vVertices.decrement( c );
		m_vVertices.decrement( d );
		m_vVertices.decrement( b, 2 );
		gDevAssert( m_vTriangles.isValid( t0 ) == false );
		gDevAssert( m_vTriangles.isValid( t1 ) == false );

		// remove edges ead, eab, eac
		m_vEdges.decrement( ead );
		m_vEdges.decrement( eab );
		m_vEdges.decrement( eac );
		gDevAssert( m_vEdges.isValid( ead ) == false );
		gDevAssert( m_vEdges.isValid( eab ) == false );
		gDevAssert( m_vEdges.isValid( eac ) == false );

		// replace t0 and t1 in edges ebd and ebc that we kept
		EdgeID ebd, ebc;
		Edge * eBD = find_edge( b, d, ebd );
		Edge * eBC = find_edge( b, c, ebc );
		DMESH3_CHECK_NEQ( eBD->replaceTriangle( t1, tad ), -1 );
		DMESH3_CHECK_NEQ( eBC->replaceTriangle( t0, tac ), -1 );

		// update tri-edge-nbrs in tad and tac
		if (tad != InvalidID) {
			DMESH3_CHECK_NEQ( m_vTriangles[tad].replaceE( ead, ebd ), -1 );
		}
		if (tac != InvalidID) {
			DMESH3_CHECK_NEQ( m_vTriangles[tac].replaceE( eac, ebc ), -1 );
		}

	} else {
		//  this is basically same code as above, just not referencing t0/d

		// remove all edges from vtx a, then remove vtx a
		edges_a.clear( &edge_list_pool );
		gDevAssert( m_vVertices.refCount( a ) == 2 );		// in t0 and initial ref
		m_vVertices.decrement( a, 2 );
		gDevAssert( m_vVertices.isValid( a ) == false );

		// remove triangle T0 and update b/c refcounts
		m_vTriangles.decrement( t0 );
		m_vVertices.decrement( c );
		m_vVertices.decrement( b );
		gDevAssert( m_vTriangles.isValid( t0 ) == false );

		// remove edges eab and eac
		m_vEdges.decrement( eab );
		m_vEdges.decrement( eac );
		gDevAssert( m_vEdges.isValid( eab ) == false );
		gDevAssert( m_vEdges.isValid( eac ) == false );

		// replace t0 in edge ebc that we kept
		EdgeID ebc;
		Edge * eBC = find_edge( b, c, ebc );
		DMESH3_CHECK_NEQ( eBC->replaceTriangle( t0, tac ), -1 );

		// update tri-edge-nbrs in tac
		if (tac != InvalidID) {
			DMESH3_CHECK_NEQ( m_vTriangles[tac].replaceE( eac, ebc ), -1 );
		}
	}

	ITimeStamped::updateTimeStamp();
	return MeshResult::Ok;
}


DMESH_TEMPLATE
bool DMESH_TYPE::isBoundaryT(TriangleID tID) const
{
    const Triangle & t = m_vTriangles[tID];
    if ( m_vEdges[t.te[0]].isBoundary() ) return true;
    if ( m_vEdges[t.te[1]].isBoundary() ) return true;
    if ( m_vEdges[t.te[2]].isBoundary() ) return true;
    return false;
}

DMESH_TEMPLATE
bool DMESH_TYPE::hasNeighbourT(TriangleID tCheck, TriangleID tNbr) const
{
    const Triangle & t = m_vTriangles[tCheck];
    if ( m_vEdges[t.te[0]].hasT(tNbr) ) return true;
    if ( m_vEdges[t.te[1]].hasT(tNbr) ) return true;
    if ( m_vEdges[t.te[2]].hasT(tNbr) ) return true;
    return false;
}

DMESH_TEMPLATE
bool DMESH_TYPE::Triangle::hasSequentialV(VertexID vA, VertexID vB) const
{
    if ( tv[0] == vA && tv[1] == vB ) return true;
    if ( tv[1] == vA && tv[2] == vB ) return true;
    if ( tv[2] == vA && tv[0] == vB ) return true;
    return false;
}

DMESH_TEMPLATE
int DMESH_TYPE::Triangle::replaceV( VertexID vOld, VertexID vNew )
{
	if ( tv[0] == vOld ) { tv[0] = vNew; return 0; }
	if ( tv[1] == vOld ) { tv[1] = vNew; return 1; }
	if ( tv[2] == vOld ) { tv[2] = vNew; return 1; }
	return -1;
}

DMESH_TEMPLATE
int DMESH_TYPE::Triangle::findNeighbourIndex(VertexID vA, VertexID vB) const
{
    if ( same_pair_unordered(tv[0], tv[1], vA, vB) ) return 0;
    if ( same_pair_unordered(tv[1], tv[2], vA, vB) ) return 1;
    if ( same_pair_unordered(tv[2], tv[0], vA, vB) ) return 2;
    return -1;
}

DMESH_TEMPLATE
EdgeID DMESH_TYPE::Triangle::findNeighbourEdge(VertexID vA, VertexID vB) const
{
    if ( same_pair_unordered(tv[0], tv[1], vA, vB) ) return te[0];
    if ( same_pair_unordered(tv[1], tv[2], vA, vB) ) return te[1];
    if ( same_pair_unordered(tv[2], tv[0], vA, vB) ) return te[2];
    return InvalidID;
}

DMESH_TEMPLATE
int DMESH_TYPE::Triangle::replaceE( EdgeID eOld, EdgeID eNew )
{
	if ( te[0] == eOld ) { te[0] = eNew; return 0; }
	if ( te[1] == eOld ) { te[1] = eNew; return 1; }
	if ( te[2] == eOld ) { te[2] = eNew; return 1; }
	return -1;
}

DMESH_TEMPLATE
int DMESH_TYPE::Edge::replaceTriangle( TriangleID tOld, TriangleID tNew )
{
    if ( t[0] == tOld ) {
		if (tNew == InvalidID) {
			t[0] = t[1]; 
			t[1] = InvalidID;
		} else
			t[0] = tNew; 
		return 0;
    } else if ( t[1] == tOld ) {
        t[1] = tNew; 
		return 1;
    } else
        return -1;
}

DMESH_TEMPLATE
int DMESH_TYPE::Edge::replaceVertex(VertexID vOld, VertexID vNew)
{
	if ( v[0] == vOld ) {
		v = { std::min(v[1],vNew), std::max(v[1],vNew) }; 
		return 0;
	} else if ( v[1] == vOld ) {
		v = { std::min(v[0],vNew), std::max(v[0],vNew) }; 
		return 1;
	} else
		return -1;
}

DMESH_TEMPLATE
EdgeID DMESH_TYPE::add_edge( VertexID vA, VertexID vB, TriangleID tA, TriangleID tB )
{
	Edge e = { {vA,vB}, {tA, tB}, 0, 0 };
	if ( vB < vA )
		std::swap(e.v[0], e.v[1]);
	EdgeID eID = m_vEdges.insert(e);

	m_vVertexEdges[vA].add(eID, &edge_list_pool);
	m_vVertexEdges[vB].add(eID, &edge_list_pool);
	return eID;
}

DMESH_TEMPLATE
void DMESH_TYPE::remove_edge(EdgeID eID)
{
    Edge e = m_vEdges[eID];
    gDevAssert(e.t[0] == InvalidID && e.t[1] == InvalidID);
    
    m_vVertexEdges[e.v[0]].remove(eID, &edge_list_pool);
    m_vVertexEdges[e.v[1]].remove(eID, &edge_list_pool);
    
    m_vEdges.decrement(eID);
    gDevAssert( m_vEdges.isValid(eID) == false );
}

DMESH_TEMPLATE 
typename DMESH_TYPE::Edge * DMESH_TYPE::find_edge(VertexID vA, VertexID vB, EdgeID & eID)
{
	VertexID vO = std::max(vA,vB);
	const edge_list & e0 = m_vVertexEdges[ std::min(vA,vB) ];
	eID = e0.find( [this,vO]( EdgeID e )-> bool { 
		return m_vEdges[e].hasV(vO);
	}, InvalidID);
	return ( eID != InvalidID ) ? &m_vEdges[eID] : nullptr;
}
DMESH_TEMPLATE 
const typename DMESH_TYPE::Edge * DMESH_TYPE::find_edge( VertexID vA, VertexID vB, EdgeID & eID ) const
{
	return const_cast<DMESH_TYPE *>(this)->find_edge(vA,vB,eID);
}



DMESH_TEMPLATE 
Vector2i DMESH_TYPE::GetEdgeOpposingV( EdgeID eID ) const
{
	const Edge * eAB = & m_vEdges[eID];
	VertexID a = eAB->v[0], b = eAB->v[1];
	TriangleID t0 = eAB->t[0];
	VertexID c = orient_tri_edge_and_find_other_vtx( a, b, m_vTriangles[t0].tv );
	TriangleID t1 = eAB->t[1];
	if (t1 != InvalidID) {
		VertexID d = find_tri_other_vtx( a, b, m_vTriangles[t1].tv );
		return Vector2i( c, d );
	} else
		return Vector2i(c, InvalidID);
}



DMESH_TEMPLATE
MeshResult DMESH_TYPE::GetVtxTriangles(VertexID vID, std::vector<int> & vTriangles, bool bUseOrientation) const
{
    if (! IsVertex(vID) )
        return MeshResult::Failed_NotAVertex;
    const edge_list & edges = m_vVertexEdges[vID];
    
    if ( bUseOrientation ) {
    
        for ( auto eid : edges ) {
            const Edge & e = m_vEdges[eid];
            VertexID vOther = e.otherV(vID);
            
            if ( m_vTriangles[e.t[0]].hasSequentialV(vID, vOther) )
                vTriangles.push_back(e.t[0]);
            if ( e.t[1] != InvalidID && m_vTriangles[e.t[1]].hasSequentialV(vID, vOther) )
                vTriangles.push_back(e.t[1]);
        }
        return MeshResult::Ok;
        
    }
    
    // brute-force method
    for ( auto eid : edges ) {
        const Edge & e = m_vEdges[eid];
        if ( std::find(vTriangles.begin(), vTriangles.end(), e.t[0]) == vTriangles.end() )
            vTriangles.push_back(e.t[0]);
        if ( e.t[1] != InvalidID && std::find(vTriangles.begin(), vTriangles.end(), e.t[1]) == vTriangles.end() )
            vTriangles.push_back(e.t[1]);
    }
    return MeshResult::Ok;
}

DMESH_TEMPLATE
EdgeID DMESH_TYPE::next_edge(TriangleID tID, VertexID vCenter, VertexID vOther) const
{
    const Triangle & t = m_vTriangles[tID];
    for ( int j = 0; j < 3; ++j ) {
        if ( t.tv[j] == vCenter ) {
            if ( t.tv[(j+1)%3] == vOther )
                return t.te[(j+2)%3];
            else if ( t.tv[(j+2)%3] == vOther )
                return t.te[j];
            return InvalidID;       // can't be any other edge!
        }
    }
    return InvalidID;
}

DMESH_TEMPLATE
MeshResult DMESH_TYPE::GetVtxOrderedTriangles(VertexID vID, std::vector<int> & vTriangles, bool bIsInteriorHint) const
{
    if (! IsVertex(vID) )
        return MeshResult::Failed_NotAVertex;
    
    const edge_list & vtx_edges = m_vVertexEdges[vID];
    
    // if we have a boundary edge we need to start there, so
    // without any information, we have to do a linear search
    EdgeID boundary_e[2] = {InvalidID,InvalidID};
    if ( bIsInteriorHint == false ) {
        for ( auto eid : vtx_edges ) {
            if ( m_vEdges[eid].isBoundary() ) {
                if ( boundary_e[0] == InvalidID )
                    boundary_e[0] = eid;
                else if ( boundary_e[1] == InvalidID )
                    boundary_e[1] = eid;
                else
                    return MeshResult::Failed_IsBowtieVertex;
            }
        }
    }
    
    EdgeID e_cur = InvalidID;
    VertexID vOther = InvalidID;
    TriangleID tNext = InvalidID;
    
    // figure out which edge/tri to start with. Prefer to walk CCW but this
    // is not always possible if tri orientation is inconsistent
    if ( boundary_e[0] == InvalidID ) {
        e_cur = vtx_edges.front();
        const Edge & e = m_vEdges[e_cur];
        vOther = e.otherV(vID);
        TriangleID tStart = ( find_tri_ordered_edge(vID, vOther, m_vTriangles[e.t[0]].tv) != InvalidID ) ? e.t[1] : e.t[0];
        vTriangles.push_back(tStart);
        tNext = e.otherT(tStart);
    } else {
        for ( int j = 0; j < 2; ++j ) {
            const Edge & e = m_vEdges[boundary_e[j]];
            vOther = e.otherV(vID);
            if ( j == 1 || find_tri_ordered_edge(vID, vOther, m_vTriangles[e.t[0]].tv) != InvalidID ) {
                tNext = e.t[0];
                e_cur = boundary_e[j];
                break;
            }
        }
    }
    
    // iterate until finished
    bool bDone = false;
    while (!bDone) {
        vTriangles.push_back(tNext);

        // find next edge
        e_cur = next_edge(tNext, vID, vOther);
        const Edge & e = m_vEdges[e_cur];

        // find next triangle
        vOther = e.otherV(vID);
        tNext = e.otherT(tNext);
        if ( tNext == vTriangles[0] || tNext == InvalidID )
            bDone = true;
    }
    
    return MeshResult::Ok;
}



DMESH_TEMPLATE
Vector3i DMESH_TYPE::GetTriTriangles(TriangleID tID) const
{
    if (! IsTriangle(tID) )
        return Vector3i(InvalidID, InvalidID, InvalidID);
    const Vector3i & te = m_vTriangles[tID].te;
    return Vector3i( m_vEdges[te[0]].otherT(tID),
                     m_vEdges[te[1]].otherT(tID),
                     m_vEdges[te[2]].otherT(tID) );
}




#define DMESH_CHECK_OR_FAIL(f) if (bAssert) { gDevAssert((f)); } else if ( (f) == false ) { return false; }

DMESH_TEMPLATE
bool DMESH_TYPE::CheckValidity(bool bAssert)
{
	std::vector<int> triToVtxRefs(GetMaxVertexID(), 0);

	for ( auto curt = m_vTriangles.begin(); curt != m_vTriangles.end(); ++curt ) {
		TriangleID tID = curt.index();
		DMESH_CHECK_OR_FAIL( IsTriangle(tID) );
		DMESH_CHECK_OR_FAIL( m_vTriangles.refCount(tID) == 1 );
		Triangle & tri = (*curt);

		// vertices must exist
		for (int j = 0; j < 3; ++j) {
			DMESH_CHECK_OR_FAIL( IsVertex( tri.tv[j] ) );
			triToVtxRefs[tri.tv[j]] += 1;
		}
        
		// edges must exist and reference this tri
		EdgeID e[3];
        for ( int j = 0; j < 3; ++j ) {
            auto a = tri.tv[j], b = tri.tv[(j+1)%3];
            e[j] = FindEdge(a,b);
            DMESH_CHECK_OR_FAIL( e[j] != InvalidID );
			DMESH_CHECK_OR_FAIL( m_vEdges[e[j]].hasT(tID) );
        }
		DMESH_CHECK_OR_FAIL( e[0] != e[1] && e[0] != e[2] && e[1] != e[2] );

		// tri nbrs must exist and reference this tri, or same edge must be boundary edge
		for (int j = 0; j < 3; ++j) {
			EdgeID eid = tri.te[j];
            DMESH_CHECK_OR_FAIL( IsEdge(eid) );
			Edge & e = m_vEdges[eid];
            TriangleID tOther = e.otherT(tID);
            if ( tOther == InvalidID ) {
                DMESH_CHECK_OR_FAIL( isBoundaryT(tID) );
                continue;
            }
            
            DMESH_CHECK_OR_FAIL( hasNeighbourT(tOther, tID) == true );

			// edge must have same two verts as tri for same index
			auto a = tri.tv[j], b = tri.tv[(j+1)%3];
			Vector2i ev = m_vEdges[tri.te[j]].v;
			DMESH_CHECK_OR_FAIL( same_pair_unordered(a, b, ev[0], ev[1]) );

            // also check that nbr edge has opposite orientation
            int found = find_tri_ordered_edge(b, a, m_vTriangles[tOther].tv);
            DMESH_CHECK_OR_FAIL(found != InvalidID);
		}
    }

	// edge verts/tris must exist
	for ( auto cure = m_vEdges.begin(); cure != m_vEdges.end(); ++cure ) {
		EdgeID eID = cure.index();
		DMESH_CHECK_OR_FAIL( IsEdge(eID) );
		DMESH_CHECK_OR_FAIL( m_vEdges.refCount(eID) == 1 );
		Edge & edge = (*cure);

		DMESH_CHECK_OR_FAIL( IsVertex(edge.v[0]) );
        DMESH_CHECK_OR_FAIL( IsVertex(edge.v[1]) );
		DMESH_CHECK_OR_FAIL( edge.v[0] < edge.v[1] );
        DMESH_CHECK_OR_FAIL( IsTriangle(edge.t[0]) );
        if ( edge.t[1] != InvalidID ) {
            DMESH_CHECK_OR_FAIL( IsTriangle(edge.t[1]) );
        }
    }

	// vertex edges must exist and reference this vert
    for ( auto curv = m_vVertices.begin(); curv != m_vVertices.end(); curv++ ) {
        VertexID vID = curv.index();
		DMESH_CHECK_OR_FAIL( IsVertex(vID) );
        edge_list & l = m_vVertexEdges[vID];
        for ( auto edgeid : l ) {
            DMESH_CHECK_OR_FAIL( IsEdge(edgeid) );
			DMESH_CHECK_OR_FAIL( m_vEdges[edgeid].hasV(vID) );

			VertexID otherV = m_vEdges[edgeid].otherV(vID);
			EdgeID e2;
			DMESH_CHECK_OR_FAIL( find_edge(vID, otherV, e2) != nullptr );
			DMESH_CHECK_OR_FAIL( e2 == edgeid );
			DMESH_CHECK_OR_FAIL( find_edge(otherV, vID, e2) != nullptr );
			DMESH_CHECK_OR_FAIL( e2 == edgeid );
		}

		std::vector<TriangleID> vTris;
		GetVtxTriangles(vID, vTris, false);
		DMESH_CHECK_OR_FAIL( m_vVertices.refCount(vID) == vTris.size() + 1 );
		DMESH_CHECK_OR_FAIL( triToVtxRefs[vID] == vTris.size() );
		for ( auto tID : vTris ) {
			DMESH_CHECK_OR_FAIL( m_vTriangles[tID].hasV(vID) );
		}
    }
    return true;
}




//----------------------------------------------------------------------------
// explicit instantiation
//----------------------------------------------------------------------------
namespace g3
{
template class DMesh3<float>;
template class DMesh3<double>;

}
//----------------------------------------------------------------------------
