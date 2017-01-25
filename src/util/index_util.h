#pragma once


namespace g3 {


// test if [a0,a1] and [b0,b1] are the same pair, ignoring order
template<typename T>
inline bool same_pair_unordered( T a0, T a1, T b0, T b1 )
{
	return ( a0 == b0 ) ? 
		(a1 == b1) : 
		(a0 == b1 && a1 == b0);
}

// return index of a in tri_verts, or InvalidID if not found
template<typename T, typename Vec>
inline int find_tri_index( T a, Vec tri_verts )
{
	if ( tri_verts[0] == a ) return 0;
	if ( tri_verts[1] == a ) return 1;
	if ( tri_verts[2] == a ) return 2;
	return InvalidID;
}

// return index of a in tri_verts, or InvalidID if not found
template<typename T, typename Vec>
inline int find_edge_index_in_tri( T a, T b, Vec & tri_verts )
{
	if ( same_pair_unordered(a, b, tri_verts[0], tri_verts[1]) ) return 0;
	if ( same_pair_unordered(a, b, tri_verts[1], tri_verts[2]) ) return 1;
	if ( same_pair_unordered(a, b, tri_verts[2], tri_verts[0]) ) return 2;
	return InvalidID;
}


// find sequence [a,b] in tri_verts (mod3) and return index of a, or InvalidID if not found
template<typename T, typename Vec>
inline int find_tri_ordered_edge( T a, T b, Vec tri_verts )
{
	if ( tri_verts[0] == a && tri_verts[1] == b )	return 0;
	if ( tri_verts[1] == a && tri_verts[2] == b )	return 1;
	if ( tri_verts[2] == a && tri_verts[0] == b )	return 2;
    return InvalidID;
}

// find sequence [a,b] in tri_verts (mod3) then return the third **value**, or InvalidID if not found
template<typename T, typename Vec>
inline int find_tri_other_vtx( T a, T b, Vec tri_verts )
{
    for (int j = 0; j < 3; ++j ) {
        if ( same_pair_unordered(a, b, tri_verts[j], tri_verts[(j+1)%3]) )
            return tri_verts[(j+2)%3];
    }
    return InvalidID;
}

// find sequence [a,b] in tri_verts (mod3) then return the third **index**, or InvalidID if not found
template<typename T, typename Vec>
inline int find_tri_other_index( T a, T b, Vec tri_verts )
{
	for (int j = 0; j < 3; ++j ) {
		if ( same_pair_unordered(a, b, tri_verts[j], tri_verts[(j+1)%3]) )
			return (j+2)%3;
	}
	return InvalidID;
}
   
// set [a,b] to order found in tri_verts (mod3), and return third **value**, or InvalidID if not found
template<typename T, typename Vec>
inline int orient_tri_edge_and_find_other_vtx( T & a, T & b, Vec tri_verts )
{
    for ( int j = 0; j < 3; ++j ) {
        if ( same_pair_unordered(a, b, tri_verts[j], tri_verts[(j+1)%3]) ) {
            a = tri_verts[j];
            b = tri_verts[(j+1)%3];
            return tri_verts[(j+2)%3];
        }
    }
    return InvalidID;
}
    
    
// set v[i] = mapper(v[i])
template<typename Func>
void apply_map( Vector3i & v, Func mapper )
{
	v[0] = mapper[v[0]];
	v[1] = mapper[v[1]];
	v[2] = mapper[v[2]];
}

// return v[i] = mapper(v[i])
template<typename Func>
Vector3i apply_map( const Vector3i & vIn, Func mapper )
{
	return Vector3i( mapper[vIn[0]], mapper[vIn[1]], mapper[vIn[2]] );
}

// set v[i] = mapper(v[i])
template<typename T, typename Func>
void apply_map( Vector3<T> & v, Func mapper )
{
	v[0] = mapper[v[0]];
	v[1] = mapper[v[1]];
	v[2] = mapper[v[2]];
}

// return v[i] = mapper(v[i])
template<typename T, typename Func>
Vector3<T> apply_map( const Vector3<T> & vIn, Func mapper )
{
	return Vector3<T>( mapper[vIn[0]], mapper[vIn[1]], mapper[vIn[2]] );
}




} // end namespace g3