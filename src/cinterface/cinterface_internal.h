#pragma once

#include <g3types.h>
#include <PackedMesh.h>
#include <DMesh3.h>

#include "cinterface.h"

using namespace g3;


struct MeshHandle
{
	static const int VALID = 0xdeadbeef;
	static const int RELEASED = 0xf00ff00f;

	int valid_token;

	PackedMesh * packed;
	DMesh3 * dynamic;

	MeshHandle( PackedMesh * pPacked ) {
		packed = pPacked;
		dynamic = nullptr;
		valid_token = VALID;
	}

	MeshHandle( DMesh3 * pDynamic ) {
		packed = nullptr;
		dynamic = pDynamic;
		valid_token = VALID;
	}

	void Release() {
		if ( valid_token != 0xdeadbeef )
			return;
		valid_token = 0xf00ff00f;
		if ( packed != nullptr )
			delete packed;
		if ( dynamic != nullptr )
			delete dynamic;
	}


	int is_valid() {
		return (valid_token == VALID  && ( packed != nullptr || dynamic != nullptr))  ? 1 : 0;
	}


	int vertex_count() {
		if ( is_valid() && packed != nullptr )
			return packed->GetVertexCount();
		else if ( is_valid() && dynamic != nullptr )
			return dynamic->VertexCount();
		return -1;
	}

	int triangle_count() {
		if ( is_valid() && packed != nullptr )
			return packed->GetTriangleCount();
		else if ( is_valid() && dynamic != nullptr )
			return dynamic->TriangleCount();
		return -1;
	}


	int get_vertices_float( int buffer_size, float * buf ) {
		int NV = vertex_count();
		if ( NV <= 0 || NV*3 > buffer_size )
			return -1;
		if (packed != nullptr) {
			memcpy_s( buf, NV*3*sizeof( float ), packed->GetPositionsBuffer(), NV*3*sizeof( float ) );
			return NV;
		} else if (dynamic != nullptr) {
			int k = 0;
			for ( int vid = 0; vid < NV; ++vid ) {
				Vector3d v = dynamic->GetVertex(vid);
				buf[k++] = (float)v[0]; buf[k++] = (float)v[1]; buf[k++] = (float)v[2];
			}
			return NV;			
		}
		return -1;
	}

	int get_vertices_double( int buffer_size, double * buf ) {
		int NV = vertex_count();
		if ( NV <= 0 || NV*3 > buffer_size )
			return -1;
		int k = 0;
		if (packed != nullptr) {
			for ( int vid = 0; vid < NV; ++vid ) {
				Vector3f v = Vector3f(packed->GetPosition(vid));
				buf[k++] = v[0]; buf[k++] = v[1]; buf[k++] = v[2];
			}
			return NV;
		} else if (dynamic != nullptr) {
			for ( int vid = 0; vid < NV; ++vid ) {
				Vector3d v = dynamic->GetVertex(vid);
				buf[k++] = v[0]; buf[k++] = v[1]; buf[k++] = v[2];
			}
			return NV;			
		}
		return -1;
	}


	int get_triangles( int buffer_size, unsigned int * buf ) {
		int NT = triangle_count();
		if ( NT <= 0 || NT*3 > buffer_size )
			return -1;
		if (packed != nullptr) {
			memcpy_s( buf, NT*3*sizeof(unsigned int), packed->GetIndicesBuffer(), NT*3*sizeof(unsigned int) );
			return NT;
		} else if (dynamic != nullptr) {
			int k = 0;
			for ( int tid = 0; tid < NT; ++tid ) {
				Vector3i t = dynamic->GetTriangle(tid);
				buf[k++] = t[0]; buf[k++] = t[1]; buf[k++] = t[2];
			}
			return NT;			
		}
		return -1;
	}


};
