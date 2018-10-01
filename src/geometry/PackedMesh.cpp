#include <geometry3PCH.h>
#include "PackedMesh.h"

#include <VectorUtil.h>
#include <g3Debug.h>
#include <g3platform.h>

using namespace g3;

unsigned int g3::PackedMesh::GetVertexCount() const
{
	return (unsigned int)m_vPositions.size() / 3;
}

unsigned int g3::PackedMesh::GetTriangleCount() const
{
	return (unsigned int)m_vIndices.size() / 3;
}

bool g3::PackedMesh::HasPositions() const
{
	return ! m_vPositions.empty();
}

const float * g3::PackedMesh::GetPositionsBuffer() const
{
	return HasPositions() ? &m_vPositions[0] : nullptr;
}

bool g3::PackedMesh::HasNormals() const
{
	return HasPositions() && m_vNormals.size() == m_vPositions.size();
}

const float * g3::PackedMesh::GetNormalsBuffer() const
{
	return HasNormals() ? &m_vNormals[0] : nullptr;
}

bool g3::PackedMesh::HasColorsFloat() const
{
	return HasPositions() && m_vPositions.size() == m_vColors.size();
}

const float * g3::PackedMesh::GetColorsFloatBuffer() const
{
	return HasColorsFloat() ? &m_vColors[0] : nullptr;
}

bool g3::PackedMesh::HasIndices() const
{
	return ! m_vIndices.empty();
}

const unsigned int * g3::PackedMesh::GetIndicesBuffer() const
{
	return HasIndices() ? &m_vIndices[0] : nullptr;
}


bool PackedMesh::HasTriColors() const
{
	return ! m_vTriColors.empty();
}
const float * PackedMesh::GetTriColorsBuffer() const
{
	return HasTriColors() ? &m_vTriColors[0] : nullptr;
}


bool PackedMesh::HasTriGroups() const
{
	return ! m_vTriGroups.empty();
}
const unsigned int * PackedMesh::GetTriGroupsBuffer() const
{
	return HasTriGroups() ? &m_vTriGroups[0] : nullptr;
}


void g3::PackedMesh::ResizeVertices( unsigned int nVertices, bool bNormals, bool bColors )
{
	m_vPositions.resize(nVertices * 3);
	if ( bNormals )
		m_vNormals.resize( m_vPositions.size() );
	if ( bColors )
		m_vColors.resize( m_vColors.size() );
	updateTimeStamp();
}

void g3::PackedMesh::ResizeTriangles( unsigned int nTriangles )
{
	m_vIndices.resize(nTriangles * 3);
	updateTimeStamp();
}

void g3::PackedMesh::CopyVertices( const std::vector<float>& vPositions, const std::vector<float>* pNormals, const std::vector<float>* pColors )
{
	m_vPositions = vPositions;
	if (pNormals) {
		gDevAssert(pNormals->size() == m_vPositions.size());
		m_vNormals = *pNormals;
	}
	if (pColors) {
		gDevAssert(pColors->size() == m_vPositions.size());
		m_vColors = *pColors;
	}
	updateTimeStamp();
}

void g3::PackedMesh::CopyVertices( const float * pPositions, size_t nVertices, const float * pNormals, const float * pColors )
{
	m_vPositions.resize( nVertices*3 );
	memcpy_s( &m_vPositions[0], nVertices*3*sizeof(float), pPositions, nVertices*3*sizeof(float) );

	if (pNormals) {
		m_vNormals.resize( nVertices*3 );
		memcpy_s( &m_vNormals[0], nVertices*3*sizeof(float), pNormals, nVertices*3*sizeof(float) );
	}
	if (pColors) {
		m_vColors.resize( nVertices*3 );
		memcpy_s( &m_vColors[0], nVertices*3*sizeof(float), pColors, nVertices*3*sizeof(float) );
	}
	updateTimeStamp();
}

void g3::PackedMesh::CopyTriangles( const std::vector<unsigned int> & vIndices,
									const std::vector<unsigned int> * pGroups )
{
	m_vIndices = vIndices;
	if ( pGroups )
		m_vTriGroups = *pGroups;
	updateTimeStamp();
}

void g3::PackedMesh::CopyTriangles( const unsigned int * pIndices, size_t nTriangles, const unsigned int * pGroups )
{
	m_vIndices.resize( nTriangles*3 );
	memcpy_s( &m_vIndices[0], nTriangles*3*sizeof(unsigned int), pIndices, nTriangles*3*sizeof(unsigned int) );
	if (pGroups) {
		m_vTriGroups.resize( nTriangles );
		memcpy_s( &m_vTriGroups[0], nTriangles*sizeof(unsigned int), pGroups, nTriangles*sizeof(unsigned int) );
	}
	updateTimeStamp();
}

void g3::PackedMesh::CopyTriangles( const int * pIndices, size_t nTriangles, const int * pGroups )
{
	CopyTriangles( (unsigned int *)pIndices, nTriangles, (unsigned int *)pGroups );
	updateTimeStamp();
}


void PackedMesh::CopyTriColors( const std::vector<float> * pColors )
{
	gDevAssert(pColors->size() == GetTriangleCount() * 3);
	m_vTriColors = *pColors;
	updateTimeStamp();
}
void PackedMesh::CopyTriColors( const float * pColors, size_t nTriangles )
{
	gDevAssert(nTriangles == GetTriangleCount());
	m_vTriColors.resize( nTriangles*3 );
	memcpy_s( &m_vTriColors[0], nTriangles*3*sizeof(float), pColors, nTriangles*3*sizeof(float) );
	updateTimeStamp();
}


void PackedMesh::CopyTriGroups( const std::vector<unsigned int> * pGroups )
{
	gDevAssert(pGroups->size() == GetTriangleCount());
	m_vTriGroups = *pGroups;
	updateTimeStamp();
}
void PackedMesh::CopyTriGroups( const unsigned int * pGroups, size_t nTriangles )
{
	gDevAssert(nTriangles == GetTriangleCount());
	m_vTriGroups.resize( nTriangles );
	memcpy_s( &m_vTriGroups[0], nTriangles*sizeof(unsigned int), pGroups, nTriangles*sizeof(unsigned int) );
	updateTimeStamp();
}

void g3::PackedMesh::SetConstantTriGroup( unsigned int nGroupID )
{
	m_vTriGroups.resize( GetTriangleCount() );
	std::fill_n(m_vTriGroups.begin(), GetTriangleCount(), nGroupID);
	updateTimeStamp();
}


void PackedMesh::GetTriangle( unsigned int k, unsigned int & i0, unsigned int & i1, unsigned int & i2 ) const
{
	const unsigned int * t = &m_vIndices[3*k];
	i0 = t[0];
	i1 = t[1];
	i2 = t[2];
}
void PackedMesh::GetTriangle( unsigned int k, unsigned int * iTri ) const
{
	const unsigned int * t = &m_vIndices[3*k];
	iTri[0] = t[0];
	iTri[1] = t[1];
	iTri[2] = t[2];
}


void PackedMesh::GetTriangle( unsigned int k, g3::Vector3f & v0, g3::Vector3f & v1, g3::Vector3f & v2 ) const
{
	unsigned int i0, i1, i2;
	GetTriangle(k, i0,i1,i2);
	v0 = g3::Vector3f(GetPosition(i0));
	v1 = g3::Vector3f(GetPosition(i1));
	v2 = g3::Vector3f(GetPosition(i2));
}
void PackedMesh::GetTriangle( unsigned int k, g3::Vector3f * pTri ) const
{
	unsigned int i0, i1, i2;
	GetTriangle(k, i0,i1,i2);
	pTri[0] = g3::Vector3f(GetPosition(i0));
	pTri[1] = g3::Vector3f(GetPosition(i1));
	pTri[2] = g3::Vector3f(GetPosition(i2));
}

g3::Vector3f PackedMesh::GetTriangleNormal( unsigned int k ) const
{
	unsigned int i0, i1, i2;
	GetTriangle(k, i0,i1,i2);
	return g3::Normal( g3::Vector3f(GetPosition(i0)), g3::Vector3f(GetPosition(i1)), g3::Vector3f(GetPosition(i2)) );
}

void PackedMesh::EstimateNormals()
{
	m_vNormals.resize(0);
	m_vNormals.resize(m_vPositions.size(), 0.0f);
	float * pNormals = &m_vNormals[0];

	unsigned int nTriangles = GetTriangleCount();
	std::vector<g3::Vector3f> vTriNormals(nTriangles);
	for (unsigned int k = 0; k < nTriangles; ++k) 
		vTriNormals[k] = GetTriangleNormal(k);

	for (unsigned int k = 0; k < nTriangles; ++k) {
		unsigned int nTri[3];
		GetTriangle( k, nTri );
		array3f_add( pNormals, nTri[0], vTriNormals[k].data() );
		array3f_add( pNormals, nTri[1], vTriNormals[k].data() );
		array3f_add( pNormals, nTri[2], vTriNormals[k].data() );
	}

	unsigned int nVertices = GetVertexCount();
	for ( unsigned int k = 0; k < nVertices; ++k )
		array3f_normalize( pNormals, k );

	updateTimeStamp();
}