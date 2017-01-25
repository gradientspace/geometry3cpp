#include <geometry3PCH.h>
#include "PackedLines.h"

#include "g3Debug.h"
#include <g3platform.h>



using namespace g3;

PackedLines::PackedLines( LinesType eType )
{
	m_eType = eType;
}

PackedLines::LinesType g3::PackedLines::GetLinesType() const
{
	return m_eType;
}

unsigned int PackedLines::GetVertexCount() const
{
	return (unsigned int)m_vPositions.size() / 3;
}

unsigned int PackedLines::GetLineCount() const
{
	switch (m_eType) {
		default:
		case Segments:
			return (unsigned int)m_vPositions.size() / 2;
		case Strip:
			return (m_vPositions.empty()) ? 0 : (unsigned int)m_vPositions.size() - 1;
		case Loop:
			return (unsigned int)m_vPositions.size();
		case IndexedSegments:
			return (unsigned int)m_vIndices.size() / 2;
	}
}

bool PackedLines::HasPositions() const
{
	return ! m_vPositions.empty();
}
const float * PackedLines::GetPositionsBuffer() const
{
	return HasPositions() ? &m_vPositions[0] : nullptr;
}



bool PackedLines::HasColorsFloat() const
{
	return HasPositions() && m_vPositions.size() == m_vColors.size();
}
const float * PackedLines::GetColorsFloatBuffer() const
{
	return HasColorsFloat() ? &m_vColors[0] : nullptr;
}


bool PackedLines::HasIndices() const
{
	return (m_eType == IndexedSegments) && (! m_vIndices.empty());
}
const unsigned int * PackedLines::GetIndicesBuffer() const
{
	return HasIndices() ? &m_vIndices[0] : nullptr;
}


void PackedLines::ResizeVertices( unsigned int nVertices, bool bNormals, bool bColors )
{
	m_vPositions.resize(nVertices * 3);
	if ( bColors )
		m_vColors.resize( m_vColors.size() );
	updateTimeStamp();
}

void PackedLines::ResizeLines( unsigned int nLines )
{
	m_vIndices.resize(nLines * 2);
	updateTimeStamp();
}

void PackedLines::CopyVertices( const std::vector<float>& vPositions, const std::vector<float>* pColors )
{
	m_vPositions = vPositions;
	if (pColors) {
		gDevAssert(pColors->size() == m_vPositions.size());
		m_vColors = *pColors;
	}
	updateTimeStamp();
}

void PackedLines::CopyVertices( const float * pPositions, size_t nVertices, const float * pColors )
{
	m_vPositions.resize( nVertices*3 );
	memcpy_s( &m_vPositions[0], nVertices*3*sizeof(float), pPositions, nVertices*3*sizeof(float) );

	if (pColors) {
		m_vColors.resize( nVertices*3 );
		memcpy_s( &m_vColors[0], nVertices*3*sizeof(float), pColors, nVertices*3*sizeof(float) );
	}
	updateTimeStamp();
}


void PackedLines::CopyLines( const std::vector<unsigned int> & vIndices )
{
	gDevAssert(m_eType == Segments || m_eType == IndexedSegments );
	m_vIndices = vIndices;
	m_eType = IndexedSegments;
	updateTimeStamp();
}

void PackedLines::CopyLines( const unsigned int * pIndices, size_t nLines )
{
	gDevAssert(m_eType == Segments || m_eType == IndexedSegments );
	m_vIndices.resize( nLines*2 );
	memcpy_s( &m_vIndices[0], nLines*2*sizeof(unsigned int), pIndices, nLines*2*sizeof(unsigned int) );
	m_eType = IndexedSegments;
	updateTimeStamp();
}

void PackedLines::CopyLines( const int * pIndices, size_t nLines )
{
	CopyLines( (unsigned int *)pIndices, nLines );
	updateTimeStamp();
}


void PackedLines::GetSegmentIndices( unsigned int k, unsigned int & i0, unsigned int & i1 ) const
{
	i0 = k; 
	i1 = (k+1) % m_vPositions.size();
	if (m_eType == Segments) {
		i0 = 2*k;
		i1 = 2*k + 1;
	} else if (m_eType == IndexedSegments) {
		i0 = m_vIndices[2*k];
		i1 = m_vIndices[2*k+1];
	}
}


void PackedLines::GetSegment( unsigned int k, g3::Vector3f & v0, g3::Vector3f & v1 ) const
{
	unsigned int i0,i1;
	GetSegmentIndices(k, i0, i1);
	v0 = g3::Vector3f(&m_vPositions[i0]);
	v1 = g3::Vector3f(&m_vPositions[i1]);
}
