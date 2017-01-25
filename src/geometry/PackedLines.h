#pragma once

#include <g3types.h>
#include <GeometryInterfaces.h>

namespace g3
{

class g3External PackedLines : public IPackedLines
{
public:
	PackedLines( LinesType eType = Segments );
	virtual ~PackedLines() {}

	/*
	* IPackedLines interface
	*/
	virtual LinesType GetLinesType() const;

	virtual unsigned int GetVertexCount() const;
	virtual unsigned int GetLineCount() const;

	virtual bool HasPositions() const;
	virtual const float * GetPositionsBuffer() const;

	virtual bool HasColorsFloat() const;
	virtual const float * GetColorsFloatBuffer() const;

	virtual bool HasIndices() const;
	virtual const unsigned int * GetIndicesBuffer() const;


	/*
	 * construction
	 */
	void ResizeVertices(unsigned int nVertices, bool bNormals, bool bColors);
	void ResizeLines(unsigned int nTriangles);

	void CopyVertices(const std::vector<float> & vPositions, 
					  const std::vector<float> * pColors = nullptr );
	void CopyVertices(const float * pPositions, size_t nVertices,
					  const float * pColors = nullptr );


	void CopyLines( const std::vector<unsigned int> & vIndices );
	void CopyLines( const unsigned int * pIndices, size_t nLines );
	void CopyLines( const int * pIndices, size_t nLines );


	/*
	 * Utility
	 */
	const float * GetPosition(unsigned int k) const { 
		return &m_vPositions[3*k]; 
	}
	const float * GetColor(unsigned int k) const { 
		return &m_vColors[3*k]; 
	}
	void GetSegmentIndices(unsigned int k, unsigned int & i0, unsigned int & i1) const;
	void GetSegment(unsigned int k, g3::Vector3f & v0, g3::Vector3f & v1) const;

protected:
	LinesType m_eType;

	std::vector<float> m_vPositions;
	std::vector<float> m_vColors;

	std::vector<unsigned int> m_vIndices;
};


}



