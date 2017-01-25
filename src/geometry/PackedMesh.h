#pragma once

#include <g3types.h>
#include <GeometryInterfaces.h>

namespace g3
{


/*
 * Simple mesh class where data is stored in internal buffers.
 * Colors are float-RGB (enforced internally w/ asserts, etc)
 */
class g3External PackedMesh : public IPackedMesh
{
public:
	virtual ~PackedMesh() {}


	/*
	* IPackedMesh interface
	*/

	virtual unsigned int GetVertexCount() const override;
	virtual unsigned int GetTriangleCount() const override;

	virtual bool HasPositions() const override;
	virtual const float * GetPositionsBuffer() const override;

	virtual bool HasNormals() const override;
	virtual const float * GetNormalsBuffer() const override;

	virtual bool HasColorsFloat() const override;
	virtual const float * GetColorsFloatBuffer() const override;

	virtual bool HasIndices() const override;
	virtual const unsigned int * GetIndicesBuffer() const override;

	virtual bool HasTriColors() const override;
	virtual const float * GetTriColorsBuffer() const override;

	virtual bool HasTriGroups() const override;
	virtual const unsigned int * GetTriGroupsBuffer() const override;


	/*
	 * construction
	 */
	void ResizeVertices(unsigned int nVertices, bool bNormals, bool bColors);
	void ResizeTriangles(unsigned int nTriangles);

	void CopyVertices(const std::vector<float> & vPositions, 
					  const std::vector<float> * pNormals = nullptr,
					  const std::vector<float> * pColors = nullptr );
	void CopyVertices(const float * pPositions, size_t nVertices,
					  const float * pNormals = nullptr,
					  const float * pColors = nullptr );


	void CopyTriangles( const std::vector<unsigned int> & vIndices,
						const std::vector<unsigned int> * pGroups = nullptr );
	void CopyTriangles( const unsigned int * pIndices, size_t nTriangles, 
						const unsigned int * pGroups = nullptr );
	void CopyTriangles( const int * pIndices, size_t nTriangles, 
						const int * pGroups = nullptr );
	void CopyTriColors( const std::vector<float> * pColors );
	void CopyTriColors( const float * pColors, size_t nTriangles );
	void CopyTriGroups( const std::vector<unsigned int> * pGroups );
	void CopyTriGroups( const unsigned int * pGroups, size_t nTriangles );

	void SetConstantTriGroup( unsigned int nGroupID );

	/*
	 * Utility
	 */
	float * GetPosition(unsigned int k) { 
		return &m_vPositions[3*k]; 
	}
	const float * GetPosition(unsigned int k) const { 
		return &m_vPositions[3*k]; 
	}
	float * GetNormal(unsigned int k) { 
		return &m_vNormals[3*k]; 
	}
	const float * GetNormal(unsigned int k) const { 
		return &m_vNormals[3*k]; 
	}
	float * GetColor(unsigned int k) { 
		return &m_vColors[3*k]; 
	}
	const float * GetColor(unsigned int k) const { 
		return &m_vColors[3*k]; 
	}
	float * GetTriColor(unsigned int k) { 
		return &m_vTriColors[3*k]; 
	}
	const float * GetTriColor(unsigned int k) const { 
		return &m_vTriColors[3*k]; 
	}
	unsigned int * GetTriGroup(unsigned int k) { 
		return &m_vTriGroups[k]; 
	}
	unsigned int GetTriGroup(unsigned int k) const { 
		return m_vTriGroups[k]; 
	}

	void GetTriangle(unsigned int k, unsigned int & i0, unsigned int & i1, unsigned int & i2) const;
	void GetTriangle(unsigned int k, unsigned int  * iTri) const;
	void GetTriangle(unsigned int k, g3::Vector3f & v0, g3::Vector3f & v1, g3::Vector3f & v2) const;
	void GetTriangle(unsigned int k, g3::Vector3f * pTri) const;
	g3::Vector3f GetTriangleNormal(unsigned int k) const;

	void EstimateNormals();

protected:
	std::vector<float> m_vPositions;
	std::vector<float> m_vNormals;
	std::vector<float> m_vColors;
	std::vector<unsigned int> m_vIndices;
	std::vector<float> m_vTriColors;
	std::vector<unsigned int> m_vTriGroups;
};


}



