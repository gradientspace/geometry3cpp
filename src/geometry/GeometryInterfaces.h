#pragma once

#include <g3types.h>
#include <BaseInterfaces.h>


namespace g3
{

class g3External IPackedLines : public ITimeStamped
{
public:
	virtual ~IPackedLines() {}

	//
	// IPackedLines interface. This interface supports various different styles
	// of line storage (via LinesType). Indices are only expected to be available
	// if the type is IndexedSegments, otherwise segments are formed from
	// the ordering of the vertices. 
	//

	enum LinesType {
		Segments = 0,
		IndexedSegments = 1,
		Strip = 2,
		Loop = 3
	};
	virtual LinesType GetLinesType() const { return Segments; }

	virtual unsigned int GetVertexCount() const  = 0;
	virtual unsigned int GetLineCount() const = 0;

	virtual bool HasPositions() const  { return true; }
	virtual const float * GetPositionsBuffer() const  = 0;

	virtual bool HasColorsFloat() const  { return false; }
	virtual const float * GetColorsFloatBuffer() const  { return nullptr; }

	virtual bool HasIndices() const { return false; }
	virtual const unsigned int * GetIndicesBuffer() const { return nullptr; }
};




class g3External IPackedMesh : public ITimeStamped
{
public:
	virtual ~IPackedMesh() {}

	//
	// IPackedMesh interface
	//

	virtual unsigned int GetVertexCount() const = 0;
	virtual unsigned int GetTriangleCount() const = 0;

	virtual bool HasPositions() const { return true; }
	virtual const float * GetPositionsBuffer() const = 0;

	virtual bool HasNormals() const { return true; }
	virtual const float * GetNormalsBuffer() const { return nullptr; }

	virtual bool HasColorsFloat() const { return false; }
	virtual const float * GetColorsFloatBuffer() const  { return nullptr; }

	virtual bool HasIndices() const { return true; }
	virtual const unsigned int * GetIndicesBuffer() const = 0;

	virtual bool HasTriColors() const { return false; }
	virtual const float * GetTriColorsBuffer() const { return nullptr; }

	virtual bool HasTriGroups() const { return false; }
	virtual const unsigned int * GetTriGroupsBuffer() const { return nullptr; }
};



template<typename Real>
class IDynamicMesh : public ITimeStamped
{
public:
	virtual ~IDynamicMesh() = default;

	//
	// IDynamicMesh interface
	//
	virtual unsigned int GetVertexCount() const = 0;
	virtual unsigned int GetTriangleCount() const = 0;

	virtual VertexID AppendVertex( const Vector3<Real> & v ) = 0;
	virtual TriangleID AppendTriangle( const Vector3i & t, GroupID gID = 0) = 0;

	// [RMS] this is mainly used for debugging, so is not strictly necessary to implement
	virtual bool CheckValidity(bool bAssert) { return true; }
};


//
//IGeometryCollector interface is used to collect up geometric sub-elements
//  from higher-level classes.
//
class IGeometryCollector
{
public:
	IGeometryCollector() = default;
	virtual ~IGeometryCollector() = default;

	//
	// IGeometryCollector interface
	//
	virtual void BeginCollection() {}
	virtual void SetCurrentLayer( int nLayer ) {}
	virtual void AddGeometry( IPackedMesh * pMesh ) = 0;
	virtual void AddGeometry( IPackedLines * pLines ) = 0;
	virtual void AddGeometry( ITimeStamped * pAny ) {}
	virtual void EndCollection() {}
};



//
//ICollectable interface provides standard API for collecting geometry
// from objects (for rendering/etc)
//  
class ICollectable
{
public:
	ICollectable() = default;
	virtual ~ICollectable() = default;

	//
	// ICollectable interface
	//
	virtual void Collect(IGeometryCollector * pCollector) = 0;
};


}

