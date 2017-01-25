#include <geometry3PCH.h>
#include "MeshUtil.h"

#include <VectorUtil.h>
#include <Wm5PolyhedralMassProperties.h>
#include <Wm5ContBox3.h>

#include <algorithm>
#include <limits>
#include <atomic_util.h>


using namespace g3;

template<class Real>
const Real * get_value( const Real * pBuffer, unsigned int nIndex, unsigned int nStride = 3 ) {
	return & pBuffer[nIndex*nStride];
}




double g3::Area( const IPackedMesh * pMesh )
{
	// was using this in tbb version but we removed tbb...
	//atomic_accumulator<double> area(0.0);

	const float * pVertices = pMesh->GetPositionsBuffer();
	const unsigned int * pTriangles = pMesh->GetIndicesBuffer();
	unsigned int nTriangles = pMesh->GetTriangleCount();
 
	double area = 0;
	for ( unsigned int ti = 0; ti < nTriangles; ++ti ) {
		auto pTri = get_value( pTriangles, (unsigned int)ti );
		double a = Area( f2d( Vector3f( get_value( pVertices, pTri[0] ) ) ),
							f2d( Vector3f( get_value( pVertices, pTri[1] ) ) ),
							f2d( Vector3f( get_value( pVertices, pTri[2] ) ) ) );
		area += a;
	}

	return area;
}


class PackedMeshVertexSource : public Wml::VertexSource<double>
{
public:
	const float * pVertices;
	virtual ~PackedMeshVertexSource() = default;
	virtual Vector3<double> operator[]( unsigned int i ) const {
		i = 3*i;
		return Vector3<double>( pVertices[i], pVertices[i+1], pVertices[i+2] );
	}
};


double g3::Volume( const IPackedMesh * pMesh )
{
	PackedMeshVertexSource tmp;
	tmp.pVertices = pMesh->GetPositionsBuffer();

	double fMass; Vector3d vCoM; Matrix3d vMoI;
	Wml::ComputeMassProperties( &tmp, 
								(int)pMesh->GetTriangleCount(), (const int *)pMesh->GetIndicesBuffer(), 
								false, fMass, vCoM, vMoI );
	return fMass;
}


Box3f g3::AxisBoundingBox( const IPackedMesh * pMesh )
{
	return Wml::ContAlignedBox( (int)pMesh->GetVertexCount(), (const Vector3f *)pMesh->GetPositionsBuffer() );
}

Box3f g3::OrientedBoundingBox( const IPackedMesh * pMesh, bool bFast )
{
	return Wml::ContOrientedBox( (int)pMesh->GetVertexCount(), (const Vector3f *)pMesh->GetPositionsBuffer() );
}




void g3::Translate( IPackedMesh * pMesh, const Vector3f & vTranslation )
{
	// AAAAHHHH 
	float * pVertices = const_cast<float *>( pMesh->GetPositionsBuffer() );
	unsigned int nVertices = pMesh->GetVertexCount();

	for ( unsigned int vi = 0; vi < nVertices; ++vi ) {
		unsigned int k = vi*3;
		pVertices[k] += vTranslation[0];
		pVertices[k+1] += vTranslation[1];
		pVertices[k+2] += vTranslation[2];
	}

	pMesh->updateTimeStamp();
}


void g3::Scale( IPackedMesh * pMesh, const Vector3f & vScale )
{
	// AAAAHHHH 
	float * pVertices = const_cast<float *>( pMesh->GetPositionsBuffer() );
	unsigned int nVertices = pMesh->GetVertexCount();

	for ( unsigned int vi = 0; vi < nVertices; ++vi ) {
		unsigned int k = vi*3;
		pVertices[k] *= vScale[0];
		pVertices[k+1] *= vScale[1];
		pVertices[k+2] *= vScale[2];
	}

	pMesh->updateTimeStamp();
}

namespace g3
{

}

