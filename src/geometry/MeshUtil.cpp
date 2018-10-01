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
	virtual Wml::Vector3<double> operator[]( unsigned int i ) const {
		i = 3*i;
		return Wml::Vector3<double>( pVertices[i], pVertices[i+1], pVertices[i+2] );
	}
};


double g3::Volume( const IPackedMesh * pMesh )
{
	PackedMeshVertexSource tmp;
	tmp.pVertices = pMesh->GetPositionsBuffer();

	double fMass; Wml::Vector3d vCoM; Wml::Matrix3d vMoI;
	Wml::ComputeMassProperties( &tmp, 
								(int)pMesh->GetTriangleCount(), (const int *)pMesh->GetIndicesBuffer(), 
								false, fMass, vCoM, vMoI );
	return fMass;
}


Box3f g3::AxisBoundingBox( const IPackedMesh * pMesh )
{
	int NV = pMesh->GetVertexCount();
	std::vector<Wml::Vector3f> pts;
	pts.reserve(NV);
	const float * p = pMesh->GetPositionsBuffer();
	for (int k = 0; k < NV; ++k)
		pts.push_back( Wml::Vector3f(p[3*k], p[3*k+1], p[3*k+2]) );
	return Wml::ContAlignedBox(NV, &pts[0]);
}

Box3f g3::OrientedBoundingBox( const IPackedMesh * pMesh, bool bFast )
{
	int NV = pMesh->GetVertexCount();
	std::vector<Wml::Vector3f> pts;
	pts.reserve(NV);
	const float * p = pMesh->GetPositionsBuffer();
	for (int k = 0; k < NV; ++k)
		pts.push_back(Wml::Vector3f(p[3 * k], p[3 * k + 1], p[3 * k + 2]));
	return Wml::ContOrientedBox(NV, &pts[0]);
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

