#pragma once

#include <GeometryInterfaces.h>
#include <DMesh3.h>

namespace g3 {

template<typename Real>
class DSurface
{
public:
	typedef DMesh3<Real, g3::dvector> MeshType;

	DSurface() {
		m_pMesh = std::shared_ptr<MeshType>(new MeshType);
	}
	virtual ~DSurface() {}

	MeshType & Mesh() { return *m_pMesh; }
	const MeshType & Mesh() const { return *m_pMesh; }

	void Collect( g3::IGeometryCollector * pCollector ) {
		pCollector->AddGeometry( m_pMesh.get() );
	}

protected:
	std::shared_ptr<MeshType> m_pMesh;
};

typedef std::shared_ptr<DSurface<float>> DSurfacePtr;


}   // end namespace g3

