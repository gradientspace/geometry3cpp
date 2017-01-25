#pragma once

#include <DMesh3.h>

namespace g3 
{

template<typename Real>
class MeshSubdivider
{
public:
	typedef DMesh3<Real, g3::dvector> MeshType;

	MeshSubdivider();
	virtual ~MeshSubdivider() = default;

	virtual void Split1to4( DMesh3<Real> & mesh );

	// [TODO]
	//    - version of 1-4 that takes ROI
	//    - optional loop subd vtx placement for 1-4 split
	//    - 1-3 split, also with 1-3+flip  (eg sqrt(3) subd)
	//    - feedback for new tris/verts
	//    - reprojection (here?)

};
typedef MeshSubdivider<float> MeshSubdividerf;
typedef MeshSubdivider<double> MeshSubdividerd;


}

