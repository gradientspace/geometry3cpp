#pragma once

#include <g3types.h>
#include <GeometryInterfaces.h>

namespace g3 
{

template<typename Real>
class SphereGenerator
{
public:
	Frame3<Real> vFrame;
	Real fRadius;
	int nSlices;
	int nStacks;

	SphereGenerator();
	virtual ~SphereGenerator() = default;

	virtual void Generate( IDynamicMesh<Real> * pMesh );
};
typedef SphereGenerator<float> SphereGeneratorf;
typedef SphereGenerator<double> SphereGeneratord;


}

