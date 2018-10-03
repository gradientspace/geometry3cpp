#pragma once

#include <g3types.h>
#include <DMesh3.h>

namespace g3 
{

class SphereGenerator
{
public:
	Frame3d vFrame;
	double fRadius;
	int nSlices;
	int nStacks;

	SphereGenerator();
	virtual ~SphereGenerator() = default;

	virtual void Generate( DMesh3 * pMesh );
};


}

