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

	SphereGenerator() {
		vFrame = Frame3d();
		fRadius = 1.0;
		nSlices = 8;
		nStacks = 8;
	}
	virtual ~SphereGenerator() = default;

	virtual void Generate(DMesh3 * pMesh) {
		GroupID gID = 0;

		pMesh->AppendVertex(vFrame.ToWorldCoords(Vector3d(0, fRadius, 0)));
		for (int j = 0; j < nStacks - 1; ++j) {
			double polar = Wml::Mathd::PI * (double)(j + 1) / (double)(nStacks);
			double sp = std::sin(polar);
			double cp = std::cos(polar);
			for (int i = 0; i < nSlices; ++i) {
				double azimuth = (double)2.0 * Wml::Mathd::PI * (double)(i) / (double)(nSlices);
				double sa = std::sin(azimuth);
				double ca = std::cos(azimuth);
				pMesh->AppendVertex(vFrame.ToWorldCoords(
					Vector3d(fRadius*sp*ca, fRadius*cp, fRadius*sp*sa)));
			}
		}
		pMesh->AppendVertex(vFrame.ToWorldCoords(Vector3d(0, -fRadius, 0)));

		for (int i = 0; i < nSlices; ++i) {
			int a = i + 1;
			int b = (i + 1) % nSlices + 1;
			pMesh->AppendTriangle({ 0, b, a }, gID);
		}

		for (int j = 0; j < nStacks - 2; ++j) {
			int aStart = j * nSlices + 1;
			int bStart = (j + 1) * nSlices + 1;
			for (int i = 0; i < nSlices; ++i)
			{
				const int a = aStart + i;
				const int a1 = aStart + (i + 1) % nSlices;
				const int b = bStart + i;
				const int b1 = bStart + (i + 1) % nSlices;
				pMesh->AppendTriangle({ a, a1, b1 }, gID);
				pMesh->AppendTriangle({ a, b1, b }, gID);
			}
		}

		int nVertices = pMesh->VertexCount();
		for (int i = 0; i < nSlices; ++i) {
			int a = i + nSlices * (nStacks - 2) + 1;
			int b = (i + 1) % nSlices + nSlices * (nStacks - 2) + 1;
			pMesh->AppendTriangle({ nVertices - 1, a, b }, gID);
		}

		// transfer to frame and compute normals
		// TODO
	}
};


}

