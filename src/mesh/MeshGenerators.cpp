#include <geometry3PCH.h>
#include <MeshGenerators.h>

using namespace g3;

template<typename Real>
SphereGenerator<Real>::SphereGenerator()
{
	vFrame = Frame3<Real>();
	fRadius = (Real)1.0;
	nSlices = 8;
	nStacks = 8;
}

template<typename Real>
void SphereGenerator<Real>::Generate( IDynamicMesh<Real> * pMesh )
{
	GroupID gID = 0;

	pMesh->AppendVertex( vFrame.ToWorldCoords( Vector3<Real>(0,fRadius,0) ) );
	for (int j = 0; j < nStacks - 1; ++j) {
		Real polar = Wml::Math<Real>::PI * (Real)(j+1) / (Real)(nStacks);
		Real sp = std::sin(polar);
		Real cp = std::cos(polar);
		for (int i = 0; i < nSlices; ++i) {
			Real azimuth = (Real)2.0 * Wml::Math<Real>::PI * (Real)(i) / (Real)(nSlices);
			Real sa = std::sin(azimuth);
			Real ca = std::cos(azimuth);
			pMesh->AppendVertex( vFrame.ToWorldCoords(
				Vector3<Real>(fRadius*sp*ca,fRadius*cp,fRadius*sp*sa) ) );
		}
	}
	pMesh->AppendVertex( vFrame.ToWorldCoords( Vector3<Real>(0,-fRadius,0) ) );

	for (int i = 0; i < nSlices; ++i) {
		int a = i + 1;
		int b = (i + 1) % nSlices + 1;
		pMesh->AppendTriangle( { 0, b, a }, gID );
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
			pMesh->AppendTriangle( { a, a1, b1 }, gID );
			pMesh->AppendTriangle( { a, b1, b }, gID );
		}
	}

	int nVertices = pMesh->GetVertexCount();
	for (int i = 0; i < nSlices; ++i) {
		int a = i + nSlices * (nStacks - 2) + 1;
		int b = (i + 1) % nSlices + nSlices * (nStacks - 2) + 1;
		pMesh->AppendTriangle( { nVertices - 1, a, b }, gID );
	}

	// transfer to frame and compute normals
	// TODO
}



namespace g3
{
template class SphereGenerator<float>;
template class SphereGenerator<double>;
}
