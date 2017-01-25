#pragma once

#include <DMesh3.h>
#include <parallel_util.h>
#include <dvector_util.h>

namespace g3 
{


template<typename Real>
class MeshContour
{
public:
	typedef DMesh3<Real, g3::dvector> MeshType;

	MeshContour() = default;
	virtual ~MeshContour() = default;

	struct ContourEdge {
		Vector2i ev;
		Vector2f ef;
		VertexID newv;
	};
	struct ContourData {
		dvector<ContourEdge> vSplitEdges;
		dvector<Vector2i> vZeroEdges;
		dvector<VertexID> vZeroVertices;
	};

	// f must implement: Real func_name(Vector3<Real> v)
	template<typename Func>
	void Insert( DMesh3<Real> & mesh, ContourData & d, const Func & f);

	// project new vertices in d.vSplitEdges onto zero-set of f
	template<typename Func>
	void Project( DMesh3<Real> & mesh, ContourData & d, unsigned int nRootFindIters, const Func & f);

};
typedef MeshContour<float> MeshContourf;
typedef MeshContour<double> MeshContourd;




template<typename Real>
template<typename Func>
void MeshContour<Real>::Insert( DMesh3<Real> & mesh, ContourData & d, const Func & f )
{
	// compute values
	int nMaxVID = mesh.GetMaxVertexID();
	std::vector<float> values(nMaxVID);
	parallel_fill(values, 
		[=](unsigned int vid)->Real { return f(mesh.GetVertex(vid)); },
		[=](unsigned int vid)->bool { return mesh.IsVertex(vid); } );

	for (auto eid : mesh.edges()) {
		Vector2i ev = mesh.GetEdgeV(eid);
		if ( ev[0] >= nMaxVID || ev[1] >= nMaxVID )
			continue;		// new edge, ignore 

		Real f0 = values[ev[0]];
		Real f1 = values[ev[1]];
			
		// If both signs are 0, this edge is on-contour
		// If one sign is 0, that vertex is on-contour
		Real fEpsilon = Wml::Math<Real>::EPSILON;
		int n0 = (fabs(f0) < fEpsilon) ? 1 : 0;
		int n1 = (fabs(f1) < fEpsilon) ? 1 : 0;
		if (n0+n1 > 0) {
			if (n0+n1 == 2)
				d.vZeroEdges.push_back(ev);
			else
				d.vZeroVertices.push_back( n0 ? ev[0] : ev[1] );
			continue;
		}

		// no crossing
		if ( f0 * f1 > (Real)0 )
			continue;

		// found a zero-crossing, need to split edge and solve for zero
		EdgeSplitInfo splitInfo;
		auto eResult = mesh.SplitEdge(eid, splitInfo);
		if (eResult != MeshResult::Ok)
			continue;

		// SplitEdge just bisects edge - use LERP to do better
		Real t = f0 / (f0-f1);
		Vector3<Real> newPos = ((Real)1 - t)*mesh.GetVertex( ev[0] ) + (t)*mesh.GetVertex( ev[1] );
		mesh.SetVertex( splitInfo.vNew, newPos );

		ContourEdge c = { ev, Vector2f(f0,f1), splitInfo.vNew };
		d.vSplitEdges.push_back(c);
	}

}



template<typename Real>
template<typename Func>
void MeshContour<Real>::Project( DMesh3<Real> & mesh, ContourData & d, unsigned int nRootFindIters, const Func & f )
{
	parallel_apply(d.vSplitEdges, 
					[=, &mesh]( ContourEdge e ) {

		// set up interval
		int iLow = (e.ef[0] < e.ef[1]) ? 0 : 1;
		Vector3<Real> v_neg = mesh.GetVertex(e.ev[iLow]), v_pos = mesh.GetVertex(e.ev[(iLow+1)%2]);
		Real f_neg = e.ef[iLow], f_pos = e.ef[(iLow+1)%2];
		Vector3<Real> v_mid = mesh.GetVertex(e.newv);

		// solve steps
		unsigned int i = 0;
		do {
			Real f_mid = f(v_mid);
			if ( fabs(f_mid) < Wml::Math<Real>::EPSILON )
				break;		// best we can hope for!

			if (f_mid < 0) {
				v_neg = v_mid;
				f_neg = f_mid;
			} else {
				v_pos = v_mid;
				f_pos = f_mid;
			}

			Real t = f_neg / (f_neg-f_pos);
			v_mid = ((Real)1 - t) * v_neg + (t) * v_pos;
		} while ( ++i < nRootFindIters );

		// update midpoint vertex
		mesh.SetVertex(e.newv, v_mid);
	} );
}

}

