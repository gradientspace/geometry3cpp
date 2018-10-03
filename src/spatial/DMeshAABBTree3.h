#pragma once

#include <DMesh3.h>


namespace g3
{

class DMeshAABBTree3
{
protected:
	const DMesh3 & mesh;
	int mesh_timestamp;
	int TopDownLeafMaxTriCount = 4;

public:
	static constexpr double DOUBLE_MAX = std::numeric_limits<double>::max();

	DMeshAABBTree3(const DMesh3 & m, bool autoBuild = false) : mesh(m)
	{
		if (autoBuild)
			Build();
	}


	void Build()
	{
		build_top_down(false);
		mesh_timestamp = mesh.ShapeTimestamp();
	}




	/// <summary>
	/// Find the triangle closest to p, and distance to it, within distance fMaxDist, or return InvalidID
	/// Use MeshQueries.TriangleDistance() to get more information
	/// </summary>
	virtual int FindNearestTriangle(const Vector3d & p, double & fNearestDistSqr, double fMaxDist = DOUBLE_MAX)
	{
		if (mesh_timestamp != mesh.ShapeTimestamp())
			throw std::exception("DMeshAABBTree3.FindNearestTriangle: mesh has been modified since tree construction");

		fNearestDistSqr = (fMaxDist < DOUBLE_MAX) ? fMaxDist * fMaxDist : DOUBLE_MAX;
		int tNearID = InvalidID;
		find_nearest_tri(root_index, p, fNearestDistSqr, tNearID);
		return tNearID;
	}
	void find_nearest_tri(int iBox, const Vector3d & p, double & fNearestSqr, int & tID)
	{
		int idx = box_to_index[iBox];
		if (idx < triangles_end) {            // triange-list case, array is [N t1 t2 ... tN]
			int num_tris = index_list[idx];
			for (int i = 1; i <= num_tris; ++i) {
				int ti = index_list[idx + i];
				//if (TriangleFilterF != null && TriangleFilterF(ti) == false)
				//	continue;
				double fTriDistSqr = tri_distance_sqr(ti, p);
				if (fTriDistSqr < fNearestSqr) {
					fNearestSqr = fTriDistSqr;
					tID = ti;
				}
			}

		} else {                                // internal node, either 1 or 2 child boxes
			int iChild1 = index_list[idx];
			if (iChild1 < 0) {                 // 1 child, descend if nearer than cur min-dist
				iChild1 = (-iChild1) - 1;
				double fChild1DistSqr = box_distance_sqr(iChild1, p);
				if (fChild1DistSqr <= fNearestSqr)
					find_nearest_tri(iChild1, p, fNearestSqr, tID);

			} else {                            // 2 children, descend closest first
				iChild1 = iChild1 - 1;
				int iChild2 = index_list[idx + 1] - 1;

				double fChild1DistSqr = box_distance_sqr(iChild1, p);
				double fChild2DistSqr = box_distance_sqr(iChild2, p);
				if (fChild1DistSqr < fChild2DistSqr) {
					if (fChild1DistSqr < fNearestSqr) {
						find_nearest_tri(iChild1, p, fNearestSqr, tID);
						if (fChild2DistSqr < fNearestSqr)
							find_nearest_tri(iChild2, p, fNearestSqr, tID);
					}
				} else {
					if (fChild2DistSqr < fNearestSqr) {
						find_nearest_tri(iChild2, p, fNearestSqr, tID);
						if (fChild1DistSqr < fNearestSqr)
							find_nearest_tri(iChild1, p, fNearestSqr, tID);
					}
				}

			}
		}
	}









protected:

	//
	// Internals - data structures, construction, etc
	//

	double box_distance_sqr(int iBox, const Vector3d & p)
	{
		const Vector3d & c = box_centers[iBox];
		const Vector3d & e = box_extents[iBox];

		// per-axis delta is max(abs(p-c) - e, 0)... ?
		return ((p - c).cwiseAbs() - e).cwiseMax(0).squaredNorm();
	}

	double tri_distance_sqr(int ti, const Vector3d & point)
	{
		Vector3d V0, V1, V2;
		mesh.GetTriVertices(ti, V0, V1, V2);

		Vector3d diff = V0 - point;
		Vector3d edge0 = V1 - V0;
		Vector3d edge1 = V2 - V0;
		double a00 = edge0.squaredNorm();
		double a01 = edge0.dot(edge1);
		double a11 = edge1.squaredNorm();
		double b0 = diff.dot(edge0);
		double b1 = diff.dot(edge1);
		double c = diff.squaredNorm();
		double det = fabs(a00 * a11 - a01 * a01);
		double s = a01 * b1 - a11 * b0;
		double t = a01 * b0 - a00 * b1;
		double sqrDistance;

		if (s + t <= det) {
			if (s < 0) {
				if (t < 0) { // region 4
					if (b0 < 0) {
						t = 0;
						if (-b0 >= a00) {
							s = 1;
							sqrDistance = a00 + (2) * b0 + c;
						} else {
							s = -b0 / a00;
							sqrDistance = b0 * s + c;
						}
					} else {
						s = 0;
						if (b1 >= 0) {
							t = 0;
							sqrDistance = c;
						} else if (-b1 >= a11) {
							t = 1;
							sqrDistance = a11 + (2) * b1 + c;
						} else {
							t = -b1 / a11;
							sqrDistance = b1 * t + c;
						}
					}
				} else { // region 3
					s = 0;
					if (b1 >= 0) {
						t = 0;
						sqrDistance = c;
					} else if (-b1 >= a11) {
						t = 1;
						sqrDistance = a11 + (2) * b1 + c;
					} else {
						t = -b1 / a11;
						sqrDistance = b1 * t + c;
					}
				}
			} else if (t < 0) { // region 5
				t = 0;
				if (b0 >= 0) {
					s = 0;
					sqrDistance = c;
				} else if (-b0 >= a00) {
					s = 1;
					sqrDistance = a00 + (2) * b0 + c;
				} else {
					s = -b0 / a00;
					sqrDistance = b0 * s + c;
				}
			} else { // region 0
				   // minimum at interior point
				double invDet = (1) / det;
				s *= invDet; 
				t *= invDet;
				sqrDistance = s * (a00 * s + a01 * t + (2) * b0) +
					t * (a01 * s + a11 * t + (2) * b1) + c;
			}
		} else {
			double tmp0, tmp1, numer, denom;
			if (s < 0) { // region 2
				tmp0 = a01 + b0;
				tmp1 = a11 + b1;
				if (tmp1 > tmp0) {
					numer = tmp1 - tmp0;
					denom = a00 - (2) * a01 + a11;
					if (numer >= denom) {
						s = 1;
						t = 0;
						sqrDistance = a00 + (2) * b0 + c;
					} else {
						s = numer / denom;
						t = 1 - s;
						sqrDistance = s * (a00 * s + a01 * t + (2) * b0) +
							t * (a01 * s + a11 * t + (2) * b1) + c;
					}
				} else {
					s = 0;
					if (tmp1 <= 0) {
						t = 1;
						sqrDistance = a11 + (2) * b1 + c;
					} else if (b1 >= 0) {
						t = 0;
						sqrDistance = c;
					} else {
						t = -b1 / a11;
						sqrDistance = b1 * t + c;
					}
				}
			} else if (t < 0) {  // region 6
				tmp0 = a01 + b1;
				tmp1 = a00 + b0;
				if (tmp1 > tmp0) {
					numer = tmp1 - tmp0;
					denom = a00 - (2) * a01 + a11;
					if (numer >= denom) {
						t = 1;
						s = 0;
						sqrDistance = a11 + (2) * b1 + c;
					} else {
						t = numer / denom;
						s = 1 - t;
						sqrDistance = s * (a00 * s + a01 * t + (2) * b0) +
							t * (a01 * s + a11 * t + (2) * b1) + c;
					}
				} else {
					t = 0;
					if (tmp1 <= 0) {
						s = 1;
						sqrDistance = a00 + (2) * b0 + c;
					} else if (b0 >= 0) {
						s = 0;
						sqrDistance = c;
					} else {
						s = -b0 / a00;
						sqrDistance = b0 * s + c;
					}
				}
			} else {  // region 1
				numer = a11 + b1 - a01 - b0;
				if (numer <= 0) {
					s = 0;
					t = 1;
					sqrDistance = a11 + (2) * b1 + c;
				} else {
					denom = a00 - (2) * a01 + a11;
					if (numer >= denom) {
						s = 1;
						t = 0;
						sqrDistance = a00 + (2) * b0 + c;
					} else {
						s = numer / denom;
						t = 1 - s;
						sqrDistance = s * (a00 * s + a01 * t + (2) * b0) +
							t * (a01 * s + a11 * t + (2) * b1) + c;
					}
				}
			}
		}

		if (sqrDistance < 0)
			sqrDistance = 0;
		return sqrDistance;
	}









	// storage for box nodes. 
	//   - box_to_index is a pointer into index_list
	//   - box_centers and box_extents are the centers/extents of the bounding boxes
	dvector<int> box_to_index;
	dvector<Vector3d> box_centers;
	dvector<Vector3d> box_extents;

	// list of indices for a given box. There is *no* marker/sentinel between
	// boxes, you have to get the starting index from box_to_index[]
	//
	// There are three kinds of records:
	//   - if i < triangles_end, then the list is a number of triangles,
	//       stored as [N t1 t2 t3 ... tN]
	//   - if i > triangles_end and index_list[i] < 0, this is a single-child
	//       internal box, with index (-index_list[i])-1     (shift-by-one in case actual value is 0!)
	//   - if i > triangles_end and index_list[i] > 0, this is a two-child
	//       internal box, with indices index_list[i]-1 and index_list[i+1]-1
	dvector<int> index_list;

	// index_list[i] for i < triangles_end is a triangle-index list, otherwise box-index pair/single
	int triangles_end = -1;

	// box_to_index[root_index] is the root node of the tree
	int root_index = -1;




	struct boxes_set
	{
		dvector<int> box_to_index;
		dvector<Vector3d> box_centers;
		dvector<Vector3d> box_extents;
		dvector<int> index_list;
		int iBoxCur;
		int iIndicesCur;
		boxes_set() {
			iBoxCur = 0; iIndicesCur = 0;
		}
	};



	void build_top_down(bool bSorted)
	{
		// build list of valid triangles & centers. We skip any
		// triangles that have infinite/garbage vertices...
		int i = 0;
		std::vector<int> triangles(mesh.TriangleCount());
		std::vector<Vector3d> centers(mesh.TriangleCount());
		for (int ti : mesh.TriangleIndices()) {
			Vector3d centroid = mesh.GetTriCentroid(ti);
			double d2 = centroid.squaredNorm();
			bool bInvalid = isnan(d2) || (isfinite(d2) == false);
			gDevAssert(bInvalid == false);
			if (bInvalid == false) {
				triangles[i] = ti;
				centers[i] = mesh.GetTriCentroid(ti);
				i++;
			} // otherwise skip this tri
		}

		boxes_set tris;
		boxes_set nodes;
		AxisAlignedBox3d rootBox;
		int rootnode = 
			//(bSorted) ? split_tri_set_sorted(triangles, centers, 0, mesh.TriangleCount, 0, TopDownLeafMaxTriCount, tris, nodes, out rootBox) :
			split_tri_set_midpoint(triangles, centers, 0, mesh.TriangleCount(), 0, TopDownLeafMaxTriCount, tris, nodes, rootBox);

		box_to_index = tris.box_to_index;
		box_centers = tris.box_centers;
		box_extents = tris.box_extents;
		index_list = tris.index_list;
		triangles_end = tris.iIndicesCur;
		int iIndexShift = triangles_end;
		int iBoxShift = tris.iBoxCur;

		// ok now append internal node boxes & index ptrs
		for (i = 0; i < nodes.iBoxCur; ++i) {
			box_centers.insertAt(nodes.box_centers[i], iBoxShift + i);
			box_extents.insertAt(nodes.box_extents[i], iBoxShift + i);
			// internal node indices are shifted
			box_to_index.insertAt(iIndexShift + nodes.box_to_index[i], iBoxShift + i);
		}

		// now append index list
		for (i = 0; i < nodes.iIndicesCur; ++i) {
			int child_box = nodes.index_list[i];
			if (child_box < 0) {        // this is a triangles box
				child_box = (-child_box) - 1;
			} else {
				child_box += iBoxShift;
			}
			child_box = child_box + 1;
			index_list.insertAt(child_box, iIndexShift + i);
		}

		root_index = rootnode + iBoxShift;
	}





	int split_tri_set_midpoint( 
		std::vector<int> & triangles, 
		std::vector<Vector3d> & centers, 
		int iStart, int iCount, int depth, int minTriCount,
		boxes_set & tris, boxes_set & nodes, AxisAlignedBox3d & box)
	{
		box = AxisAlignedBox3d::EMPTY;
		int iBox = -1;

		if (iCount < minTriCount) {
			// append new triangles box
			iBox = tris.iBoxCur++;
			tris.box_to_index.insertAt(tris.iIndicesCur, iBox);

			tris.index_list.insertAt(iCount, tris.iIndicesCur++);
			for (int i = 0; i < iCount; ++i) {
				tris.index_list.insertAt(triangles[iStart + i], tris.iIndicesCur++);
				box.Contain(mesh.GetTriBounds(triangles[iStart + i]));
			}

			tris.box_centers.insertAt(box.Center(), iBox);
			tris.box_extents.insertAt(box.Extents(), iBox);

			return -(iBox + 1);
		}

		//compute interval along an axis and find midpoint
		int axis = depth % 3;
		Wml::Interval1d interval = Interval1d::EMPTY;
		for (int i = 0; i < iCount; ++i)
			interval.Contain(centers[iStart + i][axis]);
		double midpoint = interval.Center();

		int n0, n1;
		if ( interval.Length() > Wml::Mathd::ZERO_TOLERANCE) {
			// we have to re-sort the centers & triangles lists so that centers < midpoint
			// are first, so that we can recurse on the two subsets. We walk in from each side,
			// until we find two out-of-order locations, then we swap them.
			int l = 0;
			int r = iCount - 1;
			while (l < r) {
				// [RMS] is <= right here? if v.axis == midpoint, then this loop
				//   can get stuck unless one of these has an equality test. But
				//   I did not think enough about if this is the right thing to do...
				while (centers[iStart + l][axis] <= midpoint)
					l++;
				while (centers[iStart + r][axis] > midpoint)
					r--;
				if (l >= r)
					break;      //done!
								//swap
				Vector3d tmpc = centers[iStart + l]; centers[iStart + l] = centers[iStart + r];  centers[iStart + r] = tmpc;
				int tmpt = triangles[iStart + l]; triangles[iStart + l] = triangles[iStart + r]; triangles[iStart + r] = tmpt;
			}

			n0 = l;
			n1 = iCount - n0;
			gDevAssert(n0 >= 1 && n1 >= 1);
		} else {
			// interval is near-empty, so no point trying to do sorting, just split half and half
			n0 = iCount / 2;
			n1 = iCount - n0;
		}

		// create child boxes
		AxisAlignedBox3d box1;
		int child0 = split_tri_set_midpoint(triangles, centers, iStart, n0, depth + 1, minTriCount, tris, nodes, box);
		int child1 = split_tri_set_midpoint(triangles, centers, iStart + n0, n1, depth + 1, minTriCount, tris, nodes, box1);
		box.Contain(box1);

		// append new box
		iBox = nodes.iBoxCur++;
		nodes.box_to_index.insertAt(nodes.iIndicesCur, iBox);

		nodes.index_list.insertAt(child0, nodes.iIndicesCur++);
		nodes.index_list.insertAt(child1, nodes.iIndicesCur++);

		nodes.box_centers.insertAt(box.Center(), iBox);
		nodes.box_extents.insertAt(box.Extents(), iBox);

		return iBox;
	}


};







}