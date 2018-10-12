#pragma once

#include <DMesh3.h>
#include <MeshQueries.h>
#include <SpatialInterfaces.h>


namespace g3
{


class DMeshAABBTree3 : public IMeshSpatial
{
protected:
	DMesh3Ptr mesh;
	int mesh_timestamp;
	int TopDownLeafMaxTriCount = 4;

public:
	static constexpr double DOUBLE_MAX = std::numeric_limits<double>::max();

	/// <summary>
	/// If non-null, only triangle IDs that pass this filter (ie filter is true) are considered
	/// </summary>
	std::function<bool(int)> TriangleFilterF = nullptr;


	DMeshAABBTree3(DMesh3Ptr m, bool autoBuild = false)
	{
		mesh = m;
		if (autoBuild)
			Build();
	}


	void Build()
	{
		build_top_down(false);
		mesh_timestamp = mesh->ShapeTimestamp();
	}


	virtual bool SupportsNearestTriangle() override { return true; }


	/// <summary>
	/// Find the triangle closest to p, and distance to it, within distance fMaxDist, or return InvalidID
	/// Use MeshQueries.TriangleDistance() to get more information
	/// </summary>
	virtual int FindNearestTriangle(const Vector3d & p, double & fNearestDistSqr, double fMaxDist = DOUBLE_MAX) override
	{
		if (mesh_timestamp != mesh->ShapeTimestamp())
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
				if (TriangleFilterF != nullptr && TriangleFilterF(ti) == false)
					continue;
				double fTriDistSqr = MeshQueries::TriDistanceSqr(*mesh, ti, p);
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





	virtual bool SupportsTriangleRayIntersection() override { return false; }


	virtual int FindNearestHitTriangle(const Ray3d & ray, double fMaxDist = std::numeric_limits<double>::max()) override 
	{
		return InvalidID;
	}


	virtual bool SupportsPointContainment() override { return false; }
	
	virtual bool IsInside(const Vector3d & p) override
	{
		return false;
	}









	class TreeTraversal
	{
	public:
		// return false to terminate this branch
		// arguments are box and depth in tree
		std::function<bool(const AxisAlignedBox3d &, int)> NextBoxF = [](const AxisAlignedBox3d & box, int depth) { return true; };

		std::function<void(int)> NextTriangleF = [](int tID) {};
	};


	/// <summary>
	/// Hierarchically descend through the tree nodes, calling the TreeTrversal functions at each level
	/// </summary>
	virtual void DoTraversal(TreeTraversal * traversal)
	{
		if (mesh_timestamp != mesh->ShapeTimestamp())
			throw std::exception("DMeshAABBTree3.FindNearestTriangle: mesh has been modified since tree construction");

		tree_traversal(root_index, 0, traversal);
	}

	// traversal implementation. you can override to customize this if necessary.
	virtual void tree_traversal(int iBox, int depth, TreeTraversal * traversal)
	{
		int idx = box_to_index[iBox];

		if (idx < triangles_end) {
			// triange-list case, array is [N t1 t2 ... tN]
			int n = index_list[idx];
			for (int i = 1; i <= n; ++i) {
				int ti = index_list[idx + i];
				if (TriangleFilterF != nullptr && TriangleFilterF(ti) == false)
					continue;
				traversal->NextTriangleF(ti);
			}
		}
		else {
			int i0 = index_list[idx];
			if (i0 < 0) {
				// negative index means we only have one 'child' box to descend into
				i0 = (-i0) - 1;
				if (traversal->NextBoxF(get_box(i0), depth + 1))
					tree_traversal(i0, depth + 1, traversal);
			}
			else {
				// positive index, two sequential child box indices to descend into
				i0 = i0 - 1;
				if (traversal->NextBoxF(get_box(i0), depth + 1))
					tree_traversal(i0, depth + 1, traversal);
				int i1 = index_list[idx + 1] - 1;
				if (traversal->NextBoxF(get_box(i1), depth + 1))
					tree_traversal(i1, depth + 1, traversal);
			}
		}
	}





protected:

	//
	// Internals - data structures, construction, etc
	//

	AxisAlignedBox3d get_box(int iBox)
	{
		const Vector3d & c = box_centers[iBox];
		const Vector3d & e = box_extents[iBox];
		Vector3d Min = c - e, Max = c + e;
		return AxisAlignedBox3d(Min.data(), Max.data());
	}

	AxisAlignedBox3d get_box_eps(int iBox, double epsilon = Wml::Mathd::ZERO_TOLERANCE)
	{
		const Vector3d & c = box_centers[iBox];
		Vector3d e = box_extents[iBox]; 
		e[0] += epsilon; e[1] += epsilon; e[2] += epsilon;
		Vector3d Min = c - e, Max = c + e;
		return AxisAlignedBox3d(Min.data(), Max.data());
	}


	double box_distance_sqr(int iBox, const Vector3d & v)
	{
		const Vector3d & c = box_centers[iBox];
		const Vector3d & e = box_extents[iBox];

		// per-axis delta is max(abs(p-c) - e, 0)... ?
		double dist_sqr = ((v - c).cwiseAbs() - e).cwiseMax(0).squaredNorm();
		return dist_sqr;
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
		std::vector<int> triangles(mesh->TriangleCount());
		std::vector<Vector3d> centers(mesh->TriangleCount());
		for (int ti : mesh->TriangleIndices()) {
			Vector3d centroid = mesh->GetTriCentroid(ti);
			double d2 = centroid.squaredNorm();
			bool bInvalid = isnan(d2) || (isfinite(d2) == false);
			gDevAssert(bInvalid == false);
			if (bInvalid == false) {
				triangles[i] = ti;
				centers[i] = mesh->GetTriCentroid(ti);
				i++;
			} // otherwise skip this tri
		}

		boxes_set tris;
		boxes_set nodes;
		AxisAlignedBox3d rootBox;
		int rootnode = 
			//(bSorted) ? split_tri_set_sorted(triangles, centers, 0, mesh->TriangleCount, 0, TopDownLeafMaxTriCount, tris, nodes, out rootBox) :
			split_tri_set_midpoint(triangles, centers, 0, mesh->TriangleCount(), 0, TopDownLeafMaxTriCount, tris, nodes, rootBox);

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
				box.Contain(mesh->GetTriBounds(triangles[iStart + i]));
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




public:

	// 1) make sure we can reach every tri in mesh through tree (also demo of how to traverse tree...)
	// 2) make sure that triangles are contained in parent boxes
	void TestCoverage()
	{
		std::vector<int> tri_counts(mesh->MaxTriangleID());
		std::fill(tri_counts.begin(), tri_counts.end(), 0);

		std::vector<int> parent_indices(box_to_index.size());
		std::fill(parent_indices.begin(), parent_indices.end(), 0);

		test_coverage(tri_counts, parent_indices, root_index);

		for (int ti : mesh->TriangleIndices()) {
			if (tri_counts[ti] != 1)
				gBreakToDebugger();
		}
	}

private:

	// accumulate triangle counts and track each box-parent index. 
	// also checks that triangles are contained in boxes
	void test_coverage(std::vector<int> & tri_counts, std::vector<int> & parent_indices, int iBox)
	{
		int idx = box_to_index[iBox];

		debug_check_child_tris_in_box(iBox);

		if (idx < triangles_end) {
			// triange-list case, array is [N t1 t2 ... tN]
			int n = index_list[idx];
			AxisAlignedBox3d box = get_box_eps(iBox);
			for (int i = 1; i <= n; ++i) {
				int ti = index_list[idx + i];
				tri_counts[ti]++;

				Index3i tv = mesh->GetTriangle(ti);
				for (int j = 0; j < 3; ++j) {
					Vector3d v = mesh->GetVertex(tv[j]);
					if ( ! box.Contains(v) )
						gBreakToDebugger();
				}
			}

		}
		else {
			int i0 = index_list[idx];
			if (i0 < 0) {
				// negative index means we only have one 'child' box to descend into
				i0 = (-i0) - 1;
				parent_indices[i0] = iBox;
				test_coverage(tri_counts, parent_indices, i0);
			}
			else {
				// positive index, two sequential child box indices to descend into
				i0 = i0 - 1;
				parent_indices[i0] = iBox;
				test_coverage(tri_counts, parent_indices, i0);
				int i1 = index_list[idx + 1];
				i1 = i1 - 1;
				parent_indices[i1] = iBox;
				test_coverage(tri_counts, parent_indices, i1);
			}
		}
	}
	// do full tree traversal below iBox and make sure that all triangles are further
	// than box-distance-sqr
	void debug_check_child_tri_distances(int iBox, const Vector3d & p)
	{
		double fBoxDistSqr = box_distance_sqr(iBox, p);

		// [TODO]
		TreeTraversal * t = new TreeTraversal();
		t->NextTriangleF = [&](int tID) {
			double fTriDistSqr = MeshQueries::TriDistanceSqr(*mesh, tID, p);
			if (fTriDistSqr < fBoxDistSqr)
				if (fabs(fTriDistSqr - fBoxDistSqr) > Wml::Mathd::ZERO_TOLERANCE * 100)
					gBreakToDebugger();
		};
		tree_traversal(iBox, 0, t);
	}

	// do full tree traversal below iBox to make sure that all child triangles are contained
	void debug_check_child_tris_in_box(int iBox)
	{
		AxisAlignedBox3d box = get_box_eps(iBox);
		TreeTraversal * t = new TreeTraversal(); 
		t->NextTriangleF = [&](int tID) {
			Index3i tv = mesh->GetTriangle(tID);
			for (int j = 0; j < 3; ++j) {
				Vector3d v = mesh->GetVertex(tv[j]);
				if (box.Contains(v) == false)
					gBreakToDebugger();
			}
		};
		tree_traversal(iBox, 0, t);
	}


};

}