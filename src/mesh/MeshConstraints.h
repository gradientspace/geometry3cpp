#pragma once

#include <g3types.h>

namespace g3
{

enum class EdgeRefineFlags
{
	NoConstraint = 0,
	NoFlip = 1,
	NoSplit = 2,
	NoCollapse = 4,
	FullyConstrained = NoFlip | NoSplit | NoCollapse,

	PreserveTopology = 8,	// this flag just means we want to avoid 'losing' this edge
							// when we collapse a neighbour (but allow the edge itself to be collapsed).
							// Eg when constraining edges for Reduce
};


struct EdgeConstraint
{
public:
    EdgeRefineFlags refineFlags;
	IProjectionTargetPtr Target;        // edge is associated with this projection Target.
                                            // Currently only used as information, we do not explicitly
                                            // project edges onto targets (must also set VertexConstraint)

    int TrackingSetID;               // not actually a constraint, but allows is to find descendents
                                            // of an constrained input edge


	EdgeConstraint()
	{
		refineFlags = EdgeRefineFlags::NoConstraint;
		Target = nullptr;
		TrackingSetID = -1;
	}

    EdgeConstraint(EdgeRefineFlags rflags)
    {
        refineFlags = rflags;
        Target = nullptr;
        TrackingSetID = -1;
    }

    EdgeConstraint(EdgeRefineFlags rflags, IProjectionTargetPtr target)
    {
        refineFlags = rflags;
        Target = target;
        TrackingSetID = -1;
    }

    bool CanFlip() {
        return ((int)refineFlags & (int)EdgeRefineFlags::NoFlip) == 0;
    }
    bool CanSplit() {
        return ((int)refineFlags & (int)EdgeRefineFlags::NoSplit) == 0;
    }
    bool CanCollapse() {
        return ((int)refineFlags & (int)EdgeRefineFlags::NoCollapse) == 0;
    }
    bool NoModifications() {
        return ((int)refineFlags & (int)EdgeRefineFlags::FullyConstrained) == (int)EdgeRefineFlags::FullyConstrained;
    }

    bool IsUnconstrained() {
        return refineFlags == EdgeRefineFlags::NoConstraint && Target == nullptr;
    }

	static EdgeConstraint Unconstrained() { return EdgeConstraint(EdgeRefineFlags::NoConstraint); }
    static EdgeConstraint NoFlips() { return EdgeConstraint(EdgeRefineFlags::NoFlip); }
    static EdgeConstraint FullyConstrained() { return EdgeConstraint(EdgeRefineFlags::FullyConstrained); }
};



struct VertexConstraint
{
public:
	static constexpr int InvalidSetID = -1;      // clients should interpret negative values as invalid
												// (in case you wanted to use negative values for something else...)

    bool Fixed;
    int FixedSetID;      // in Remesher, we can allow two Fixed vertices with 
                                // same FixedSetID to be collapsed together

	IProjectionTargetPtr Target;    // vertex is constrained to lie on this projection Target.
                                        // Fixed and Target are mutually exclusive
        
	VertexConstraint() {
		Fixed = false;
		FixedSetID = InvalidSetID;
		Target = nullptr;
	}

    VertexConstraint(bool isFixed, int setID = InvalidSetID)
    {
        Fixed = isFixed;
        FixedSetID = setID;
        Target = nullptr;
    }

    VertexConstraint(IProjectionTargetPtr target)
    {
        Fixed = false;
        FixedSetID = InvalidSetID;
        Target = target;
    }

	static VertexConstraint Unconstrained() { return VertexConstraint(false); }
	static VertexConstraint Pinned() { return VertexConstraint(true); }
};





class MeshConstraints
{
public:
	std::map<int, EdgeConstraint> Edges;

	int set_id_counter;         // use this to allocate FixedSetIDs

	MeshConstraints()
	{
		set_id_counter = 0;
	}

	int AllocateSetID() {
		return set_id_counter++;
	}

	bool HasEdgeConstraint(int eid)
	{
		return Edges.find(eid) != Edges.end();
	}

	EdgeConstraint GetEdgeConstraint(int eid)
	{
		auto found = Edges.find(eid);
		if (found == Edges.end())
			return EdgeConstraint::Unconstrained();
		else
			return found->second;
	}

	void SetOrUpdateEdgeConstraint(int eid, EdgeConstraint ec)
	{
		Edges[eid] = ec;
	}

	void ClearEdgeConstraint(int eid)
	{
		Edges.erase(eid);
	}

	void FindConstrainedEdgesBySetID(int setID, std::vector<int> & result)
	{
		for ( auto pair : Edges ) {
			if ( pair.second.TrackingSetID == setID )
				result.push_back(pair.first);
		}
	}


	std::map<int, VertexConstraint> Vertices;

	bool HasVertexConstraint(int vid)
	{
		return Vertices.find(vid) != Vertices.end();
	}

	VertexConstraint GetVertexConstraint(int vid)
	{
		auto found = Vertices.find(vid);
		if (found == Vertices.end())
			return VertexConstraint::Unconstrained();
		else
			return found->second;
	}

	bool GetVertexConstraint(int vid, VertexConstraint & vc)
	{
		auto found = Vertices.find(vid);
		if (found == Vertices.end()) 
			return false;
		vc = found->second;
		return true;
	}

	void SetOrUpdateVertexConstraint(int vid, VertexConstraint vc)
	{
		Vertices[vid] = vc;
	}

	void ClearVertexConstraint(int vid)
	{
		Vertices.erase(vid);
	}


	//System.Collections.IEnumerable VertexConstraintsItr() {
	//	foreach(KeyValuePair<int, VertexConstraint> v in Vertices)
	//		yield return v;
	//}


	bool HasConstraints() {
		return Edges.size() > 0 || Vertices.size() > 0;
	}

};




}