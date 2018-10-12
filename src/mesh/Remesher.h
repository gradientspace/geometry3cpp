#pragma once

#include <MeshRefinerBase.h>
#include <SpatialInterfaces.h>
#include <BasicProjectionTargets.h>
#include <MeshUtil.h>

namespace g3
{


class Remesher : public MeshRefinerBase
{
protected:
    IProjectionTargetPtr target = nullptr;

public:
	bool EnableFlips = true;
	bool EnableCollapses = true;
	bool EnableSplits = true;
	bool EnableSmoothing = true;

    bool PreventNormalFlips = false;

	double MinEdgeLength = 0.001f;
	double MaxEdgeLength = 0.1f;

	double SmoothSpeedT = 0.1f;

	enum class SmoothTypes {
		Uniform, Cotan, MeanValue
	};
	SmoothTypes SmoothType = SmoothTypes::Uniform;

    // this overrides default smoothing if provided
    std::function<Vector3d(const DMesh3 &, int, double)> CustomSmoothF;

        
    // Sometimes we need to have very granular control over what happens to
    // specific vertices. This function allows client to specify such behavior.
    // Somewhat redundant w/ VertexConstraints, but simpler to code.
	enum class VertexControl
	{
		AllowAll = 0,
		NoSmooth = 1,
		NoProject = 2,
		NoMovement = NoSmooth | NoProject
	};
	std::function<VertexControl(int)> VertexControlF = nullptr;


    // other options

    // [RMS] this is a debugging aid, will break to debugger if these edges are touched, in debug builds
    std::vector<int> DebugEdges;

    // if Target is set, we can project onto it in different ways
	enum class TargetProjectionMode
	{
		NoProjection,           // disable projection
		AfterRefinement,        // do all projection after the refine/smooth pass
		Inline                  // project after each vertex update. Better results but more
								// expensive because eg we might create a vertex with
								// split, then project, then smooth, then project again.
	};
    TargetProjectionMode ProjectionMode = TargetProjectionMode::AfterRefinement;

    // this just lets us write more concise tests below
    bool EnableInlineProjection() { return ProjectionMode == TargetProjectionMode::Inline; }

    // enable parallel projection. Only applied in AfterRefinement mode
    bool EnableParallelProjection = true;

    // Enable parallel smoothing. This will produce slightly different results
    // across runs because we smooth in-place and hence there will be order side-effects.
    bool EnableParallelSmooth = true;

    // if smoothing is done in-place, we don't need an extra buffer, but also 
    // there will some randomness introduced in results. Probably worse.
    bool EnableSmoothInPlace = false;

	Remesher(DMesh3Ptr m) : MeshRefinerBase(m)
    {
	}
protected:
    Remesher()        // for subclasses that extend our behavior
    {
    }
public:



    IProjectionTargetPtr ProjectionTarget() { return this->target; }
    void SetProjectionTarget(IProjectionTargetPtr target) { this->target = target; }
	void SetProjectionTarget(MeshProjectionTargetPtr target) { this->target = std::dynamic_pointer_cast<IProjectionTarget>(target); }

        
    /// <summary>
    /// Set min/max edge-lengths to sane values for given target edge length
    /// </summary>
    void SetTargetEdgeLength(double fLength)
    {
        // from Botsch paper
        //MinEdgeLength = fLength * (4.0/5.0);
        //MaxEdgeLength = fLength * (4.0/3.0);
        // much nicer!! makes sense as when we split, edges are both > min !
        MinEdgeLength = fLength * 0.66;
        MaxEdgeLength = fLength * 1.33;
    }


    bool ENABLE_PROFILING = false;


    // glboal mesh info that, if known, lets us avoid work in remesh
    bool MeshIsClosed = false;

    /// <summary>
    /// we can vastly speed things up if we precompute some invariants. 
    /// You need to re-run this if you are changing the mesh externally
    /// between remesh passes, otherwise you will get weird results.
    /// But you will probably still come out ahead, computation-time-wise
    /// </summary>
    virtual void Precompute()
    {
        // if we know mesh is closed, we can skip is-boundary checks, which makes
        // the flip-valence tests much faster!
        MeshIsClosed = true;
        for (int eid : mesh->EdgeIndices() ) {
            if (mesh->IsBoundaryEdge(eid)) {
                MeshIsClosed = false;
                break;
            }
        }
    }



    /// <summary>
    /// Number of edges that were modified in previous Remesh pass.
    /// If this number gets small relative to edge count, you have probably converged (ish)
    /// </summary>
    int ModifiedEdgesLastPass = 0;





    /// <summary>
    /// Linear edge-refinement pass, followed by smoothing and projection
    /// - Edges are processed in prime-modulo-order to break symmetry
    /// - smoothing is done in parallel if EnableParallelSmooth = true
    /// - Projection pass if ProjectionMode == AfterRefinement
    /// - number of modified edges returned in ModifiedEdgesLastPass
    /// </summary>
	virtual void BasicRemeshPass() {
        if (mesh->TriangleCount() == 0)    // badness if we don't catch this...
            return;

        begin_pass();

        // Iterate over all edges in the mesh at start of pass.
        // Some may be removed, so we skip those.
        // However, some old eid's may also be re-used, so we will touch
        // some new edges. Can't see how we could efficiently prevent this.
        //
        begin_ops();

        int cur_eid = start_edges();
        bool done = false;
        ModifiedEdgesLastPass = 0;
        do {
            if (mesh->IsEdge(cur_eid)) {
                ProcessResult result = ProcessEdge(cur_eid);
                if (result == ProcessResult::Ok_Collapsed || result == ProcessResult::Ok_Flipped || result == ProcessResult::Ok_Split)
                    ModifiedEdgesLastPass++;
            }
            if (Cancelled())        // expensive to check every iter?
                return;
            cur_eid = next_edge(cur_eid, done);
        } while (done == false);
        end_ops();

        if (Cancelled())
            return;

        begin_smooth();
        if (EnableSmoothing && SmoothSpeedT > 0) {
			if (EnableSmoothInPlace)
				//FullSmoothPass_InPlace(EnableParallelSmooth);
				gBreakToDebugger();
            else
                FullSmoothPass_Buffer(EnableParallelSmooth);
            DoDebugChecks();
        }
        end_smooth();

        if (Cancelled())
            return;

        begin_project();
        if (target != nullptr && ProjectionMode == TargetProjectionMode::AfterRefinement) {
            FullProjectionPass();
            DoDebugChecks();
        }
        end_project();

        if (Cancelled())
            return;

        end_pass();
	}



    // subclasses can override these to implement custom behavior...

    virtual void OnEdgeSplit(int edgeID, int va, int vb, const DMesh3::EdgeSplitInfo & splitInfo)
    {
        // this is for subclasses...
    }

    virtual void OnEdgeCollapse(int edgeID, int va, int vb, const DMesh3::EdgeCollapseInfo & collapseInfo)
    {
        // this is for subclasses...
    }


protected:


    // start_edges() and next_edge() control the iteration over edges that will be refined.
    // Default here is to iterate over entire mesh->
    // Subclasses can override these two functions to restrict the affected edges (eg EdgeLoopRemesher)


    // We are using a modulo-index loop to break symmetry/pathological conditions. 
    // For example in a highly tessellated minimal cylinder, if the top/bottom loops have
    // sequential edge IDs, and all edges are < min edge length, then we can easily end
    // up successively collapsing each tiny edge, and eroding away the entire mesh!
    // By using modulo-index loop we jump around and hence this is unlikely to happen.
    const int nPrime = 31337;     // any prime will do...
    int nMaxEdgeID;
    virtual int start_edges()
    {
        nMaxEdgeID = mesh->MaxEdgeID();
        return 0;
    }

    virtual int next_edge(int cur_eid, bool & bDone)
    {
        int new_eid = (cur_eid + nPrime) % nMaxEdgeID;
        bDone = (new_eid == 0);
        return new_eid;
    }





	enum class ProcessResult {
		Ok_Collapsed,
		Ok_Flipped,
		Ok_Split,
		Ignored_EdgeIsFine,
        Ignored_EdgeIsFullyConstrained,
		Failed_OpNotSuccessful,
		Failed_NotAnEdge
	};

	virtual ProcessResult ProcessEdge(int edgeID) 
	{
        RuntimeDebugCheck(edgeID);

        EdgeConstraint constraint =
            (constraints == nullptr) ? EdgeConstraint::Unconstrained() : constraints->GetEdgeConstraint(edgeID);
        if (constraint.NoModifications())
            return ProcessResult::Ignored_EdgeIsFullyConstrained;

		// look up verts and tris for this edge
		int a = 0, b = 0, t0 = 0, t1 = 0;
		if ( mesh->GetEdge(edgeID, a, b, t0, t1) == false )
			return ProcessResult::Failed_NotAnEdge;
		bool bIsBoundaryEdge = (t1 == InvalidID);

        // look up 'other' verts c (from t0) and d (from t1, if it exists)
        Index2i ov = mesh->GetEdgeOpposingV(edgeID);
        int c = ov[0], d = ov[1];

		Vector3d vA = mesh->GetVertex(a);
		Vector3d vB = mesh->GetVertex(b);
        double edge_len_sqr = (vA-vB).squaredNorm();

        begin_collapse();

        // check if we should collapse, and also find which vertex we should collapse to,
        // in cases where we have constraints/etc
        int collapse_to = -1;
        bool bCanCollapse = EnableCollapses
                            && constraint.CanCollapse()
                            && edge_len_sqr < MinEdgeLength*MinEdgeLength
                            && can_collapse_constraints(edgeID, a, b, c, d, t0, t1, collapse_to);

		// optimization: if edge cd exists, we cannot collapse or flip. look that up here?
		//  funcs will do it internally...
		//  (or maybe we can collapse if cd exists? edge-collapse doesn't check for it explicitly...)

		// if edge length is too short, we want to collapse it
		bool bTriedCollapse = false;
		if ( bCanCollapse ) {

            int iKeep = b, iCollapse = a;
            Vector3d vNewPos = (vA + vB) * 0.5;

            // if either vtx is fixed, collapse to that position
            if ( collapse_to == b ) {
                vNewPos = vB;
            } else if ( collapse_to == a ) {
                iKeep = a; iCollapse = b;
                vNewPos = vA;
            } else
                vNewPos = get_projected_collapse_position(iKeep, vNewPos);

            // if new position would flip normal of one of the existing triangles
            // either one-ring, don't allow it
            if (PreventNormalFlips) {
                if (collapse_creates_flip_or_invalid(a, b, vNewPos, t0, t1) || collapse_creates_flip_or_invalid(b, a, vNewPos, t0, t1)) {
                    goto abort_collapse;
                }
            }

            // TODO be smart about picking b (keep vtx). 
            //    - swap if one is bdry vtx, for example?
            // lots of cases where we cannot collapse, but we should just let
            // mesh sort that out, right?
            COUNT_COLLAPSES++;
			DMesh3::EdgeCollapseInfo collapseInfo;
			MeshResult result = mesh->CollapseEdge(iKeep, iCollapse, collapseInfo);
			if ( result == MeshResult::Ok ) {
				mesh->SetVertex(iKeep, vNewPos);
                if (constraints != nullptr) {
                    constraints->ClearEdgeConstraint(edgeID);
                    constraints->ClearEdgeConstraint(collapseInfo.eRemoved0);
                    if ( collapseInfo.eRemoved1 != InvalidID )
                        constraints->ClearEdgeConstraint(collapseInfo.eRemoved1);
                    constraints->ClearVertexConstraint(iCollapse);
                }
                OnEdgeCollapse(edgeID, iKeep, iCollapse, collapseInfo);
                DoDebugChecks();

				return ProcessResult::Ok_Collapsed;
			} else 
				bTriedCollapse = true;
		}
        abort_collapse:

        end_collapse();
        begin_flip();

		// if this is not a boundary edge, maybe we want to flip
		bool bTriedFlip = false;
		if ( EnableFlips && constraint.CanFlip() && bIsBoundaryEdge == false ) {

            // can we do this more efficiently somehow?
            bool a_is_boundary_vtx = (MeshIsClosed) ? false : (bIsBoundaryEdge || mesh->IsBoundaryVertex(a));
			bool b_is_boundary_vtx = (MeshIsClosed) ? false : (bIsBoundaryEdge || mesh->IsBoundaryVertex(b));
			bool c_is_boundary_vtx = (MeshIsClosed) ? false : mesh->IsBoundaryVertex(c);
			bool d_is_boundary_vtx = (MeshIsClosed) ? false :  mesh->IsBoundaryVertex(d);
			int valence_a = mesh->GetVtxEdgeCount(a), valence_b = mesh->GetVtxEdgeCount(b);
			int valence_c = mesh->GetVtxEdgeCount(c), valence_d = mesh->GetVtxEdgeCount(d);
			int valence_a_target = (a_is_boundary_vtx) ? valence_a : 6;
			int valence_b_target = (b_is_boundary_vtx) ? valence_b : 6;
			int valence_c_target = (c_is_boundary_vtx) ? valence_c : 6;
			int valence_d_target = (d_is_boundary_vtx) ? valence_d : 6;


			// if total valence error improves by flip, we want to do it
			int curr_err = abs(valence_a-valence_a_target) + abs(valence_b-valence_b_target)
				                + abs(valence_c-valence_c_target) + abs(valence_d-valence_d_target);
			int flip_err = abs((valence_a-1)-valence_a_target) + abs((valence_b-1)-valence_b_target)
				                + abs((valence_c+1)-valence_c_target) + abs((valence_d+1)-valence_d_target);

            bool bTryFlip = flip_err < curr_err;
            if (bTryFlip && PreventNormalFlips && flip_inverts_normals(a, b, c, d, t0))
                bTryFlip = false;

            if (bTryFlip) {
                DMesh3::EdgeFlipInfo flipInfo;
                COUNT_FLIPS++;
				MeshResult result = mesh->FlipEdge(edgeID, flipInfo);
				if ( result == MeshResult::Ok ) {
                    DoDebugChecks();
					return ProcessResult::Ok_Flipped;
				} else 
					bTriedFlip = true;

			}
		}

        end_flip();
        begin_split();

		// if edge length is too long, we want to split it
		bool bTriedSplit = false;
		if ( EnableSplits && constraint.CanSplit() && edge_len_sqr > MaxEdgeLength*MaxEdgeLength ) {

			DMesh3::EdgeSplitInfo splitInfo;
            COUNT_SPLITS++;
			MeshResult result = mesh->SplitEdge(edgeID, splitInfo);
			if ( result == MeshResult::Ok ) {
                update_after_split(edgeID, a, b, splitInfo);
                OnEdgeSplit(edgeID, a, b, splitInfo);
                DoDebugChecks();
				return ProcessResult::Ok_Split;
			} else
				bTriedSplit = true;
		}

        end_split();


		if ( bTriedFlip || bTriedSplit || bTriedCollapse )
			return ProcessResult::Failed_OpNotSuccessful;
		else
			return ProcessResult::Ignored_EdgeIsFine;
	}



    // After we split an edge, we have created a new edge and a new vertex.
    // The edge needs to inherit the constraint on the other pre-existing edge that we kept.
    // In addition, if the edge vertices were both constrained, then we /might/
    // want to also constrain this new vertex, possibly project to constraint target. 
    virtual void update_after_split(int edgeID, int va, int vb, const DMesh3::EdgeSplitInfo & splitInfo)
    {
        bool bPositionFixed = false;
        if (constraints != nullptr && constraints->HasEdgeConstraint(edgeID)) {
            // inherit edge constraint
            constraints->SetOrUpdateEdgeConstraint(splitInfo.eNewBN, constraints->GetEdgeConstraint(edgeID));

            // [RMS] update vertex constraints-> Note that there is some ambiguity here.
            //   Both verts being constrained doesn't inherently mean that the edge is on
            //   a constraint, that's why these checks are only applied if edge is constrained.
            //   But constrained edge doesn't necessarily mean we want to inherit vert constraints!!
            //
            //   although, pretty safe to assume that we would at least disable flips
            //   if both vertices are constrained to same line/curve. So, maybe this makes sense...
            //
            //   (perhaps edge constraint should be explicitly tagged to resolve this ambiguity??)

            // vert inherits Fixed if both orig edge verts Fixed, and both tagged with same SetID
            VertexConstraint ca = constraints->GetVertexConstraint(va);
            VertexConstraint cb = constraints->GetVertexConstraint(vb);
            if (ca.Fixed && cb.Fixed) {
                int nSetID = (ca.FixedSetID > 0 && ca.FixedSetID == cb.FixedSetID) ?
                    ca.FixedSetID : VertexConstraint::InvalidSetID;
                constraints->SetOrUpdateVertexConstraint(splitInfo.vNew,
                    new VertexConstraint(true, nSetID));
                bPositionFixed = true;
            }

            // vert inherits Target if:
            //  1) both source verts and edge have same Target, and is same as edge target
            //  2) either vert has same target as edge, and other vert is fixed
            if ( ca.Target != nullptr || cb.Target != nullptr) {
				IProjectionTargetPtr edge_target = constraints->GetEdgeConstraint(edgeID).Target;
                IProjectionTargetPtr set_target = nullptr;
                if (ca.Target == cb.Target && ca.Target == edge_target)
                    set_target = edge_target;
                else if (ca.Target == edge_target && cb.Fixed)
                    set_target = edge_target;
                else if (cb.Target == edge_target && ca.Fixed)
                    set_target = edge_target;
                if ( set_target != nullptr ) {
                    constraints->SetOrUpdateVertexConstraint(splitInfo.vNew,
                        new VertexConstraint(set_target));
                    project_vertex(splitInfo.vNew, set_target);
                    bPositionFixed = true;
                }
            }
        }

        if ( EnableInlineProjection() && bPositionFixed == false && target != nullptr ) {
            project_vertex(splitInfo.vNew, target);
        }
    }


    virtual void project_vertex(int vID, IProjectionTargetPtr targetIn)
    {
        Vector3d curpos = mesh->GetVertex(vID);
        Vector3d projected = targetIn->Project(curpos, vID);
        mesh->SetVertex(vID, projected);
    }

    // used by collapse-edge to get projected position for new vertex
    virtual Vector3d get_projected_collapse_position(int vid, const Vector3d & vNewPos)
    {
        if (constraints != nullptr) {
            VertexConstraint vc = constraints->GetVertexConstraint(vid);
            if (vc.Target != nullptr)
                return vc.Target->Project(vNewPos, vid);
            if (vc.Fixed)
                return vNewPos;
        }
        // no constraint applied, so if we have a target surface, project to that
        if ( EnableInlineProjection() && target != nullptr ) {
            if (VertexControlF == nullptr || ((int)VertexControlF(vid) & (int)VertexControl::NoProject) == 0)
                return target->Project(vNewPos, vid);
        }
        return vNewPos;
    }







	virtual void apply_to_smooth_vertices(const std::function<void(int)> & apply_f) {
		for (int vid : mesh->VertexIndices())
			apply_f(vid);
	}
	//virtual IEnumerable<int> smooth_vertices()
	//{
	//	return mesh->VertexIndices();
	//}



	//protected virtual void FullSmoothPass_InPlace(bool bParallel) {
 //       Func<DMesh3, int, double, Vector3d> smoothFunc = MeshUtil.UniformSmooth;
 //       if (CustomSmoothF != null) {
 //           smoothFunc = CustomSmoothF;
 //       } else {
 //           if (SmoothType == SmoothTypes.MeanValue)
 //               smoothFunc = MeshUtil.MeanValueSmooth;
 //           else if (SmoothType == SmoothTypes.Cotan)
 //               smoothFunc = MeshUtil.CotanSmooth;
 //       }

 //       Action<int> smooth = (vID) => {
 //           bool bModified = false;
 //           Vector3d vSmoothed = ComputeSmoothedVertexPos(vID, smoothFunc, out bModified);
 //           if ( bModified )
 //               mesh->SetVertex(vID, vSmoothed);
 //       };

 //       if (bParallel) {
 //           gParallel.ForEach<int>(smooth_vertices(), smooth);
 //       } else {
 //           foreach ( int vID in smooth_vertices() )
 //               smooth(vID);
 //       }
	//}




    virtual void FullSmoothPass_Buffer(bool bParallel)
    {
        InitializeVertexBufferForPass();

		std::function<Vector3d(const DMesh3 &, int, double)> smoothFunc = UniformSmooth;
        //Func<DMesh3, int, double, Vector3d> smoothFunc = MeshUtil.UniformSmooth;
        if (CustomSmoothF != nullptr) {
            smoothFunc = CustomSmoothF;
        } else {
			if (SmoothType == SmoothTypes::MeanValue || SmoothType == SmoothTypes::Cotan)
				gBreakToDebugger(); // not implemented!
            //if (SmoothType == SmoothTypes.MeanValue)
            //    smoothFunc = MeshUtil.MeanValueSmooth;
            //else if (SmoothType == SmoothTypes.Cotan)
            //    smoothFunc = MeshUtil.CotanSmooth;
        }

		auto smooth = [&](int vID) {
            bool bModified = false;
            Vector3d vSmoothed = ComputeSmoothedVertexPos(vID, smoothFunc, bModified);
            if (bModified) {
                vModifiedV[vID] = true;
                vBufferV[vID] = vSmoothed;
            }
        };

        //if (bParallel) {
        //    gParallel.ForEach<int>(smooth_vertices(), smooth);
        //} else {
        //    foreach (int vID in smooth_vertices())
        //        smooth(vID);
        //}
		apply_to_smooth_vertices(smooth);

        ApplyVertexBuffer(bParallel);
    }



	dvector<Vector3d> vBufferV;
	std::vector<bool> vModifiedV;

    virtual void InitializeVertexBufferForPass()
    {
        if (vBufferV.size() < mesh->MaxVertexID())
            vBufferV.resize(mesh->MaxVertexID() + mesh->MaxVertexID() / 5);
        if (vModifiedV.size() < mesh->MaxVertexID()) {
            vModifiedV.resize(2 * mesh->MaxVertexID());
        } 
		
		vModifiedV.assign(vModifiedV.size(), false);
    }

    virtual void ApplyVertexBuffer(bool bParallel)
    {
		for (int vid : mesh->VertexIndices()) {
			if (vModifiedV[vid])
				mesh->SetVertex(vid, vBufferV[vid]);
		}

        // [TODO] can probably use block-parallel here...
        //if (bParallel) {
        //    gParallel.BlockStartEnd(0, mesh->MaxVertexID-1, (a,b) => {
        //        for (int vid = a; vid <= b; vid++) {
        //            if (vModifiedV[vid])
        //                mesh->SetVertex(vid, vBufferV[vid]);
        //        }
        //    });
        //} else {
        //    foreach (int vid in mesh->VertexIndices()) {
        //        if (vModifiedV[vid])
        //            mesh->SetVertex(vid, vBufferV[vid]);
        //    }
        //}
    }




    /// <summary>
    /// This computes smoothed positions w/ proper constraints/etc.
    /// Does not modify mesh->
    /// </summary>
    virtual Vector3d ComputeSmoothedVertexPos(int vID, 
		std::function<Vector3d(const DMesh3 &, int, double)> smoothFunc, bool & bModified)
    {
        bModified = false;
        VertexConstraint vConstraint = VertexConstraint::Unconstrained();
        get_vertex_constraint(vID, vConstraint);
        if (vConstraint.Fixed)
            return mesh->GetVertex(vID);
        VertexControl vControl = (VertexControlF == nullptr) ? VertexControl::AllowAll : VertexControlF(vID);
        if ( ((int)vControl & (int)VertexControl::NoSmooth) != 0 )
            return mesh->GetVertex(vID);

        Vector3d vSmoothed = smoothFunc(*mesh, vID, SmoothSpeedT);
		gDevAssert(IsFinite(vSmoothed));     // this will really catch a lot of bugs...

        // project onto either vtx constraint target, or surface target
        if (vConstraint.Target != nullptr) {
            vSmoothed = vConstraint.Target->Project(vSmoothed, vID);
        } else if (EnableInlineProjection() && target != nullptr) {
            if ( ((int)vControl & (int)VertexControl::NoProject) == 0)
                vSmoothed = target->Project(vSmoothed, vID);
        }

        bModified = true;
        return vSmoothed;
    }





	virtual void apply_to_project_vertices(const std::function<void(int)> & apply_f) {
		for (int vid : mesh->VertexIndices())
			apply_f(vid);
	}

	// [RMS] how would we do something like this in C++?
	//virtual IEnumerable<int> project_vertices()
	//{
	//	return mesh->VertexIndices();
	//}



    // Project vertices onto projection target. 
    // We can do projection in parallel if we have .net 
    virtual void FullProjectionPass()
    {
        auto project = [&](int vID) {
            if (vertex_is_constrained(vID))
                return;
            if (VertexControlF != nullptr && ((int)VertexControlF(vID) & (int)VertexControl::NoProject) != 0)
                return;
            Vector3d curpos = mesh->GetVertex(vID);
            Vector3d projected = target->Project(curpos, vID);
            mesh->SetVertex(vID, projected);
        };
		
		apply_to_project_vertices(project);

		// [RMS] not sure how to do this...
        //if (EnableParallelProjection) {
        //    gParallel.ForEach<int>(project_vertices(), project);
        //} else {
        //    foreach (int vid in project_vertices())
        //        project(vid);
        //}
    }




	/*
	 * testing/debug/profiling stuff
	 */

protected:


    void RuntimeDebugCheck(int eid)
    {
//#ifdef DEBUG
		if (Contains(DebugEdges, eid))
			gBreakToDebugger();
//#endif
    }


    bool ENABLE_DEBUG_CHECKS = false;
    virtual void DoDebugChecks()
    {
        if (ENABLE_DEBUG_CHECKS == false)
            return;

        DebugCheckVertexConstraints();

        // [RMS] keeping this for now, is useful in testing that we are preserving group boundaries
        //foreach ( int eid in mesh->EdgeIndices() ) {
        //    if (mesh->IsGroupBoundaryEdge(eid))
        //        if (constraints->GetEdgeConstraint(eid).CanFlip) {
        //            Util.gBreakToDebugger();
        //            throw new Exception("fuck");
        //        }
        //}
        //foreach ( int vid in mesh->VertexIndices() ) {
        //    if (mesh->IsGroupBoundaryVertex(vid))
        //        if (constraints->GetVertexConstraint(vid).Target == null)
        //            Util.gBreakToDebugger();
        //}
    }

    void DebugCheckVertexConstraints()
    {
        if (constraints == nullptr)
            return;
		for ( std::pair<int,VertexConstraint> vc : constraints->Vertices ) {
			int vid = vc.first;
            if (vc.second.Target != nullptr) {
                Vector3d curpos = mesh->GetVertex(vid);
                Vector3d projected = vc.second.Target->Project(curpos, vid);
				if ((curpos - projected).squaredNorm() > 0.0001f)
					gBreakToDebugger();
            }
        }
    }




    //
    // profiling functions, turn on ENABLE_PROFILING to see output in console
    // 
    int COUNT_SPLITS, COUNT_COLLAPSES, COUNT_FLIPS;
    //Stopwatch AllOpsW, SmoothW, ProjectW, FlipW, SplitW, CollapseW;

    virtual void begin_pass() {
        if ( ENABLE_PROFILING ) {
            COUNT_SPLITS = COUNT_COLLAPSES = COUNT_FLIPS = 0;
            //AllOpsW = new Stopwatch();
            //SmoothW = new Stopwatch();
            //ProjectW = new Stopwatch();
            //FlipW = new Stopwatch();
            //SplitW = new Stopwatch();
            //CollapseW = new Stopwatch();
        }
    }

    virtual void end_pass() {
        if ( ENABLE_PROFILING ) {
            //System.Console.WriteLine(string.Format(
            //    "RemeshPass: T {0} V {1} splits {2} flips {3} collapses {4}", mesh->TriangleCount, mesh->VertexCount, COUNT_SPLITS, COUNT_FLIPS, COUNT_COLLAPSES
            //    ));
            //System.Console.WriteLine(string.Format(
            //    "           Timing1:  ops {0} smooth {1} project {2}", Util.ToSecMilli(AllOpsW.Elapsed), Util.ToSecMilli(SmoothW.Elapsed), Util.ToSecMilli(ProjectW.Elapsed)
            //    ));
            //System.Console.WriteLine(string.Format(
            //    "           Timing2:  collapse {0} flip {1} split {2}", Util.ToSecMilli(CollapseW.Elapsed), Util.ToSecMilli(FlipW.Elapsed), Util.ToSecMilli(SplitW.Elapsed)
            //    ));
        }
    }

    virtual void begin_ops() {
        //if ( ENABLE_PROFILING ) AllOpsW.Start();
    }
    virtual void end_ops() {
        //if ( ENABLE_PROFILING ) AllOpsW.Stop();
    }
    virtual void begin_smooth() {
        //if ( ENABLE_PROFILING ) SmoothW.Start();
    }
    virtual void end_smooth() {
        //if ( ENABLE_PROFILING ) SmoothW.Stop();
    }
    virtual void begin_project() {
        //if ( ENABLE_PROFILING ) ProjectW.Start();
    }
    virtual void end_project() {
        //if ( ENABLE_PROFILING ) ProjectW.Stop();
    }

    virtual void begin_collapse() {
        //if ( ENABLE_PROFILING ) CollapseW.Start();
    }
    virtual void end_collapse() {
        //if ( ENABLE_PROFILING ) CollapseW.Stop();
    }
    virtual void begin_flip() {
        //if ( ENABLE_PROFILING ) FlipW.Start();
    }
    virtual void end_flip() {
        //if ( ENABLE_PROFILING ) FlipW.Stop();
    }
    virtual void begin_split() {
        //if ( ENABLE_PROFILING ) SplitW.Start();
    }
    virtual void end_split() {
        //if ( ENABLE_PROFILING ) SplitW.Stop();
    }

};


}