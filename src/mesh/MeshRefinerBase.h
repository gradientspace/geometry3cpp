#pragma once

#include <g3types.h>
#include <ProgressCancel.h>
#include <DMesh3.h>
#include <index_util.h>
#include <VectorUtil.h>
#include <MeshConstraints.h>

namespace g3
{

class MeshRefinerBase
{
protected:
	DMesh3Ptr mesh;
    MeshConstraintsPtr constraints = nullptr;


public:

    // if true, then when two Fixed vertices have the same non-invalid SetID,
    // we treat them as not fixed and allow collapse
    bool AllowCollapseFixedVertsWithSameSetID = true;


    /// <summary>
    /// If normals dot product is less than this, we consider it a normal flip. default = 0
    /// </summary>
    double GetEdgeFlipTolerance() { return edge_flip_tol; }
	void GetEdgeFlipTolerance(double value) { edge_flip_tol = Wml::Mathd::Clamp(value, -1.0, 1.0); }
private:
	double edge_flip_tol = 0.0f;
public:


    MeshRefinerBase(DMesh3Ptr mesh) {
        this->mesh = mesh;
    }

    MeshRefinerBase() {
    }

	virtual ~MeshRefinerBase() {}

    DMesh3Ptr Mesh() { return mesh; }
	MeshConstraintsPtr Constraints() { return constraints; }
        
    //! This object will be modified !!!
    void SetExternalConstraints(MeshConstraintsPtr cons) {
        constraints = cons;
    }


    /// <summary>
    /// Set this to be able to cancel running remesher
    /// </summary>
	ProgressCancelPtr Progress = nullptr;

    /// <summary>
    /// if this returns true, abort computation. 
    /// </summary>
    virtual bool Cancelled() {
        return (Progress == nullptr) ? false : Progress->Cancelled();
    }


    double edge_flip_metric(const Vector3d & n0, const Vector3d & n1)
    {
        if (edge_flip_tol == 0) {
            return n0.dot(n1);
        } else {
            return n0.normalized().dot(n1.normalized());
        }
    }


    /// <summary>
    /// check if edge collapse will create a face-normal flip. 
    /// Also checks if collapse would violate link condition, since
    /// we are iterating over one-ring anyway.
    /// 
    /// This only checks one-ring of vid, so you have to call it twice,
    /// with vid and vother reversed, to check both one-rings
    /// </summary>
    bool collapse_creates_flip_or_invalid(int vid, int vother, Vector3d & newv, int tc, int td)
    {
        Vector3d va = Vector3d::Zero(), vb = Vector3d::Zero(), vc = Vector3d::Zero();
        for (int tid : mesh->VtxTrianglesItr(vid)) {
            if (tid == tc || tid == td)
                continue;
            Index3i curt = mesh->GetTriangle(tid);
            if (curt[0] == vother || curt[1] == vother || curt[2] == vother)
                return true;		// invalid nbrhood for collapse
            mesh->GetTriVertices(tid, va, vb, vc);
            Vector3d ncur = (vb - va).cross(vc - va);
            double sign = 0;
            if (curt[0] == vid) {
                Vector3d nnew = (vb - newv).cross(vc - newv);
                sign = edge_flip_metric(ncur, nnew);
            } else if (curt[1] == vid) {
                Vector3d nnew = (newv - va).cross(vc - va);
                sign = edge_flip_metric(ncur, nnew);
            } else if (curt[2] == vid) {
                Vector3d nnew = (vb - va).cross(newv - va);
                sign = edge_flip_metric(ncur, nnew);
            } else
                throw new std::exception("should never be here!");
            if (sign <= edge_flip_tol)
                return true;
        }
        return false;
    }



    /// <summary>
    /// Check if edge flip might reverse normal direction. 
    /// 
    /// Not entirely clear on how to implement this test. 
    /// Currently checking if any normal-pairs are reversed.
    /// </summary>
    bool flip_inverts_normals(int a, int b, int c, int d, int t0)
    {
        Vector3d vC = mesh->GetVertex(c), vD = mesh->GetVertex(d);
        Index3i tri_v = mesh->GetTriangle(t0);
        int oa = a, ob = b;
        orient_tri_edge(oa, ob, tri_v);
        Vector3d vOA = mesh->GetVertex(oa), vOB = mesh->GetVertex(ob);
        Vector3d n0 = FastNormalDirection(vOA, vOB, vC);
        Vector3d n1 = FastNormalDirection(vOB, vOA, vD);
        Vector3d f0 = FastNormalDirection(vC, vD, vOB);
        if ( edge_flip_metric(n0, f0) <= edge_flip_tol || edge_flip_metric(n1, f0) <= edge_flip_tol)
            return true;
        Vector3d f1 = FastNormalDirection(vD, vC, vOA);
        if (edge_flip_metric(n0, f1) <= edge_flip_tol || edge_flip_metric(n1, f1) <= edge_flip_tol)
            return true;

        // this only checks if output faces are pointing towards eachother, which seems 
        // to still result in normal-flips in some cases
        //if (f0.Dot(f1) < 0)
        //    return true;

        return false;
    }






    // Figure out if we can collapse edge eid=[a,b] under current constraint set.
    // First we resolve vertex constraints using can_collapse_vtx(). However this
    // does not catch some topological cases at the edge-constraint level, which 
    // which we will only be able to detect once we know if we are losing a or b.
    // See comments on can_collapse_vtx() for what collapse_to is for.
    bool can_collapse_constraints(int eid, int a, int b, int c, int d, int tc, int td, int & collapse_to)
    {
        collapse_to = -1;
        if (constraints == nullptr)
            return true;
        bool bVtx = can_collapse_vtx(eid, a, b, collapse_to);
        if (bVtx == false)
            return false;

        // when we lose a vtx in a collapse, we also lose two edges [iCollapse,c] and [iCollapse,d].
        // If either of those edges is constrained, we would lose that constraint.
        // This would be bad.
        int iCollapse = (collapse_to == a) ? b : a;
        if (c != InvalidID) {
            int ec = mesh->FindEdgeFromTri(iCollapse, c, tc);
            if (constraints->GetEdgeConstraint(ec).IsUnconstrained() == false)
                return false;
        }
        if (d != InvalidID) {
            int ed = mesh->FindEdgeFromTri(iCollapse, d, td);
            if (constraints->GetEdgeConstraint(ed).IsUnconstrained() == false)
                return false;
        }

        return true;
    }






    // resolve vertex constraints for collapsing edge eid=[a,b]. Generally we would
    // collapse a to b, and set the new position as 0.5*(v_a+v_b). However if a *or* b
    // are constrained, then we want to keep that vertex and collapse to its position.
    // This vertex (a or b) will be returned in collapse_to, which is -1 otherwise.
    // If a *and* b are constrained, then things are complicated (and documented below).
    bool can_collapse_vtx(int eid, int a, int b, int & collapse_to)
    {
        collapse_to = -1;
        if (constraints == nullptr)
            return true;
		VertexConstraint ca, cb;
		constraints->GetVertexConstraint(a, ca);
        constraints->GetVertexConstraint(b, cb);

        // no constraint at all
        if (ca.Fixed == false && cb.Fixed == false && ca.Target == nullptr && cb.Target == nullptr)
            return true;

        // handle a or b fixed
        if (ca.Fixed == true && cb.Fixed == false) {
            // if b is fixed to a target, and it is different than a's target, we can't collapse
            if (cb.Target != nullptr && cb.Target != ca.Target)
                return false;
            collapse_to = a;
            return true;
        }
        if (cb.Fixed == true && ca.Fixed == false) {
            if (ca.Target != nullptr && ca.Target != cb.Target)
                return false;
            collapse_to = b;
            return true;
        }
        // if both fixed, and options allow, treat this edge as unconstrained (eg collapse to midpoint)
        // [RMS] tried picking a or b here, but something weird happens, where
        //   eg cylinder cap will entirely erode away. Somehow edge lengths stay below threshold??
        if (AllowCollapseFixedVertsWithSameSetID
                && ca.FixedSetID >= 0
                && ca.FixedSetID == cb.FixedSetID) {
            return true;
        }

        // handle a or b w/ target
        if (ca.Target != nullptr && cb.Target == nullptr) {
            collapse_to = a;
            return true;
        }
        if (cb.Target != nullptr && ca.Target == nullptr) {
            collapse_to = b;
            return true;
        }
        // if both vertices are on the same target, and the edge is on that target,
        // then we can collapse to either and use the midpoint (which will be projected
        // to the target). *However*, if the edge is not on the same target, then we 
        // cannot collapse because we would be changing the constraint topology!
        if (cb.Target != nullptr && ca.Target != nullptr && ca.Target == cb.Target) {
            if (constraints->GetEdgeConstraint(eid).Target == ca.Target)
                return true;
        }

        return false;
    }




    bool vertex_is_fixed(int vid)
    {
        if (constraints != nullptr && constraints->GetVertexConstraint(vid).Fixed)
            return true;
        return false;
    }
    bool vertex_is_constrained(int vid)
    {
        if (constraints != nullptr) {
            VertexConstraint vc = constraints->GetVertexConstraint(vid);
            if (vc.Fixed || vc.Target != nullptr)
                return true;
        }
        return false;
    }

    VertexConstraint get_vertex_constraint(int vid)
    {
        if (constraints != nullptr)
            return constraints->GetVertexConstraint(vid);
        return VertexConstraint::Unconstrained();
    }
    bool get_vertex_constraint(int vid, VertexConstraint & vc)
    {
        return (constraints == nullptr) ? false :
            constraints->GetVertexConstraint(vid, vc);
    }

};

}