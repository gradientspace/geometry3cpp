#pragma once

#include <g3types.h>


/*
* Frame3 represents an orthonormal 3D coordinate frame, ie 3 orthogonal axes and an origin point.
* Although it is possible to pack all that information into a Matrix4, or a Matrix3+Vector3, it
* is not very convenient to work with. Making the axes available explicitly allows for much
* clearer code. 
*
* Frame3 also provides functions to:
*    - manipulate the frame by (eg) setting one axis to point in a particular direction
*    - transform vectors and frames to/from frame-relative coordinates
*    - (TODO: ray-intersection, project to plane, ...)
*
* Internally, Frame3 stores the inverse of the rotation matrix that would take the XYZ axes
* to the Frame axes. This is so that the axes are the rows of the matrix, and hence in
* row-major storage the 3 floats of each axis are contiguous. 
*
*/
namespace g3
{

template<class Real>
class Frame3
{
public:
	//! origin point of reference frame
	Vector3<Real> Origin;

	//! transpose of matrix that takes canonical axes to frame axes. Rows are frame axes.
	Matrix3<Real> InverseRotation;

public:
	//! create a reference frame at a specific origin
	Frame3( const Vector3<Real> & Origin = Vector3<Real>::Zero() )
		: Origin(Origin), InverseRotation(Matrix3<Real>::Identity()) {}

	//! create an orthogonal frame at an origin with a given z axis vector
	Frame3( const Vector3<Real> & vOrigin, const Vector3<Real> & vNormalizedAxis, int nAxis = 2, bool bUseFastPerpVectors = true );

	//! create frame at an origin with a given x/y/z vectors
	Frame3( const Vector3<Real> & vOrigin, const Vector3<Real> & vXAxis, const Vector3<Real> & vYAxis, const Vector3<Real> & vZAxis );

	//! create frame with a given rotation
	Frame3(const Matrix3<Real> & mRotation);

	//! create frame at an origin with a given rotation
	Frame3(const Vector3<Real> & vOrigin, const Matrix3<Real> & mRotation );

	//! copy a reference frame
	Frame3( const Frame3 & copy )
		: Origin(copy.Origin), InverseRotation(copy.InverseRotation) {}

	~Frame3();


	//! set the reference frame axes
	void SetFrame(const Vector3<Real> & vAxisX, const Vector3<Real> & vAxisY, const Vector3<Real> & vAxisZ);

	//! allow external setting of matrix...
	void SetFrame(const Matrix3<Real> & matRotate);


	//! get an axis of the frame
	Vector3<Real> Axis(unsigned int nAxis) const {
		return InverseRotation.row(nAxis);
	}

	//! access the axes of the frame
	Vector3<Real> X() const { 
		return InverseRotation.row( 0 ); }
	Vector3<Real> Y() const { 
		return InverseRotation.row( 1 ); }
	Vector3<Real> Z() const { 
		return InverseRotation.row( 2 ); }

	//! matrix that rotates canonical/unit XYZ axes to frame axes
	Matrix3<Real> GetRotation() const {
		return InverseRotation.transpose();
	}



	//! treat v as begin in frame-relative coords, transform to absolute/world coords
	void SetToWorldCoords( Vector3<Real> & v ) const { 
		v = Origin + InverseRotation.transpose() * v;
	}
	//! treat v as begin in frame-relative coords, transform to absolute/world coords
	Vector3<Real> ToWorldCoords(const Vector3<Real> & v) const {
		return Origin + InverseRotation.transpose() * v;
	}


	//! treat v as begin in absolute/world coords, transform to frame-relative coords
	void SetToAxisCoords( Vector3<Real> & v ) const {
		v = InverseRotation * (v - Origin);
	}
	//! treat v as begin in absolute/world coords, transform to frame-relative coords
	Vector3<Real> ToAxisCoords(const Vector3<Real> & v ) const {
		return InverseRotation * (v - Origin);
	}


	//! treat f as being relative to this frame, transform to absolute/world frame
	Frame3<Real> ToWorldCoords( const Frame3<Real> & f, bool bPreserveOrigin = false ) const {	
		Vector3<Real> x = InverseRotation.transpose() * f.X();
		Vector3<Real> y = InverseRotation.transpose() * f.Y();
		Vector3<Real> z = InverseRotation.transpose() * f.Z();
		Frame3<Real> l( (bPreserveOrigin) ? f.Origin : ToWorldCoords(f.Origin) );
		l.SetFrame(x,y,z);
		return l;
	}

	//! treat f as being an absolute/world frame, transform to be relative to this frame
	Frame3<Real> ToAxisCoords( const Frame3<Real> & f, bool bPreserveOrigin = false) const {
		Vector3<Real> x = InverseRotation * f.X();
		Vector3<Real> y = InverseRotation * f.Y();
		Vector3<Real> z = InverseRotation * f.Z();
		Frame3<Real> w( (bPreserveOrigin) ? f.Origin : ToAxisCoords(f.Origin) );
		w.SetFrame(x,y,z);
		return w;
	}



	//! translate the reference frame origin
	void Translate(const Vector3<Real> & vTranslate, bool bRelative = true);

	//! rotate the reference frame
	void Rotate(const Matrix3<Real> & mRotation, bool bReNormalize = true);

	//! rotate selected axis of this frame into toAxis
	void AlignAxis(int nAxis, const Vector3<Real> & toAxis, bool bNormalized = false);


	//! compute matrix that rotates this frame into toFrame
	void ComputeAlignmentMatrix( const Frame3<Real> & toFrame, Matrix3<Real> & matRotate );

	//! make axes perpendicular. can "preserve" one axis by setting nPreserveAxis (0=X,1=Y,2=Z)
	void ReNormalize(int nPreserveAxis = -1);

	//! intersect ray with plane defined by one of the axis vectors
	bool IntersectRay( const Vector3<Real> & vRayOrigin, const Vector3<Real> & vRayDirection,
					   Vector3<Real> & vRayHit, int nPlaneNormalAxis = 2 );
	Vector3<Real> IntersectRay( const Vector3<Real> & vRayOrigin, const Vector3<Real> & vRayDirection,
							    int nPlaneNormalAxis = 2 );

};


typedef Frame3<float> Frame3f;
typedef Frame3<double> Frame3d;






}	// namespace g3

