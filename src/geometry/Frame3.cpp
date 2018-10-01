#include <geometry3PCH.h>
#include "Frame3.h"

#include <VectorUtil.h>

using namespace g3;

template <class Real>
Frame3<Real>::~Frame3()
{
}

template <class Real>
Frame3<Real>::Frame3( const Vector3<Real> & vOrigin, const Vector3<Real> & vNormalizedAxis, int nAxis, bool bUseFastPerpVectors)
	: Origin(vOrigin), InverseRotation(Matrix3<Real>::Identity())
{
	if (bUseFastPerpVectors && nAxis == 2) {
		Vector3<Real> kTangent0, kTangent1;
		g3::ComputePerpVectors(vNormalizedAxis, kTangent0, kTangent1, true);
		InverseRotation.row(0) = kTangent0;
		InverseRotation.row(1) = kTangent1;
		InverseRotation.row(2) = vNormalizedAxis;
	} else { 
		AlignAxis(nAxis, vNormalizedAxis);
	}
}


//! create frame at an origin with a given x/y/z vectors
template <class Real>
Frame3<Real>::Frame3( const Vector3<Real> & vOrigin, const Vector3<Real> & vXAxis, const Vector3<Real> & vYAxis, const Vector3<Real> & vZAxis )
	: Origin(vOrigin), InverseRotation(Matrix3<Real>::Identity())
{
	SetFrame( vXAxis, vYAxis, vZAxis );
}


//! create frame with a given rotation
template <class Real>
Frame3<Real>::Frame3(const Matrix3<Real> & mRotation)
	: Origin(Vector3<Real>::Zero()), InverseRotation(mRotation.transpose())
{
}

//! create frame at an origin with a given rotation
template <class Real>
Frame3<Real>::Frame3(const Vector3<Real> & vOrigin, const Matrix3<Real> & mRotation)
	: Origin(vOrigin), InverseRotation(mRotation.transpose())
{
}

template <class Real>
void Frame3<Real>::SetFrame(const Matrix3<Real> & matRotate)
{
	InverseRotation = matRotate.transpose();
}

//! set the reference frame axes
template <class Real>
void Frame3<Real>::SetFrame(const Vector3<Real> & vAxisX, const Vector3<Real> & vAxisY, const Vector3<Real> & vAxisZ)
{
	InverseRotation.row(0) = vAxisX.normalized();
	InverseRotation.row(1) = vAxisY.normalized();
	InverseRotation.row(2) = vAxisZ.normalized();
}


//! rotate selected axis of this frame into toAxis
template <class Real>
void Frame3<Real>::AlignAxis(int nAxis, const Vector3<Real> & toAxis, bool bNormalized)
{
	Matrix3<Real> matAlign;
	g3::ComputeAlignAxisMatrix(
		Axis(nAxis), (bNormalized) ? toAxis : toAxis.normalized(), matAlign);
	Rotate(matAlign, false);
}


//! compute matrix that rotates this frame into toFrame
template <class Real>
void Frame3<Real>::ComputeAlignmentMatrix( const Frame3<Real> & toFrame, Matrix3<Real> & matRotate )
{
	// align vCurFrame.Z() with vDestFrame.Z()
	Matrix3<Real> matAlignZ;
	ComputeAlignAxisMatrix( this->Z(), toFrame.Z(), matAlignZ );
	Frame3<Real> vCopy( *this );
	vCopy.Rotate( matAlignZ );

	// compute rotation angle around vDestFrame.Z()
	Vector3<Real> vX1( toFrame.X() );
	Vector3<Real> vX2( vCopy.X() );
	Matrix3<Real> matAlignX;
	ComputeAlignAxisMatrix( vX2, vX1, matAlignX );

	matRotate = matAlignX * matAlignZ;

	vCopy = Frame3<Real>( *this );
	vCopy.Rotate( matRotate );
	Real fDotX = vCopy.X().dot( toFrame.X() );
	Real fDotY = vCopy.Y().dot( toFrame.Y() );
	Real fDotZ = vCopy.Z().dot( toFrame.Z() );

	// check if Y is flipped - if it is, flip Y...
	//if ( vCopy.Y().dot( toFrame.Y() ) < 0 ) {
	//	Matrix3<Real> matFlip( 1,0,0, 0,-1,0, 0,0,1 );
	//	matRotate = matRotate * matFlip;
	//}
	//// check if Z is flipped - if it is, flip Z...
	if ( vCopy.Z().dot( toFrame.Z() ) < 0 ) {
		Matrix3<Real> matFlip;
		matFlip.row(0) = Vector3<Real>(1, 0, 0);
		matFlip.row(1) = Vector3<Real>(0, 1, 0);
		matFlip.row(2) = Vector3<Real>(0, 0, -1);
		matRotate = matRotate * matFlip;
	}

	// [RMS] does this do anything??
	vCopy = Frame3<Real>( *this );
	vCopy.Rotate( matRotate );
	fDotX = vCopy.X().dot( toFrame.X() );
	fDotY = vCopy.Y().dot( toFrame.Y() );
	fDotZ = vCopy.Z().dot( toFrame.Z() );
}



template <class Real>
void Frame3<Real>::Translate(const Vector3<Real> & vTranslate, bool bRelative) 
{
	if (bRelative)
		Origin += vTranslate;
	else
		Origin = vTranslate;
}



template <class Real>
void Frame3<Real>::Rotate(const Matrix3<Real> & mRotation, bool bReNormalize)
{
	// [RMS] does not seem to work as expected...??
	Vector3<Real> rowX(mRotation * InverseRotation.row(0).transpose());
	Vector3<Real> rowY(mRotation * InverseRotation.row(1).transpose());
	Vector3<Real> rowZ(mRotation * InverseRotation.row(2).transpose());

	if (bReNormalize) {
		Vector3<Real> vCrossY = rowZ.cross(rowX);
		vCrossY.normalize();
		rowY = (vCrossY.dot(rowY) < 0) ? -vCrossY : vCrossY;
		Vector3<Real> vCrossX = rowZ.cross(rowY);
		vCrossX.normalize();
		rowX = (vCrossX.dot(rowX) < 0) ? -vCrossX : vCrossX;
		Vector3<Real> vCrossZ = rowX.cross(rowY);
		vCrossZ.normalize();
		rowZ = (vCrossZ.dot(rowZ) < 0) ? -vCrossZ : vCrossZ;
	}

	InverseRotation.row(0) = rowX;
	InverseRotation.row(1) = rowY;
	InverseRotation.row(2) = rowZ;
}



template <class Real>
void Frame3<Real>::ReNormalize(int nPreserveAxis)
{
	Vector3<Real> rowX(InverseRotation.row(0));
	Vector3<Real> rowY(InverseRotation.row(1));
	Vector3<Real> rowZ(InverseRotation.row(2));

	switch ( nPreserveAxis ) {
		default:
		case -1: {
			Vector3<Real> vCrossY = rowZ.cross(rowX);
			vCrossY.normalize();
			rowY = ( vCrossY.dot(rowY) < 0 ) ? -vCrossY : vCrossY;
			Vector3<Real> vCrossX = rowZ.cross(rowY);
			vCrossX.normalize();
			rowX = ( vCrossX.dot(rowX) < 0 ) ? -vCrossX : vCrossX;
			Vector3<Real> vCrossZ = rowX.cross(rowY);
			vCrossZ.normalize();
			rowZ = ( vCrossZ.dot(rowZ) < 0 ) ? -vCrossZ : vCrossZ;
		} break;

		case 0: {
			rowX.normalize();
			Vector3<Real> vCrossY( rowX.cross(rowZ) );
			vCrossY.normalize();
			rowY = ( vCrossY.dot(rowY) < 0 ) ? -vCrossY : vCrossY;
			Vector3<Real> vCrossZ( rowX.cross(rowY) );
			vCrossZ.normalize();
			rowZ = ( vCrossZ.dot(rowZ) < 0 ) ? -vCrossZ : vCrossZ;
		} break;

		case 1: {
			rowY.normalize();
			Vector3<Real> vCrossX( rowY.cross(rowZ) );
			vCrossX.normalize();
			rowX = ( vCrossX.dot(rowX) < 0 ) ? -vCrossX : vCrossX;
			Vector3<Real> vCrossZ( rowX.cross(rowY) );
			vCrossZ.normalize();
			rowZ = ( vCrossZ.dot(rowZ) < 0 ) ? -vCrossZ : vCrossZ;
		} break;

		case 2: {
			rowZ.normalize();
			Vector3<Real> vCrossX( rowY.cross(rowZ) );
			vCrossX.normalize();
			rowX = ( vCrossX.dot(rowX) < 0 ) ? -vCrossX : vCrossX;
			Vector3<Real> vCrossY( rowX.cross(rowZ) );
			vCrossY.normalize();
			rowY = ( vCrossY.dot(rowY) < 0 ) ? -vCrossY : vCrossY;
		} break;
	}

	InverseRotation.row(0) = rowX;
	InverseRotation.row(1) = rowY;
	InverseRotation.row(2) = rowZ;
}


template <class Real>
bool Frame3<Real>::IntersectRay( const Vector3<Real> & vRayOrigin, const Vector3<Real> & vRayDirection,
							     Vector3<Real> & vRayHit, int nPlaneNormalAxis )
{
	Vector3<Real> N = Axis(nPlaneNormalAxis);
	Real d = -( Origin.dot(N) );
	Real fDenom = vRayDirection.dot(N);
	if ( fDenom < Math<Real>::ZERO_TOLERANCE )
		return false;

	Real t = - ( vRayOrigin.dot(N) + d ) / fDenom;
	vRayHit = vRayOrigin + t * vRayDirection;
	return true;
}
template <class Real>
Vector3<Real> Frame3<Real>::IntersectRay( const Vector3<Real> & vRayOrigin, const Vector3<Real> & vRayDirection,
							int nPlaneNormalAxis )
{
	Vector3<Real> vHit = Vector3<Real>::Zero();
	bool bOK = IntersectRay(vRayOrigin, vRayDirection, vHit, nPlaneNormalAxis);
    // [TODO] return invalid vector here instead of ZERO
    return (bOK) ? vHit : Vector3<Real>::Zero();
}






//----------------------------------------------------------------------------
// explicit instantiation
//----------------------------------------------------------------------------
namespace g3
{
template class Frame3<float>;
template class Frame3<double>;

}
//----------------------------------------------------------------------------
