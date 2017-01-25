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
	: Origin(vOrigin), InverseRotation(false)
{
	if (bUseFastPerpVectors && nAxis == 2) {
		Vector3<Real> kTangent0, kTangent1;
		g3::ComputePerpVectors(vNormalizedAxis, kTangent0, kTangent1, true);
		InverseRotation.SetRow(0, kTangent0);
		InverseRotation.SetRow(1, kTangent1);
		InverseRotation.SetRow(2, vNormalizedAxis);
	} else { 
		AlignAxis(nAxis, vNormalizedAxis);
	}
}


//! create frame at an origin with a given x/y/z vectors
template <class Real>
Frame3<Real>::Frame3( const Vector3<Real> & vOrigin, const Vector3<Real> & vXAxis, const Vector3<Real> & vYAxis, const Vector3<Real> & vZAxis )
	: Origin(vOrigin), InverseRotation(false)
{
	SetFrame( vXAxis, vYAxis, vZAxis );
}


//! create frame with a given rotation
template <class Real>
Frame3<Real>::Frame3(const Matrix3<Real> & mRotation)
	: Origin(Vector3<Real>::ZERO), InverseRotation(mRotation.Transpose())
{
}

//! create frame at an origin with a given rotation
template <class Real>
Frame3<Real>::Frame3(const Vector3<Real> & vOrigin, const Matrix3<Real> & mRotation)
	: Origin(vOrigin), InverseRotation(mRotation.Transpose())
{
}

template <class Real>
void Frame3<Real>::SetFrame(const Matrix3<Real> & matRotate)
{
	InverseRotation = matRotate.Transpose();
}


//! set the reference frame axes
template <class Real>
void Frame3<Real>::SetFrame(const Vector3<Real> & vAxisX, const Vector3<Real> & vAxisY, const Vector3<Real> & vAxisZ)
{
	InverseRotation.SetRow(0, vAxisX.Normalized() );
	InverseRotation.SetRow(1, vAxisY.Normalized());
	InverseRotation.SetRow(2, vAxisZ.Normalized());
}


//! rotate selected axis of this frame into toAxis
template <class Real>
void Frame3<Real>::AlignAxis(int nAxis, const Vector3<Real> & toAxis, bool bNormalized)
{
	Matrix3<Real> matAlign;
	g3::ComputeAlignAxisMatrix(
		Axis(nAxis), (bNormalized) ? toAxis : toAxis.Normalized(), matAlign);
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
	Real fDotX = vCopy.X().Dot( toFrame.X() );
	Real fDotY = vCopy.Y().Dot( toFrame.Y() );
	Real fDotZ = vCopy.Z().Dot( toFrame.Z() );

	// check if Y is flipped - if it is, flip Y...
	//if ( vCopy.Y().Dot( toFrame.Y() ) < 0 ) {
	//	Matrix3<Real> matFlip( 1,0,0, 0,-1,0, 0,0,1 );
	//	matRotate = matRotate * matFlip;
	//}
	//// check if Z is flipped - if it is, flip Z...
	if ( vCopy.Z().Dot( toFrame.Z() ) < 0 ) {
		Matrix3<Real> matFlip( 1,0,0, 0,1,0, 0,0,-1 );
		matRotate = matRotate * matFlip;
	}

	// [RMS] does this do anything??
	vCopy = Frame3<Real>( *this );
	vCopy.Rotate( matRotate );
	fDotX = vCopy.X().Dot( toFrame.X() );
	fDotY = vCopy.Y().Dot( toFrame.Y() );
	fDotZ = vCopy.Z().Dot( toFrame.Z() );
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

	Vector3<Real> rowX(mRotation * InverseRotation.GetRow(0));
	Vector3<Real> rowY(mRotation * InverseRotation.GetRow(1));
	Vector3<Real> rowZ(mRotation * InverseRotation.GetRow(2));

	if (bReNormalize) {
		Vector3<Real> vCrossY = rowZ.Cross(rowX);
		vCrossY.Normalize();
		rowY = (vCrossY.Dot(rowY) < 0) ? -vCrossY : vCrossY;
		Vector3<Real> vCrossX = rowZ.Cross(rowY);
		vCrossX.Normalize();
		rowX = (vCrossX.Dot(rowX) < 0) ? -vCrossX : vCrossX;
		Vector3<Real> vCrossZ = rowX.Cross(rowY);
		vCrossZ.Normalize();
		rowZ = (vCrossZ.Dot(rowZ) < 0) ? -vCrossZ : vCrossZ;
	}

	InverseRotation.SetRow(0, rowX);
	InverseRotation.SetRow(1, rowY);
	InverseRotation.SetRow(2, rowZ);
}



template <class Real>
void Frame3<Real>::ReNormalize(int nPreserveAxis)
{
	Vector3<Real> rowX(InverseRotation.GetRow(0));
	Vector3<Real> rowY(InverseRotation.GetRow(1));
	Vector3<Real> rowZ(InverseRotation.GetRow(2));

	switch ( nPreserveAxis ) {
		default:
		case -1: {
			Vector3<Real> vCrossY = rowZ.Cross(rowX);
			vCrossY.Normalize();
			rowY = ( vCrossY.Dot(rowY) < 0 ) ? -vCrossY : vCrossY;
			Vector3<Real> vCrossX = rowZ.Cross(rowY);
			vCrossX.Normalize();
			rowX = ( vCrossX.Dot(rowX) < 0 ) ? -vCrossX : vCrossX;
			Vector3<Real> vCrossZ = rowX.Cross(rowY);
			vCrossZ.Normalize();
			rowZ = ( vCrossZ.Dot(rowZ) < 0 ) ? -vCrossZ : vCrossZ;
		} break;

		case 0: {
			rowX.Normalize();
			Vector3<Real> vCrossY( rowX.Cross(rowZ) );
			vCrossY.Normalize();
			rowY = ( vCrossY.Dot(rowY) < 0 ) ? -vCrossY : vCrossY;
			Vector3<Real> vCrossZ( rowX.Cross(rowY) );
			vCrossZ.Normalize();
			rowZ = ( vCrossZ.Dot(rowZ) < 0 ) ? -vCrossZ : vCrossZ;
		} break;

		case 1: {
			rowY.Normalize();
			Vector3<Real> vCrossX( rowY.Cross(rowZ) );
			vCrossX.Normalize();
			rowX = ( vCrossX.Dot(rowX) < 0 ) ? -vCrossX : vCrossX;
			Vector3<Real> vCrossZ( rowX.Cross(rowY) );
			vCrossZ.Normalize();
			rowZ = ( vCrossZ.Dot(rowZ) < 0 ) ? -vCrossZ : vCrossZ;
		} break;

		case 2: {
			rowZ.Normalize();
			Vector3<Real> vCrossX( rowY.Cross(rowZ) );
			vCrossX.Normalize();
			rowX = ( vCrossX.Dot(rowX) < 0 ) ? -vCrossX : vCrossX;
			Vector3<Real> vCrossY( rowX.Cross(rowZ) );
			vCrossY.Normalize();
			rowY = ( vCrossY.Dot(rowY) < 0 ) ? -vCrossY : vCrossY;
		} break;
	}

	InverseRotation.SetRow(0, rowX);
	InverseRotation.SetRow(1, rowY);
	InverseRotation.SetRow(2, rowZ);
}


template <class Real>
bool Frame3<Real>::IntersectRay( const Vector3<Real> & vRayOrigin, const Vector3<Real> & vRayDirection,
							     Vector3<Real> & vRayHit, int nPlaneNormalAxis )
{
	Vector3<Real> N = Axis(nPlaneNormalAxis);
	Real d = -( Origin.Dot(N) );
	Real fDenom = vRayDirection.Dot(N);
	if ( fDenom < Math<Real>::ZERO_TOLERANCE )
		return false;

	Real t = - ( vRayOrigin.Dot(N) + d ) / fDenom;
	vRayHit = vRayOrigin + t * vRayDirection;
	return true;
}
template <class Real>
Vector3<Real> Frame3<Real>::IntersectRay( const Vector3<Real> & vRayOrigin, const Vector3<Real> & vRayDirection,
							int nPlaneNormalAxis )
{
	Vector3<Real> vHit = Vector3<Real>::ZERO;
	bool bOK = IntersectRay(vRayOrigin, vRayDirection, vHit, nPlaneNormalAxis);
    // [TODO] return invalid vector here instead of ZERO
    return (bOK) ? vHit : Vector3<Real>::ZERO;
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
