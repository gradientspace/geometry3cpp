#include <geometry3PCH.h>
#include "VectorUtil.h"

using namespace g3;


#include <algorithm>
#include <limits>


template <class Real>
void g3::ComputeAlignZAxisMatrix( const Vector3<Real> & vAlignWith,
								  Matrix3<Real> & matrix, bool bInvert )
{
	// compute cosine of angle between vectors
	Real axisDot = vAlignWith.Dot( Vector3<Real>::UNIT_Z );

	// compute rotation axis
	Vector3<Real> axisCross( Vector3<Real>::UNIT_Z.Cross( vAlignWith ) );

	Real fInverter = (bInvert) ? (Real)-1 : (Real)1;

	// apply rotation if necessary
	if (axisCross.SquaredLength() > Math<Real>::EPSILON) {

		// compute normalized axis and angle, then create rotation around axis
		axisCross.Normalize();
		Real fAngle = Math<Real>::ACos( axisDot / vAlignWith.Length() );
		matrix = Matrix3<Real>( axisCross, fAngle * fInverter );

	} else if (axisDot < (Real)0) {
		matrix = Matrix3<Real>( Vector3<Real>::UNIT_X, (Real)180 * Math<Real>::DEG_TO_RAD * fInverter );
	} else {
		matrix = Matrix3<Real>::IDENTITY;
	}
}


template <class Real>
void g3::ComputeAlignAxisMatrix( const Vector3<Real> & vInitial,
								  const Vector3<Real> & vAlignWith, Matrix3<Real> & matrix )
{
	// compute cosine of angle between vectors
	Real axisDot = vAlignWith.Dot( vInitial );

	// compute rotation axis
	Vector3<Real> axisCross( vInitial.Cross( vAlignWith ) );

	// apply rotation if necessary
	if (axisCross.SquaredLength() > Math<Real>::EPSILON) {

		// compute normalized axis and angle, then create rotation around axis
		axisCross.Normalize();
		Real fAngle = Math<Real>::ACos( axisDot / vAlignWith.Length() );
		matrix = Matrix3<Real>( axisCross, fAngle );

	} else if (axisDot < (Real)0) {

		// find some perpendicular vectors
		Vector3<Real> vPerp1, vPerp2;
		ComputePerpVectors( vInitial, vPerp1, vPerp2 );

		matrix = Matrix3<Real>( vPerp1, (Real)180 * Math<Real>::DEG_TO_RAD );
	} else {
		matrix = Matrix3<Real>::IDENTITY;
	}
}


template <class Real>
void g3::ComputePerpVectors( const Vector3<Real> & vIn,
							  Vector3<Real> & vOut1, Vector3<Real> & vOut2, bool bInIsNormalized )
{
	Vector3<Real> vPerp(vIn);
	if ( ! bInIsNormalized )
		vPerp.Normalize();

	if ( Math<Real>::FAbs(vPerp.X()) >= Math<Real>::FAbs(vPerp.Y())
		 &&   Math<Real>::FAbs(vPerp.X()) >= Math<Real>::FAbs(vPerp.Z()) )
    {
        vOut1.X() = -vPerp.Y();
        vOut1.Y() = vPerp.X();
        vOut1.Z() = (Real)0.0;
    }
    else
    {
        vOut1.X() = (Real)0.0;
        vOut1.Y() = vPerp.Z();
        vOut1.Z() = -vPerp.Y();
    }

    vOut1.Normalize();
    vOut2 = vPerp.Cross(vOut1);	
}


template <class Real>
void g3::ComputePerpVectors( const Vector3<Real> & vNormal,  const Vector3<Real> & vEstX,
							  Vector3<Real> & vOut1, Vector3<Real> & vOut2,
							  bool bInputIsNormalized )
{
	Vector3<Real> n( vNormal );
	Vector3<Real> tan2( vEstX );
	if ( ! bInputIsNormalized ){
		n.Normalize();
		tan2.Normalize();
	}
	Vector3<Real> tan1 = n.Cross(tan2.Cross(n));
	tan1.Normalize();
	tan2 = n.Cross(tan1);

	vOut1 = tan2;
	vOut2 = tan1;
}










template <class Real>
Real g3::VectorAngle( const Vector2<Real> & v1, const Vector2<Real> & v2 )
{
	Real fDot = Clamp(v1.Dot(v2), (Real)-1.0, (Real)1.0);
	return (Real)acos(fDot);
}

template <class Real>
Real g3::VectorAngle( const Vector3<Real> & v1, const Vector3<Real> & v2 )
{
	Real fDot = Clamp(v1.Dot(v2), (Real)-1.0, (Real)1.0);
	return (Real)acos(fDot);
}

template <class Real>
Real g3::VectorCot( const Vector3<Real> & v1, const Vector3<Real> & v2 )
{
	Real fDot = v1.Dot(v2);
	return fDot / (Real)sqrt( v1.Dot(v1) * v2.Dot(v2) - fDot*fDot );
}





template <class Real>
Vector2<Real> g3::ToUV( const Vector3<Real> & vec, int nUIndex, int nVIndex )
{
	return Vector2<Real>( vec[nUIndex], vec[nVIndex] );
}

template <class Real>
Vector3<Real> g3::To3D( const Vector2<Real> & vec, int nUIndex, int nVIndex )
{
	Vector3<Real> tmp = Vector3<Real>::ZERO;
	tmp[nUIndex] = vec.X();
	tmp[nVIndex] = vec.Y();
	return tmp;
}





template <class Real>
void g3::ToGLMatrix( const Matrix3<Real> & matrix, Real glMatrix[16] )
{
	for (int r = 0; r < 3; ++r)
		for (int c = 0; c < 4; ++c)
			glMatrix[c*4 + r] = (c < 3) ? matrix(r,c) : 0;
	glMatrix[3] = glMatrix[7] = glMatrix[11] = 0;
	glMatrix[15] = 1;
}








template <class Real>
void g3::BarycentricCoords( const Vector3<Real> & vTriVtx1, 
							 const Vector3<Real> & vTriVtx2,
							 const Vector3<Real> & vTriVtx3,
							 const Vector3<Real> & vVertex,
							 Real & fBary1, Real & fBary2, Real & fBary3 )
{

	Vector3<Real> kV02 = vTriVtx1 - vTriVtx3;
    Vector3<Real> kV12 = vTriVtx2 - vTriVtx3;
    Vector3<Real> kPV2 = vVertex - vTriVtx3;

    Real fM00 = kV02.Dot(kV02);
    Real fM01 = kV02.Dot(kV12);
    Real fM11 = kV12.Dot(kV12);
    Real fR0 = kV02.Dot(kPV2);
    Real fR1 = kV12.Dot(kPV2);
    Real fDet = fM00*fM11 - fM01*fM01;
//    lgASSERT( Math<Real>::FAbs(fDet) > (Real)0.0 );
    Real fInvDet = ((Real)1.0)/fDet;

    fBary1 = (fM11*fR0 - fM01*fR1)*fInvDet;
    fBary2 = (fM00*fR1 - fM01*fR0)*fInvDet;
    fBary3 = (Real)1.0 - fBary1 - fBary2;
}

template <class Real>
Real g3::Area( const Vector3<Real> & vTriVtx1, 
				 const Vector3<Real> & vTriVtx2,
				 const Vector3<Real> & vTriVtx3 )
{
	Vector3<Real> edge1( vTriVtx2 - vTriVtx1 );
	Vector3<Real> edge2( vTriVtx3 - vTriVtx1 );
	Vector3<Real> vCross( edge1.Cross(edge2) );

	return (Real)0.5 * vCross.Length();	
}


template <class Real>
void g3::BarycentricCoords( const Vector2<Real> & vTriVtx1, 
							 const Vector2<Real> & vTriVtx2,
							 const Vector2<Real> & vTriVtx3,
							 const Vector2<Real> & vVertex,
							 Real & fBary1, Real & fBary2, Real & fBary3 )
{

	Vector2<Real> kV02 = vTriVtx1 - vTriVtx3;
    Vector2<Real> kV12 = vTriVtx2 - vTriVtx3;
    Vector2<Real> kPV2 = vVertex - vTriVtx3;

    Real fM00 = kV02.Dot(kV02);
    Real fM01 = kV02.Dot(kV12);
    Real fM11 = kV12.Dot(kV12);
    Real fR0 = kV02.Dot(kPV2);
    Real fR1 = kV12.Dot(kPV2);
    Real fDet = fM00*fM11 - fM01*fM01;
//    lgASSERT( Math<Real>::FAbs(fDet) > (Real)0.0 );
    Real fInvDet = ((Real)1.0)/fDet;

    fBary1 = (fM11*fR0 - fM01*fR1)*fInvDet;
    fBary2 = (fM00*fR1 - fM01*fR0)*fInvDet;
    fBary3 = (Real)1.0 - fBary1 - fBary2;
}


template <class Real>
Real g3::Area( const Vector2<Real> & vTriVtx1, 
				 const Vector2<Real> & vTriVtx2,
				 const Vector2<Real> & vTriVtx3 )
{
	Vector2<Real> edge1( vTriVtx2 - vTriVtx1 );
	Vector2<Real> edge2( vTriVtx3 - vTriVtx1 );
	Real fDot = edge1.Dot(edge2);
	return (Real)0.5 * sqrt( edge1.SquaredLength()*edge2.SquaredLength() - fDot*fDot );
}





template <class Real>
Vector3<Real> g3::Normal( const Vector3<Real> & vTriVtx1, 
								const Vector3<Real> & vTriVtx2,
								const Vector3<Real> & vTriVtx3, Real * pArea )
{
	Vector3<Real> edge1( vTriVtx2 - vTriVtx1 );			
	Vector3<Real> edge2( vTriVtx3 - vTriVtx1 );			
	if ( pArea ) {
		Real fDot = edge1.Dot(edge2);
		*pArea = (Real)0.5 * sqrt( edge1.SquaredLength()*edge2.SquaredLength() - fDot*fDot );
	}
	edge1.Normalize(); edge2.Normalize();
	Vector3<Real> vCross( edge1.Cross(edge2) );
	vCross.Normalize();
	return vCross;
}




template <class Real>
Vector3<Real> g3::InterpNormal( const Vector3<Real> & vTriVtx1, 
									  const Vector3<Real> & vTriVtx2,
									  const Vector3<Real> & vTriVtx3, 
									  const Vector3<Real> & vTriNorm1, 
									  const Vector3<Real> & vTriNorm2,
									  const Vector3<Real> & vTriNorm3,
									  const Vector3<Real> & vPointInTri )
{
	Real fBary[3];
	g3::BarycentricCoords(vTriVtx1, vTriVtx2, vTriVtx3, vPointInTri, fBary[0], fBary[1], fBary[2]);
	Vector3<Real> vNormal( fBary[0]*vTriNorm1 + fBary[1]*vTriNorm1 + fBary[2]*vTriNorm1 );
	vNormal.Normalize();
	return vNormal;
}






template <class Real>
void g3::StretchMetric1( const Vector3<Real> & q1, 
						 const Vector3<Real> & q2,
						 const Vector3<Real> & q3,
						 const Vector2<Real> & p1,
						 const Vector2<Real> & p2,
						 const Vector2<Real> & p3,
						 Real & MaxSV, Real & MinSV, Real & L2Norm, Real & LInfNorm )
{
	Real s1 = p1.X();
	Real t1 = p1.Y();
	Real s2 = p2.X();
	Real t2 = p2.Y();
	Real s3 = p3.X();
	Real t3 = p3.Y();

	Real A = (Real)0.5 * ( (s2 - s1) * (t3 - t1) - (s3 - s1) * (t2 - t1));
	if ( A > 0 ) {

		Vector3<Real> Ss = 
			(q1 * (t2-t3) + q2 * (t3-t1) + q3 * (t1-t2)) / (2*A);
		Vector3<Real> St = 
			(q1 * (s3-s2) + q2 * (s1-s3) + q3 * (s2-s1)) / (2*A);

		Real a = Ss.Dot(Ss);
		Real b = Ss.Dot(St);
		Real c = St.Dot(St);

		Real discrim = (Real)sqrt( (a-c)*(a-c) + 4*b*b );

		MaxSV = (Real)sqrt( (Real)0.5 * ( (a+c) + discrim ) );
		MinSV = (Real)sqrt( (Real)0.5 * ( (a+c) - discrim ) );

		L2Norm = (Real)sqrt( (Real)0.5 * (a+c)  );
		LInfNorm = MaxSV;
	} else {
		MaxSV = MinSV = L2Norm = LInfNorm = std::numeric_limits<Real>::max();
	}

}




template <class Real>
void g3::StretchMetric3( const Vector3<Real> & q1, 
						 const Vector3<Real> & q2,
						 const Vector3<Real> & q3,
						 const Vector3<Real> & p1_3D,
						 const Vector3<Real> & p2_3D,
						 const Vector3<Real> & p3_3D,
						 Real & MaxSV, Real & MinSV, Real & L2Norm, Real & LInfNorm )
{
	// compute plane containing p1/2/3
	Vector3<Real> e1(p2_3D-p1_3D);  e1.Normalize();
	Vector3<Real> e2(p3_3D-p1_3D);  e2.Normalize();
	Vector3<Real> n(e1.Cross(e2));  n.Normalize();
	e2 = n.Cross(e1);   e2.Normalize();
	
	Vector2<Real> p1( Vector2<Real>::ZERO );
	Vector2<Real> p2( (p2_3D-p1_3D).Dot(e1), (p2_3D-p1_3D).Dot(e2) );
	Vector2<Real> p3( (p3_3D-p1_3D).Dot(e1), (p3_3D-p1_3D).Dot(e2) );

	Real s1 = p1.X();
	Real t1 = p1.Y();
	Real s2 = p2.X();
	Real t2 = p2.Y();
	Real s3 = p3.X();
	Real t3 = p3.Y();

	Real A = (Real)0.5 * ( (s2 - s1) * (t3 - t1) - (s3 - s1) * (t2 - t1));
	if ( A > 0 ) {

		Vector3<Real> Ss = 
			(q1 * (t2-t3) + q2 * (t3-t1) + q3 * (t1-t2)) / (2*A);
		Vector3<Real> St = 
			(q1 * (s3-s2) + q2 * (s1-s3) + q3 * (s2-s1)) / (2*A);

		Real a = Ss.Dot(Ss);
		Real b = Ss.Dot(St);
		Real c = St.Dot(St);

		Real discrim = (Real)sqrt( (a-c)*(a-c) + 4*b*b );

		MaxSV = (Real)sqrt( (Real)0.5 * ( (a+c) + discrim ) );
		MinSV = (Real)sqrt( (Real)0.5 * ( (a+c) - discrim ) );

		L2Norm = (Real)sqrt( (Real)0.5 * (a+c)  );
		LInfNorm = MaxSV;
	} else {
		MaxSV = MinSV = L2Norm = LInfNorm = std::numeric_limits<Real>::max();
	}

}




template <class Real>
bool g3::IsObtuse( const Vector2<Real> & v1, const Vector2<Real> & v2, const Vector2<Real> & v3 )
{
	// from http://mathworld.wolfram.com/ObtuseTriangle.html
	Real a2 = (v1-v2).SquaredLength();
	Real b2 = (v1-v3).SquaredLength();
	Real c2 = (v2-v3).SquaredLength();
	return (a2+b2 < c2) || (b2+c2 < a2) || (c2+a2 < b2);
}
template <class Real>
bool g3::IsObtuse( const Vector3<Real> & v1, const Vector3<Real> & v2, const Vector3<Real> & v3 )
{
	Real a2 = (v1-v2).SquaredLength();
	Real b2 = (v1-v3).SquaredLength();
	Real c2 = (v2-v3).SquaredLength();
	return (a2+b2 < c2) || (b2+c2 < a2) || (c2+a2 < b2);
}




namespace g3
{
template void ToGLMatrix( const Matrix3<float> & matrix, float glMatrix[16] );
template void ToGLMatrix( const Matrix3<double> & matrix, double glMatrix[16] );

template void ComputeAlignZAxisMatrix( const Vector3<float> & vAlignWith, Matrix3<float> & matrix, bool bInvert );
template void ComputeAlignZAxisMatrix( const Vector3<double> & vAlignWith, Matrix3<double> & matrix, bool bInvert );

template void ComputeAlignAxisMatrix( const Vector3<float> & vInitial, const Vector3<float> & vAlignWith, Matrix3<float> & matrix );
template void ComputeAlignAxisMatrix( const Vector3<double> & vInitial, const Vector3<double> & vAlignWith, Matrix3<double> & matrix );

template void ComputePerpVectors( const Vector3<float> & vIn, Vector3<float> & vOut1, Vector3<float> & vOut2, bool bInIsNormalized );
template void ComputePerpVectors( const Vector3<double> & vIn, Vector3<double> & vOut1, Vector3<double> & vOut2, bool bInIsNormalized );

template void ComputePerpVectors( const Vector3<float> & vNormal,  const Vector3<float> & vEstX, Vector3<float> & vOut1, Vector3<float> & vOut2, bool bInputIsNormalized );
template void ComputePerpVectors( const Vector3<double> & vNormal,  const Vector3<double> & vEstX, Vector3<double> & vOut1, Vector3<double> & vOut2, bool bInputIsNormalized );



template float g3::VectorAngle( const Vector2<float> & v1, const Vector2<float> & v2 );
template double g3::VectorAngle( const Vector2<double> & v1, const Vector2<double> & v2 );
template float g3::VectorAngle( const Vector3<float> & v1, const Vector3<float> & v2 );
template double g3::VectorAngle( const Vector3<double> & v1, const Vector3<double> & v2 );

template float g3::VectorCot( const Vector3<float> & v1, const Vector3<float> & v2 );
template double g3::VectorCot( const Vector3<double> & v1, const Vector3<double> & v2 );



template Vector2<float> ToUV( const Vector3<float> & vec, int nUIndex, int nVIndex );
template Vector2<double> ToUV( const Vector3<double> & vec, int nUIndex, int nVIndex );
template Vector3<float> To3D( const Vector2<float> & vec, int nUIndex, int nVIndex );
template Vector3<double> To3D( const Vector2<double> & vec, int nUIndex, int nVIndex );


template void BarycentricCoords( const Vector3<float> & TriVtx1, const Vector3<float> & TriVtx2,
											 const Vector3<float> & TriVtx3, const Vector3<float> & vVertex,
											 float & fWeight1, float & fWeight2, float & fWeight3 );
template  void BarycentricCoords( const Vector3<double> & TriVtx1, const Vector3<double> & TriVtx2,
											 const Vector3<double> & TriVtx3, const Vector3<double> & vVertex,
											 double & fWeight1, double & fWeight2, double & fWeight3 );
template  void BarycentricCoords( const Vector2<float> & TriVtx1, const Vector2<float> & TriVtx2,
											 const Vector2<float> & TriVtx3, const Vector2<float> & vVertex,
											 float & fWeight1, float & fWeight2, float & fWeight3 );
template  void BarycentricCoords( const Vector2<double> & TriVtx1, const Vector2<double> & TriVtx2,
											 const Vector2<double> & TriVtx3, const Vector2<double> & vVertex,
											 double & fWeight1, double & fWeight2, double & fWeight3 );


template  float Area( const Vector3<float> & TriVtx1, const Vector3<float> & TriVtx2,
											 const Vector3<float> & TriVtx3 );
template  double Area( const Vector3<double> & TriVtx1, const Vector3<double> & TriVtx2,
											 const Vector3<double> & TriVtx3 );
template  float  Area( const Vector2<float> & TriVtx1, const Vector2<float> & TriVtx2,
											 const Vector2<float> & TriVtx3 );
template  double Area( const Vector2<double> & TriVtx1, const Vector2<double> & TriVtx2,
											 const Vector2<double> & TriVtx3 );


template  Vector3<float> Normal( const Vector3<float> & TriVtx1, const Vector3<float> & TriVtx2,
											 const Vector3<float> & TriVtx3, float * pArea );
template  Vector3<double> Normal( const Vector3<double> & TriVtx1, const Vector3<double> & TriVtx2,
											 const Vector3<double> & TriVtx3, double * pArea );




template Vector3<float> InterpNormal( const Vector3<float> & vTriVtx1, const Vector3<float> & vTriVtx2,
									  const Vector3<float> & vTriVtx3, const Vector3<float> & vTriNorm1, 
									  const Vector3<float> & vTriNorm2,const Vector3<float> & vTriNorm3,
									  const Vector3<float> & vPointInTri );
template Vector3<double> InterpNormal( const Vector3<double> & vTriVtx1, const Vector3<double> & vTriVtx2,
									  const Vector3<double> & vTriVtx3, const Vector3<double> & vTriNorm1, 
									  const Vector3<double> & vTriNorm2,const Vector3<double> & vTriNorm3,
									  const Vector3<double> & vPointInTri );



template  void StretchMetric1( const Vector3<float> & q1, const Vector3<float> & q2, const Vector3<float> & q3,
										  const Vector2<float> & p1, const Vector2<float> & p2, const Vector2<float> & p3,
										  float & MaxSV, float & MinSV, float & L2Norm, float & LInfNorm );
template  void StretchMetric1( const Vector3<double> & q1, const Vector3<double> & q2, const Vector3<double> & q3,
										  const Vector2<double> & p1, const Vector2<double> & p2, const Vector2<double> & p3,
										  double & MaxSV, double & MinSV, double & L2Norm, double & LInfNorm );


template  void StretchMetric3( const Vector3<float> & q1, const Vector3<float> & q2, const Vector3<float> & q3,
										  const Vector3<float> & p1, const Vector3<float> & p2, const Vector3<float> & p3,
										  float & MaxSV, float & MinSV, float & L2Norm, float & LInfNorm );
template  void StretchMetric3( const Vector3<double> & q1, const Vector3<double> & q2, const Vector3<double> & q3,
										  const Vector3<double> & p1, const Vector3<double> & p2, const Vector3<double> & p3,
										  double & MaxSV, double & MinSV, double & L2Norm, double & LInfNorm );

template bool IsObtuse( const Vector2<float> & v1, const Vector2<float> & v2, const Vector2<float> & v3 );
template bool IsObtuse( const Vector2<double> & v1, const Vector2<double> & v2, const Vector2<double> & v3 );
template bool IsObtuse( const Vector3<float> & v1, const Vector3<float> & v2, const Vector3<float> & v3 );
template bool IsObtuse( const Vector3<double> & v1, const Vector3<double> & v2, const Vector3<double> & v3 );
}

