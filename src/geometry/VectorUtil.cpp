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
	Real axisDot = vAlignWith.dot( Vector3<Real>::UnitZ() );

	// compute rotation axis
	Vector3<Real> axisCross( Vector3<Real>::UnitZ().cross( vAlignWith ) );

	Real fInverter = (bInvert) ? (Real)-1 : (Real)1;

	// apply rotation if necessary
	if (axisCross.squaredNorm() > Math<Real>::EPSILON) {

		// compute normalized axis and angle, then create rotation around axis
		axisCross.normalize();
		Real fAngle = Math<Real>::ACos( axisDot / vAlignWith.norm() );
		matrix = Wml::Matrix3<Real>( axisCross, fAngle * fInverter );

	} else if (axisDot < (Real)0) {
		axisCross = Vector3<Real>::UnitX();
		Real fAngle = (Real)180 * Math<Real>::DEG_TO_RAD * fInverter;
		matrix = Wml::Matrix3<Real>( axisCross, fAngle );
	} else {
		matrix = Matrix3<Real>::Identity();
	}
}


template <class Real>
void g3::ComputeAlignAxisMatrix( const Vector3<Real> & vInitial,
								 const Vector3<Real> & vAlignWith, Matrix3<Real> & matrix )
{
	// compute cosine of angle between vectors
	Real axisDot = vAlignWith.dot( vInitial );

	// compute rotation axis
	Vector3<Real> axisCross( vInitial.cross( vAlignWith ) );

	// apply rotation if necessary
	if (axisCross.squaredNorm() > Math<Real>::EPSILON) {

		// compute normalized axis and angle, then create rotation around axis
		axisCross.normalize();
		Real fAngle = Math<Real>::ACos( axisDot / vAlignWith.norm() );
		matrix = Wml::Matrix3<Real>( axisCross, fAngle );

	} else if (axisDot < (Real)0) {

		// find some perpendicular vectors
		Vector3<Real> vPerp1, vPerp2;
		ComputePerpVectors( vInitial, vPerp1, vPerp2 );

		matrix = Wml::Matrix3<Real>( vPerp1, (Real)180 * Math<Real>::DEG_TO_RAD );
	} else {
		matrix = Matrix3<Real>::Identity();
	}
}


template <class Real>
void g3::ComputePerpVectors( const Vector3<Real> & vIn,
							  Vector3<Real> & vOut1, Vector3<Real> & vOut2, bool bInIsNormalized )
{
	Vector3<Real> vPerp(vIn);
	if ( ! bInIsNormalized )
		vPerp.normalize();

	if ( Math<Real>::FAbs(vPerp.x()) >= Math<Real>::FAbs(vPerp.y())
		 &&   Math<Real>::FAbs(vPerp.x()) >= Math<Real>::FAbs(vPerp.z()) )
    {
        vOut1.x() = -vPerp.y();
        vOut1.y() = vPerp.x();
        vOut1.z() = (Real)0.0;
    }
    else
    {
        vOut1.x() = (Real)0.0;
        vOut1.y() = vPerp.z();
        vOut1.z() = -vPerp.y();
    }
	
    vOut1.normalize();
    vOut2 = vPerp.cross(vOut1);	
}


template <class Real>
void g3::ComputePerpVectors( const Vector3<Real> & vNormal,  const Vector3<Real> & vEstX,
							  Vector3<Real> & vOut1, Vector3<Real> & vOut2,
							  bool bInputIsNormalized )
{
	Vector3<Real> n( vNormal );
	Vector3<Real> tan2( vEstX );
	if ( ! bInputIsNormalized ){
		n.normalize();
		tan2.normalize();
	}
	Vector3<Real> tan1 = n.cross(tan2.cross(n));
	tan1.normalize();
	tan2 = n.cross(tan1);

	vOut1 = tan2;
	vOut2 = tan1;
}










template <class Real>
Real g3::VectorAngle( const Vector2<Real> & v1, const Vector2<Real> & v2 )
{
	Real fDot = Clamp(v1.dot(v2), (Real)-1.0, (Real)1.0);
	return (Real)acos(fDot);
}

template <class Real>
Real g3::VectorAngle( const Vector3<Real> & v1, const Vector3<Real> & v2 )
{
	Real fDot = Clamp(v1.dot(v2), (Real)-1.0, (Real)1.0);
	return (Real)acos(fDot);
}

template <class Real>
Real g3::VectorCot( const Vector3<Real> & v1, const Vector3<Real> & v2 )
{
	Real fDot = v1.dot(v2);
	return fDot / (Real)sqrt( v1.dot(v1) * v2.dot(v2) - fDot*fDot );
}





template <class Real>
Vector2<Real> g3::ToUV( const Vector3<Real> & vec, int nUIndex, int nVIndex )
{
	return Vector2<Real>( vec[nUIndex], vec[nVIndex] );
}

template <class Real>
Vector3<Real> g3::To3D( const Vector2<Real> & vec, int nUIndex, int nVIndex )
{
	Vector3<Real> tmp = Vector3<Real>::Zero();
	tmp[nUIndex] = vec.x();
	tmp[nVIndex] = vec.y();
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

    Real fM00 = kV02.dot(kV02);
    Real fM01 = kV02.dot(kV12);
    Real fM11 = kV12.dot(kV12);
    Real fR0 = kV02.dot(kPV2);
    Real fR1 = kV12.dot(kPV2);
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
	Vector3<Real> vCross( edge1.cross(edge2) );

	return (Real)0.5 * vCross.norm();	
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

    Real fM00 = kV02.dot(kV02);
    Real fM01 = kV02.dot(kV12);
    Real fM11 = kV12.dot(kV12);
    Real fR0 = kV02.dot(kPV2);
    Real fR1 = kV12.dot(kPV2);
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
	Real fDot = edge1.dot(edge2);
	return (Real)0.5 * sqrt( edge1.squaredNorm()*edge2.squaredNorm() - fDot*fDot );
}





template <class Real>
Vector3<Real> g3::Normal( const Vector3<Real> & vTriVtx1, 
								const Vector3<Real> & vTriVtx2,
								const Vector3<Real> & vTriVtx3, Real * pArea )
{
	Vector3<Real> edge1( vTriVtx2 - vTriVtx1 );			
	Vector3<Real> edge2( vTriVtx3 - vTriVtx1 );			
	if ( pArea ) {
		Real fDot = edge1.dot(edge2);
		*pArea = (Real)0.5 * sqrt( edge1.squaredNorm()*edge2.squaredNorm() - fDot*fDot );
	}
	edge1.normalize(); edge2.normalize();
	Vector3<Real> vCross( edge1.cross(edge2) );
	vCross.normalize();
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
	vNormal.normalize();
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
	Real s1 = p1.x();
	Real t1 = p1.y();
	Real s2 = p2.x();
	Real t2 = p2.y();
	Real s3 = p3.x();
	Real t3 = p3.y();

	Real A = (Real)0.5 * ( (s2 - s1) * (t3 - t1) - (s3 - s1) * (t2 - t1));
	if ( A > 0 ) {

		Vector3<Real> Ss = 
			(q1 * (t2-t3) + q2 * (t3-t1) + q3 * (t1-t2)) / (2*A);
		Vector3<Real> St = 
			(q1 * (s3-s2) + q2 * (s1-s3) + q3 * (s2-s1)) / (2*A);

		Real a = Ss.dot(Ss);
		Real b = Ss.dot(St);
		Real c = St.dot(St);

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
	Vector3<Real> e1(p2_3D-p1_3D);  e1.normalize();
	Vector3<Real> e2(p3_3D-p1_3D);  e2.normalize();
	Vector3<Real> n(e1.cross(e2));  n.normalize();
	e2 = n.cross(e1);   e2.normalize();
	
	Vector2<Real> p1( Vector2<Real>::Zero() );
	Vector2<Real> p2( (p2_3D-p1_3D).dot(e1), (p2_3D-p1_3D).dot(e2) );
	Vector2<Real> p3( (p3_3D-p1_3D).dot(e1), (p3_3D-p1_3D).dot(e2) );

	Real s1 = p1.x();
	Real t1 = p1.y();
	Real s2 = p2.x();
	Real t2 = p2.y();
	Real s3 = p3.x();
	Real t3 = p3.y();

	Real A = (Real)0.5 * ( (s2 - s1) * (t3 - t1) - (s3 - s1) * (t2 - t1));
	if ( A > 0 ) {

		Vector3<Real> Ss = 
			(q1 * (t2-t3) + q2 * (t3-t1) + q3 * (t1-t2)) / (2*A);
		Vector3<Real> St = 
			(q1 * (s3-s2) + q2 * (s1-s3) + q3 * (s2-s1)) / (2*A);

		Real a = Ss.dot(Ss);
		Real b = Ss.dot(St);
		Real c = St.dot(St);

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
	Real a2 = (v1-v2).squaredNorm();
	Real b2 = (v1-v3).squaredNorm();
	Real c2 = (v2-v3).squaredNorm();
	return (a2+b2 < c2) || (b2+c2 < a2) || (c2+a2 < b2);
}
template <class Real>
bool g3::IsObtuse( const Vector3<Real> & v1, const Vector3<Real> & v2, const Vector3<Real> & v3 )
{
	Real a2 = (v1-v2).squaredNorm();
	Real b2 = (v1-v3).squaredNorm();
	Real c2 = (v2-v3).squaredNorm();
	return (a2+b2 < c2) || (b2+c2 < a2) || (c2+a2 < b2);
}




namespace g3
{



template g3External void ToGLMatrix( const Matrix3<float> & matrix, float glMatrix[16] );
template g3External void ToGLMatrix( const Matrix3<double> & matrix, double glMatrix[16] );

template g3External void ComputeAlignZAxisMatrix( const Vector3<float> & vAlignWith, Matrix3<float> & matrix, bool bInvert );
template g3External void ComputeAlignZAxisMatrix( const Vector3<double> & vAlignWith, Matrix3<double> & matrix, bool bInvert );

template g3External void ComputeAlignAxisMatrix( const Vector3<float> & vInitial, const Vector3<float> & vAlignWith, Matrix3<float> & matrix );
template g3External void ComputeAlignAxisMatrix( const Vector3<double> & vInitial, const Vector3<double> & vAlignWith, Matrix3<double> & matrix );

template g3External void ComputePerpVectors( const Vector3<float> & vIn, Vector3<float> & vOut1, Vector3<float> & vOut2, bool bInIsNormalized );
template g3External void ComputePerpVectors( const Vector3<double> & vIn, Vector3<double> & vOut1, Vector3<double> & vOut2, bool bInIsNormalized );

template g3External void ComputePerpVectors( const Vector3<float> & vNormal,  const Vector3<float> & vEstX, Vector3<float> & vOut1, Vector3<float> & vOut2, bool bInputIsNormalized );
template g3External void ComputePerpVectors( const Vector3<double> & vNormal,  const Vector3<double> & vEstX, Vector3<double> & vOut1, Vector3<double> & vOut2, bool bInputIsNormalized );

template g3External float g3::VectorAngle( const Vector2<float> & v1, const Vector2<float> & v2 );
template g3External double g3::VectorAngle( const Vector2<double> & v1, const Vector2<double> & v2 );
template g3External float g3::VectorAngle( const Vector3<float> & v1, const Vector3<float> & v2 );
template g3External double g3::VectorAngle( const Vector3<double> & v1, const Vector3<double> & v2 );

template g3External float g3::VectorCot( const Vector3<float> & v1, const Vector3<float> & v2 );
template g3External double g3::VectorCot( const Vector3<double> & v1, const Vector3<double> & v2 );



template g3External Vector2<float> ToUV( const Vector3<float> & vec, int nUIndex, int nVIndex );
template g3External Vector2<double> ToUV( const Vector3<double> & vec, int nUIndex, int nVIndex );
template g3External Vector3<float> To3D( const Vector2<float> & vec, int nUIndex, int nVIndex );
template g3External Vector3<double> To3D( const Vector2<double> & vec, int nUIndex, int nVIndex );


template g3External void BarycentricCoords( const Vector3<float> & TriVtx1, const Vector3<float> & TriVtx2,
											 const Vector3<float> & TriVtx3, const Vector3<float> & vVertex,
											 float & fWeight1, float & fWeight2, float & fWeight3 );
template g3External void BarycentricCoords( const Vector3<double> & TriVtx1, const Vector3<double> & TriVtx2,
											 const Vector3<double> & TriVtx3, const Vector3<double> & vVertex,
											 double & fWeight1, double & fWeight2, double & fWeight3 );
template g3External void BarycentricCoords( const Vector2<float> & TriVtx1, const Vector2<float> & TriVtx2,
											 const Vector2<float> & TriVtx3, const Vector2<float> & vVertex,
											 float & fWeight1, float & fWeight2, float & fWeight3 );
template g3External void BarycentricCoords( const Vector2<double> & TriVtx1, const Vector2<double> & TriVtx2,
											 const Vector2<double> & TriVtx3, const Vector2<double> & vVertex,
											 double & fWeight1, double & fWeight2, double & fWeight3 );


template g3External float Area( const Vector3<float> & TriVtx1, const Vector3<float> & TriVtx2,
											 const Vector3<float> & TriVtx3 );
template g3External double Area( const Vector3<double> & TriVtx1, const Vector3<double> & TriVtx2,
											 const Vector3<double> & TriVtx3 );
template g3External float Area( const Vector2<float> & TriVtx1, const Vector2<float> & TriVtx2,
											 const Vector2<float> & TriVtx3 );
template g3External double Area( const Vector2<double> & TriVtx1, const Vector2<double> & TriVtx2,
											 const Vector2<double> & TriVtx3 );


template g3External Vector3<float> Normal( const Vector3<float> & TriVtx1, const Vector3<float> & TriVtx2,
											 const Vector3<float> & TriVtx3, float * pArea );
template g3External Vector3<double> Normal( const Vector3<double> & TriVtx1, const Vector3<double> & TriVtx2,
											 const Vector3<double> & TriVtx3, double * pArea );




template g3External Vector3<float> InterpNormal( const Vector3<float> & vTriVtx1, const Vector3<float> & vTriVtx2,
									  const Vector3<float> & vTriVtx3, const Vector3<float> & vTriNorm1, 
									  const Vector3<float> & vTriNorm2,const Vector3<float> & vTriNorm3,
									  const Vector3<float> & vPointInTri );
template g3External Vector3<double> InterpNormal( const Vector3<double> & vTriVtx1, const Vector3<double> & vTriVtx2,
									  const Vector3<double> & vTriVtx3, const Vector3<double> & vTriNorm1, 
									  const Vector3<double> & vTriNorm2,const Vector3<double> & vTriNorm3,
									  const Vector3<double> & vPointInTri );



template g3External void StretchMetric1( const Vector3<float> & q1, const Vector3<float> & q2, const Vector3<float> & q3,
										  const Vector2<float> & p1, const Vector2<float> & p2, const Vector2<float> & p3,
										  float & MaxSV, float & MinSV, float & L2Norm, float & LInfNorm );
template g3External void StretchMetric1( const Vector3<double> & q1, const Vector3<double> & q2, const Vector3<double> & q3,
										  const Vector2<double> & p1, const Vector2<double> & p2, const Vector2<double> & p3,
										  double & MaxSV, double & MinSV, double & L2Norm, double & LInfNorm );


template g3External  void StretchMetric3( const Vector3<float> & q1, const Vector3<float> & q2, const Vector3<float> & q3,
										  const Vector3<float> & p1, const Vector3<float> & p2, const Vector3<float> & p3,
										  float & MaxSV, float & MinSV, float & L2Norm, float & LInfNorm );
template g3External  void StretchMetric3( const Vector3<double> & q1, const Vector3<double> & q2, const Vector3<double> & q3,
										  const Vector3<double> & p1, const Vector3<double> & p2, const Vector3<double> & p3,
										  double & MaxSV, double & MinSV, double & L2Norm, double & LInfNorm );

template g3External bool IsObtuse( const Vector2<float> & v1, const Vector2<float> & v2, const Vector2<float> & v3 );
template g3External bool IsObtuse( const Vector2<double> & v1, const Vector2<double> & v2, const Vector2<double> & v3 );
template g3External bool IsObtuse( const Vector3<float> & v1, const Vector3<float> & v2, const Vector3<float> & v3 );
template g3External bool IsObtuse( const Vector3<double> & v1, const Vector3<double> & v2, const Vector3<double> & v3 );
}

