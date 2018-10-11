#pragma once

#include <g3types.h>

/*
 * VectorUtil contains a set of utility functions for matrix/vector math.
 */


namespace g3
{

	/*
	 * if bToZAxis is false, compute matrix that rotates Z axis into vAlignWith
	 * if bToZAxis is true, compute matrix that rotates vAlignWith into Z axis
	 */
	template <class Real>
	void ComputeAlignZAxisMatrix( const Vector3<Real> & vAlignWith,
								  Matrix3<Real> & matrix, bool bToZAxis = false );

	template <class Real>
	void ComputeAlignAxisMatrix( const Vector3<Real> & vInitial,
								 const Vector3<Real> & vAlignWith, Matrix3<Real> & matrix );

	//! compute vectors in a plane perpendicular to vIn
	template <class Real>
	void ComputePerpVectors( const Vector3<Real> & vIn,
							 Vector3<Real> & vOut1, Vector3<Real> & vOut2,
							 bool bInIsNormalized = false);

	//! compute tangent vectors in plane perp to vNormal, using non-orthogonal vEstX as estimate of vOut1
	template <class Real>
	void ComputePerpVectors( const Vector3<Real> & vNormal,  const Vector3<Real> & vEstX,
							 Vector3<Real> & vOut1, Vector3<Real> & vOut2,
							 bool bInputIsNormalized = false);


	template <class Real>
	void ToGLMatrix( const Matrix3<Real> & matrix, Real glMatrix[16] );

	template<class Real>
	Vector2<Real> ToUV( const Vector3<Real> & vec, int nUIndex, int nVIndex );
	template<class Real>
	Vector3<Real> To3D( const Vector2<Real> & vec, int nUIndex, int nVIndex );


	template <class Real>
	Real VectorAngleR( const Vector2<Real> & v1, const Vector2<Real> & v2 );
	template <class Real>
	Real VectorAngleR( const Vector3<Real> & v1, const Vector3<Real> & v2 );

	template <class Real>
	Real VectorAngleD(const Vector2<Real> & v1, const Vector2<Real> & v2);
	template <class Real>
	Real VectorAngleD(const Vector3<Real> & v1, const Vector3<Real> & v2);

	template <class Real>
	Real VectorCot( const Vector3<Real> & v1, const Vector3<Real> & v2 );


	template <class Real>
	Vector2<Real> Lerp(const Vector2<Real> & v1, const Vector2<Real> & v2, Real t);
	template <class Real>
	Vector3<Real> Lerp(const Vector3<Real> & v1, const Vector3<Real> & v2, Real t);

	template <class Real>
	void BarycentricCoords( const Vector3<Real> & vTriVtx1, 
										const Vector3<Real> & vTriVtx2,
										const Vector3<Real> & vTriVtx3,
										const Vector3<Real> & vVertex,
										Real & fBary1, Real & fBary2, Real & fBary3 );

	template <class Real>
	Real Area( const Vector3<Real> & vTriVtx1, 
						   const Vector3<Real> & vTriVtx2,
						   const Vector3<Real> & vTriVtx3 );

	template <class Real>
	void BarycentricCoords( const Vector2<Real> & vTriVtx1, 
										const Vector2<Real> & vTriVtx2,
										const Vector2<Real> & vTriVtx3,
										const Vector2<Real> & vVertex,
										Real & fBary1, Real & fBary2, Real & fBary3 );

	template <class Real>
	Real Area( const Vector2<Real> & vTriVtx1, 
						   const Vector2<Real> & vTriVtx2,
						   const Vector2<Real> & vTriVtx3 );


	template <class Real>
	Vector3<Real> Normal( const Vector3<Real> & vTriVtx1, 
						  const Vector3<Real> & vTriVtx2,
						  const Vector3<Real> & vTriVtx3, Real * pArea = nullptr );


	template <class Real>
	Vector3<Real> InterpNormal( const Vector3<Real> & vTriVtx1, 
									 const Vector3<Real> & vTriVtx2,
									 const Vector3<Real> & vTriVtx3, 
									 const Vector3<Real> & vTriNorm1, 
									 const Vector3<Real> & vTriNorm2,
									 const Vector3<Real> & vTriNorm3,
									 const Vector3<Real> & vPointInTri );
	


	//! This metric is from Texture Mapping Progressive Meshes, Sander et al, Siggraph 2001
	template <class Real>
	void StretchMetric1( const Vector3<Real> & vTriVtx1, 
									 const Vector3<Real> & vTriVtx2,
									 const Vector3<Real> & vTriVtx3,
									 const Vector2<Real> & vVtxParam1,
									 const Vector2<Real> & vVtxParam2,
									 const Vector2<Real> & vVtxParam3,
									 Real & MaxSV, Real & MinSV, Real & L2Norm, Real & LInfNorm );

	template <class Real>
	void StretchMetric3( const Vector3<Real> & vTriVtx1, 
									 const Vector3<Real> & vTriVtx2,
									 const Vector3<Real> & vTriVtx3,
									 const Vector3<Real> & vVtxParam1,
									 const Vector3<Real> & vVtxParam2,
									 const Vector3<Real> & vVtxParam3,
									 Real & MaxSV, Real & MinSV, Real & L2Norm, Real & LInfNorm );


	template <class Real>
	bool IsObtuse( const Vector2<Real> & v1, const Vector2<Real> & v2, const Vector2<Real> & v3 );
	template <class Real>
	bool IsObtuse( const Vector3<Real> & v1, const Vector3<Real> & v2, const Vector3<Real> & v3 );



	/*
	 * inline utilities
	 */

	inline bool IsFinite(const Vector3d & v) {
		return _finite(v.x()) && _finite(v.y()) && _finite(v.z());
	}

	inline Vector2f d2f(const Vector2d & v) { 
		return Vector2f((float)v[0], (float)v[1]); 
	}
	inline Vector3f d2f(const Vector3d & v) {
		return Vector3f((float)v[0], (float)v[1], (float)v[2]); 
	}
	inline Vector2d f2d(const Vector2f & v) {
		return Vector2d((double)v[0], (double)v[1]); 
	}
	inline Vector3d f2d(const Vector3f & v) {
		return Vector3d((double)v[0], (double)v[1], (double)v[2]); 
	}


	inline Matrix2f d2f(const Matrix2d & v) {
		return v.cast<float>();
	}
	inline Matrix3f d2f(const Matrix3d & v) {
		return v.cast<float>();
	}
	inline Matrix3d f2d(const Matrix3f & v) {
		return v.cast<double>();
	}
	inline Matrix2d f2d(const Matrix2f & v) {
		return v.cast<double>();
	}


	inline AxisAlignedBox2f d2f(const AxisAlignedBox2d & v) {
		return AxisAlignedBox2f((float)v.Min[0], (float)v.Max[0], (float)v.Min[1], (float)v.Max[1] ); 
	}
	inline AxisAlignedBox3f d2f(const AxisAlignedBox3d & v) {
		return AxisAlignedBox3f((float)v.Min[0], (float)v.Max[0], (float)v.Min[1], (float)v.Max[1], (float)v.Min[2], (float)v.Max[2] ); 
	}
	inline AxisAlignedBox2d f2d(const AxisAlignedBox2f & v) {
		return AxisAlignedBox2d((double)v.Min[0], (double)v.Max[0], (double)v.Min[1], (double)v.Max[1]); 
	}
	inline AxisAlignedBox3d f2d(const AxisAlignedBox3f & v) {
		return AxisAlignedBox3d((double)v.Min[0], (double)v.Max[0], (double)v.Min[1], (double)v.Max[1], (double)v.Min[2], (double)v.Max[2]);
	}


	template<class Real>
	inline Real Clamp(const Real & fValue, const Real & fMin, const Real & fMax) {
		if (fValue < fMin)
			return fMin;
		else if (fValue > fMax)
			return fMax;
		else
			return fValue;
	}


	inline void array3f_add( float * pBuffer, unsigned int nIndex, const float * pAdd ) {
		pBuffer[3*nIndex] += pAdd[0]; pBuffer[3*nIndex+1] += pAdd[1]; pBuffer[3*nIndex+2] += pAdd[2];
	}
	inline void array3f_normalize( float * pBuffer, unsigned int nIndex, float fEpsilon = 0.0f ) {
		auto v = Vector3f(&pBuffer[3 * nIndex]).normalized();
		pBuffer[3*nIndex] = v[0]; pBuffer[3*nIndex+1] = v[1]; pBuffer[3*nIndex+2] = v[2];
	}
	inline void vectorf_push( std::vector<float> & v, const g3::Vector3f & p ) {
		v.push_back(p[0]); v.push_back(p[1]); v.push_back(p[2]);
	}
	inline void vectori_push( std::vector<unsigned int> & v, const g3::Vector2i & p ) {
		v.push_back(p[0]); v.push_back(p[1]);
	}
	inline void vectori_push( std::vector<unsigned int> & v, const g3::Vector3i & p ) {
		v.push_back(p[0]); v.push_back(p[1]); v.push_back(p[2]);
	}

	template <typename DerivedA, typename DerivedB>
	inline bool EpsilonEqual(const Eigen::MatrixBase<DerivedA> & m1, const Eigen::MatrixBase<DerivedB> & m2, double eps)
	{
		return (m1 - m2).cwiseAbs().maxCoeff() <= eps;
	}

	template <typename DerivedA, typename T>
	inline bool Contains(const Eigen::MatrixBase<DerivedA> & m1, T value)
	{
		int n = (int)m1.array().size();
		for (int k = 0; k < n; ++k) {
			if (m1.array()[k] == value)
				return true;
		}
		return false;
	}

	template<typename T1, typename T2>
	inline bool ContainsKey(const std::map<T1,T2> & dict, T1 value) {
		return dict.find(value) != dict.end();
	}

	template<typename T>
	inline bool Contains(const std::vector<T> & vec, T value) {
		return std::find(vec.begin(), vec.end(), value) != vec.end();
	}
	template<typename T>
	inline bool Remove(std::vector<T> & vec, T value) {
		auto itr = std::find(vec.begin(), vec.end(), value);
		if (itr == vec.end())
			return false;
		vec.erase(itr);
		return true;
	}

}  // namespace g3


