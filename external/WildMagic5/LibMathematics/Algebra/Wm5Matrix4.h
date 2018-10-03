// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.1 (2010/10/01)

#ifndef WM5MATRIX4_H
#define WM5MATRIX4_H

#include "Wm5MathematicsLIB.h"
#include "Wm5Table.h"
#include "Wm5Vector3.h"
#include "Wm5Vector4.h"

namespace Wm5
{

template <typename Real>
class Matrix4 : public Table<4,4,Real>
{
public:
    // If makeZero is 'true', create the zero matrix; otherwise, create the
    // identity matrix.
    Matrix4 (bool makeZero = true);

    // Copy constructor.
    Matrix4 (const Matrix4& mat);

    // Input mrc is in row r, column c.
    Matrix4 (
        Real m00, Real m01, Real m02, Real m03,
        Real m10, Real m11, Real m12, Real m13,
        Real m20, Real m21, Real m22, Real m23,
        Real m30, Real m31, Real m32, Real m33);

    // Create a matrix from an array of numbers.  The input array is
    // interpreted based on the bool input as
    //   true:  entry[0..15]={m00,m01,m02,m03,m10,m11,m12,m13,m20,m21,m22,
    //                        m23,m30,m31,m32,m33} [row major]
    //   false: entry[0..15]={m00,m10,m20,m30,m01,m11,m21,m31,m02,m12,m22,
    //                        m32,m03,m13,m23,m33} [col major]
    Matrix4 (const Real entry[16], bool rowMajor);

    // Assignment.
    Matrix4& operator= (const Matrix4& mat);

    // Create various matrices.
    void MakeZero ();
    void MakeIdentity ();

    // Arithmetic operations.
    Matrix4 operator+ (const Matrix4& mat) const;
    Matrix4 operator- (const Matrix4& mat) const;
    Matrix4 operator* (Real scalar) const;
    Matrix4 operator/ (Real scalar) const;
    Matrix4 operator- () const;

    // Arithmetic updates.
    Matrix4& operator+= (const Matrix4& mat);
    Matrix4& operator-= (const Matrix4& mat);
    Matrix4& operator*= (Real scalar);
    Matrix4& operator/= (Real scalar);

    // M*vec
    Vector4<Real> operator* (const Vector4<Real>& vec) const;

    // u^T*M*v
    Real QForm (const Vector4<Real>& u, const Vector4<Real>& v) const;

    // M^T
    Matrix4 Transpose () const;

    // M*mat
    Matrix4 operator* (const Matrix4& mat) const;

    // M^T*mat
    Matrix4 TransposeTimes (const Matrix4& mat) const;

    // M*mat^T
    Matrix4 TimesTranspose (const Matrix4& mat) const;

    // M^T*mat^T
    Matrix4 TransposeTimesTranspose (const Matrix4& mat) const;

    // Other operations.
    Matrix4 Inverse (const Real epsilon = (Real)0) const;
    Matrix4 Adjoint () const;
    Real Determinant () const;

    // Projection matrices onto a specified plane (containing an 'origin'
    // point and a unit-length 'normal').
    void MakeObliqueProjection (const Vector3<Real>& normal,
        const Vector3<Real>& origin, const Vector3<Real>& direction);

    void MakePerspectiveProjection (const Vector3<Real>& normal,
        const Vector3<Real>& origin, const Vector3<Real>& eye);

    // Reflection matrix through a specified plane.
    void MakeReflection (const Vector3<Real>& normal,
        const Vector3<Real>& origin);

    // Special matrices.
    WM5_MATHEMATICS_ITEM static const Matrix4 ZERO;
    WM5_MATHEMATICS_ITEM static const Matrix4 IDENTITY;

protected:
    using Table<4,4,Real>::mEntry;

	// [geometry3]
public:
	using EMatrix4 = Eigen::Matrix<Real, 4, 4>;
	operator EMatrix4() const {
		return Eigen::Map<EMatrix4>((Real *)this);
	}
	Matrix4(const EMatrix4 & mat) {
		const Real * p = mat.data();
		mEntry[0] = p[0];
		mEntry[1] = p[1];
		mEntry[2] = p[2];
		mEntry[3] = p[3];
		mEntry[4] = p[4];
		mEntry[5] = p[5];
		mEntry[6] = p[6];
		mEntry[7] = p[7];
		mEntry[8] = p[8];
		mEntry[9] = p[9];
		mEntry[10] = p[10];
		mEntry[11] = p[11];
		mEntry[12] = p[12];
		mEntry[13] = p[13];
		mEntry[14] = p[14];
		mEntry[15] = p[15];
	}
	// [geometry3]

};

// c * M
template <typename Real>
inline Matrix4<Real> operator* (Real scalar, const Matrix4<Real>& mat);

// v^T * M
template <typename Real>
inline Vector4<Real> operator* (const Vector4<Real>& vec,
    const Matrix4<Real>& mat);

#include "Wm5Matrix4.inl"

typedef Matrix4<float> Matrix4f;
typedef Matrix4<double> Matrix4d;

}

#endif
