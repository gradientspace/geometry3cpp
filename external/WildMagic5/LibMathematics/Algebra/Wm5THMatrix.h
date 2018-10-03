// Geometric Tools, LLC
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.16.0 (2017/08/24)

#ifndef WM5THMATRIX_H
#define WM5THMATRIX_H

#include "Wm5MathematicsLIB.h"
#include "Wm5TAPoint.h"
#include "Wm5Matrix3.h"

namespace Wm5
{

template <typename Real>
class THMatrix
{
public:
    // Construction and destruction.
    THMatrix()
    {
        // uninitialized members
    }

    THMatrix(const THMatrix& mat)
    {
        for (int i = 0; i < 16; ++i)
        {
            mEntry[i] = mat.mEntry[i];
        }
    }

    THMatrix(Matrix3<Real> const& mat)
    {
        mEntry[0] = mat[0][0];
        mEntry[1] = mat[0][1];
        mEntry[2] = mat[0][2];
        mEntry[3] = (Real)0;
        mEntry[4] = mat[1][0];
        mEntry[5] = mat[1][1];
        mEntry[6] = mat[1][2];
        mEntry[7] = (Real)0;
        mEntry[8] = mat[2][0];
        mEntry[9] = mat[2][1];
        mEntry[10] = mat[2][2];
        mEntry[11] = (Real)0;
        mEntry[12] = (Real)0;
        mEntry[13] = (Real)0;
        mEntry[14] = (Real)0;
        mEntry[15] = (Real)1;
    }

    // If makeZero is 'true', create the zero matrix; otherwise, create the
    // identity matrix.
    THMatrix(bool makeZero)
    {
        if (makeZero)
        {
            MakeZero();
        }
        else
        {
            MakeIdentity();
        }
    }

    // Input mrc is in row r, column c.
    THMatrix(
        Real m00, Real m01, Real m02, Real m03,
        Real m10, Real m11, Real m12, Real m13,
        Real m20, Real m21, Real m22, Real m23,
        Real m30, Real m31, Real m32, Real m33)
    {
        mEntry[0] = m00;
        mEntry[1] = m01;
        mEntry[2] = m02;
        mEntry[3] = m03;
        mEntry[4] = m10;
        mEntry[5] = m11;
        mEntry[6] = m12;
        mEntry[7] = m13;
        mEntry[8] = m20;
        mEntry[9] = m21;
        mEntry[10] = m22;
        mEntry[11] = m23;
        mEntry[12] = m30;
        mEntry[13] = m31;
        mEntry[14] = m32;
        mEntry[15] = m33;
    }

    // Create a matrix from an array of numbers.  The input array is
    // interpreted based on the bool input as
    //   true:  entry[0..15]={m00,m01,m02,m03,m10,m11,m12,m13,m20,m21,m22,
    //                        m23,m30,m31,m32,m33} [row major]
    //   false: entry[0..15]={m00,m10,m20,m30,m01,m11,m21,m31,m02,m12,m22,
    //                        m32,m03,m13,m23,m33} [col major]
    THMatrix(const Real* entry, bool rowMajor)
    {
        if (rowMajor)
        {
            mEntry[0] = entry[0];
            mEntry[1] = entry[1];
            mEntry[2] = entry[2];
            mEntry[3] = entry[3];
            mEntry[4] = entry[4];
            mEntry[5] = entry[5];
            mEntry[6] = entry[6];
            mEntry[7] = entry[7];
            mEntry[8] = entry[8];
            mEntry[9] = entry[9];
            mEntry[10] = entry[10];
            mEntry[11] = entry[11];
            mEntry[12] = entry[12];
            mEntry[13] = entry[13];
            mEntry[14] = entry[14];
            mEntry[15] = entry[15];
        }
        else
        {
            mEntry[0] = entry[0];
            mEntry[1] = entry[4];
            mEntry[2] = entry[8];
            mEntry[3] = entry[12];
            mEntry[4] = entry[1];
            mEntry[5] = entry[5];
            mEntry[6] = entry[9];
            mEntry[7] = entry[13];
            mEntry[8] = entry[2];
            mEntry[9] = entry[6];
            mEntry[10] = entry[10];
            mEntry[11] = entry[14];
            mEntry[12] = entry[3];
            mEntry[13] = entry[7];
            mEntry[14] = entry[11];
            mEntry[15] = entry[15];
        }
    }

    // Create matrices based on point and vector input.  The bool is
    // interpreted as
    //   true: inputs are columns of the matrix
    //   false: inputs are rows of the matrix
    THMatrix(const Real* tuple0, const Real* tuple1, const Real* tuple2,
        const Real* tuple3, bool columns)
    {
        if (columns)
        {
            mEntry[0] = tuple0[0];
            mEntry[1] = tuple1[0];
            mEntry[2] = tuple2[0];
            mEntry[3] = tuple3[0];
            mEntry[4] = tuple0[1];
            mEntry[5] = tuple1[1];
            mEntry[6] = tuple2[1];
            mEntry[7] = tuple3[1];
            mEntry[8] = tuple0[2];
            mEntry[9] = tuple1[2];
            mEntry[10] = tuple2[2];
            mEntry[11] = tuple3[2];
            mEntry[12] = tuple0[3];
            mEntry[13] = tuple1[3];
            mEntry[14] = tuple2[3];
            mEntry[15] = tuple3[3];
        }
        else
        {
            mEntry[0] = tuple0[0];
            mEntry[1] = tuple0[1];
            mEntry[2] = tuple0[2];
            mEntry[3] = tuple0[3];
            mEntry[4] = tuple1[0];
            mEntry[5] = tuple1[1];
            mEntry[6] = tuple1[2];
            mEntry[7] = tuple1[3];
            mEntry[8] = tuple2[0];
            mEntry[9] = tuple2[1];
            mEntry[10] = tuple2[2];
            mEntry[11] = tuple2[3];
            mEntry[12] = tuple3[0];
            mEntry[13] = tuple3[1];
            mEntry[14] = tuple3[2];
            mEntry[15] = tuple3[3];
        }
    }

    // Create a diagonal matrix.
    THMatrix(Real m00, Real m11, Real m22)
    {
        MakeDiagonal(m00, m11, m22);
    }

    // Create a rotation matrix (positive angle -> counterclockwise).  The
    // angle must be in radians, not degrees.
    THMatrix(const TAVector<Real>& axis, Real angle)
    {
        MakeRotation(axis, angle);
    }

    ~THMatrix()
    {
    }

    // Implicit conversions.  The upper 3x3 block of THMatrix is copied to
    // the Matrix3<Real> object.
    operator Matrix3<Real>() const
    {
        return Matrix3<Real>
        (
            mEntry[0], mEntry[1], mEntry[2],
            mEntry[4], mEntry[5], mEntry[6],
            mEntry[8], mEntry[9], mEntry[10]
        );
    }

    // Coordinate access.
    inline operator const Real* () const
    {
        return mEntry;
    }

    inline operator Real* ()
    {
        return mEntry;
    }

    inline const Real* operator[] (int row) const
    {
        return &mEntry[4 * row];
    }

    inline Real* operator[] (int row)
    {
        return &mEntry[4 * row];
    }

    inline const Real& operator() (int row, int column) const
    {
        return mEntry[column + 4 * row];
    }

    inline Real& operator() (int row, int column)
    {
        return mEntry[column + 4 * row];
    }

    void SetRow(int row, const THPoint<Real>& hpnt)
    {
        int i = 4 * row;
        Real* entry = &mEntry[i];
        *entry = hpnt[0];  ++entry;
        *entry = hpnt[1];  ++entry;
        *entry = hpnt[2];  ++entry;
        *entry = hpnt[3];
    }

    void GetRow(int row, THPoint<Real>& hpnt) const
    {
        int i = 4 * row;
        const Real* entry = &mEntry[i];
        hpnt[0] = *entry;  ++entry;
        hpnt[1] = *entry;  ++entry;
        hpnt[2] = *entry;  ++entry;
        hpnt[3] = *entry;
    }

    void SetColumn(int column, const THPoint<Real>& hpnt)
    {
        int i = column;
        Real* entry = &mEntry[i];
        *entry = hpnt[0];  entry += 4;
        *entry = hpnt[1];  entry += 4;
        *entry = hpnt[2];  entry += 4;
        *entry = hpnt[3];
    }

    void GetColumn(int column, THPoint<Real>& hpnt) const
    {
        int i = column;
        const Real* entry = &mEntry[i];
        hpnt[0] = *entry;  entry += 4;
        hpnt[1] = *entry;  entry += 4;
        hpnt[2] = *entry;  entry += 4;
        hpnt[3] = *entry;
    }

    // The matrix is stored in row-major order.  Store these in column-major
    // order in the specified array, which must have at least 16 slots.
    void GetColumnMajor(Real* columnMajor) const
    {
        columnMajor[0] = mEntry[0];
        columnMajor[1] = mEntry[4];
        columnMajor[2] = mEntry[8];
        columnMajor[3] = mEntry[12];
        columnMajor[4] = mEntry[1];
        columnMajor[5] = mEntry[5];
        columnMajor[6] = mEntry[9];
        columnMajor[7] = mEntry[13];
        columnMajor[8] = mEntry[2];
        columnMajor[9] = mEntry[6];
        columnMajor[10] = mEntry[10];
        columnMajor[11] = mEntry[14];
        columnMajor[12] = mEntry[3];
        columnMajor[13] = mEntry[7];
        columnMajor[14] = mEntry[11];
        columnMajor[15] = mEntry[15];
    }

    // Assignment.
    THMatrix& operator= (const THMatrix& mat)
    {
        mEntry[0] = mat.mEntry[0];
        mEntry[1] = mat.mEntry[1];
        mEntry[2] = mat.mEntry[2];
        mEntry[3] = mat.mEntry[3];
        mEntry[4] = mat.mEntry[4];
        mEntry[5] = mat.mEntry[5];
        mEntry[6] = mat.mEntry[6];
        mEntry[7] = mat.mEntry[7];
        mEntry[8] = mat.mEntry[8];
        mEntry[9] = mat.mEntry[9];
        mEntry[10] = mat.mEntry[10];
        mEntry[11] = mat.mEntry[11];
        mEntry[12] = mat.mEntry[12];
        mEntry[13] = mat.mEntry[13];
        mEntry[14] = mat.mEntry[14];
        mEntry[15] = mat.mEntry[15];
        return *this;
    }

    // Comparison (for use by STL containers).
    bool operator== (const THMatrix& mat) const
    {
        for (int i = 0; i < 16; ++i)
        {
            if (mEntry[i] != mat.mEntry[i])
            {
                return false;
            }
        }
        return true;
    }

    bool operator!= (const THMatrix& mat) const
    {
        return !operator==(mat);
    }

    bool operator< (const THMatrix& mat) const
    {
        // lexicographical ordering
        for (int i = 0; i < 16; ++i)
        {
            if (mEntry[i] < mat.mEntry[i])
            {
                return true;
            }
            if (mEntry[i] > mat.mEntry[i])
            {
                return false;
            }
        }
        return false;
    }

    bool operator<= (const THMatrix& mat) const
    {
        // (x <= y) <=> !(y < x)
        return !(mat.operator<(*this));
    }

    bool operator> (const THMatrix& mat) const
    {
        // (x > y) <=> (y < x)
        return mat.operator<(*this);
    }

    bool operator>= (const THMatrix& mat) const
    {
        // (x >= y) <=> !(x < y)
        return !operator<(mat);
    }

    // Arithmetic operations.
    THMatrix operator+ (const THMatrix& mat) const
    {
        THMatrix result;
        for (int i = 0; i < 16; ++i)
        {
            result.mEntry[i] = mEntry[i] + mat.mEntry[i];
        }
        return result;
    }

    THMatrix operator- (const THMatrix& mat) const
    {
        THMatrix result;
        for (int i = 0; i < 16; ++i)
        {
            result.mEntry[i] = mEntry[i] - mat.mEntry[i];
        }
        return result;
    }

    THMatrix operator* (Real scalar) const
    {
        THMatrix result;
        for (int i = 0; i < 16; ++i)
        {
            result.mEntry[i] = scalar*mEntry[i];
        }
        return result;
    }

    THMatrix operator/ (Real scalar) const
    {
        THMatrix result;
        int i;

        if (scalar != (Real)0)
        {
            Real invScalar = (Real)1 / scalar;
            for (i = 0; i < 16; ++i)
            {
                result.mEntry[i] = invScalar*mEntry[i];
            }
        }
        else
        {
            Real infinity = std::numeric_limits<Real>::infinity();
            for (i = 0; i < 16; ++i)
            {
                result.mEntry[i] = infinity;
            }
        }

        return result;
    }

    THMatrix operator- () const
    {
        THMatrix result;
        for (int i = 0; i < 16; ++i)
        {
            result.mEntry[i] = -mEntry[i];
        }
        return result;
    }

    // Arithmetic updates.
    THMatrix& operator+= (const THMatrix& mat)
    {
        for (int i = 0; i < 16; ++i)
        {
            mEntry[i] += mat.mEntry[i];
        }
        return *this;
    }

    THMatrix& operator-= (const THMatrix& mat)
    {
        for (int i = 0; i < 16; ++i)
        {
            mEntry[i] -= mat.mEntry[i];
        }
        return *this;
    }

    THMatrix& operator*= (Real scalar)
    {
        for (int i = 0; i < 16; ++i)
        {
            mEntry[i] *= scalar;
        }
        return *this;
    }

    THMatrix& operator/= (Real scalar)
    {
        for (int i = 0; i < 16; ++i)
        {
            mEntry[i] /= scalar;
        }
        return *this;
    }

    // Operations on matrices.
    void MakeZero()  // Z
    {
        mEntry[0] = (Real)0;
        mEntry[1] = (Real)0;
        mEntry[2] = (Real)0;
        mEntry[3] = (Real)0;
        mEntry[4] = (Real)0;
        mEntry[5] = (Real)0;
        mEntry[6] = (Real)0;
        mEntry[7] = (Real)0;
        mEntry[8] = (Real)0;
        mEntry[9] = (Real)0;
        mEntry[10] = (Real)0;
        mEntry[11] = (Real)0;
        mEntry[12] = (Real)0;
        mEntry[13] = (Real)0;
        mEntry[14] = (Real)0;
        mEntry[15] = (Real)0;
    }

    void MakeIdentity()  // I
    {
        mEntry[0] = (Real)1;
        mEntry[1] = (Real)0;
        mEntry[2] = (Real)0;
        mEntry[3] = (Real)0;
        mEntry[4] = (Real)0;
        mEntry[5] = (Real)1;
        mEntry[6] = (Real)0;
        mEntry[7] = (Real)0;
        mEntry[8] = (Real)0;
        mEntry[9] = (Real)0;
        mEntry[10] = (Real)1;
        mEntry[11] = (Real)0;
        mEntry[12] = (Real)0;
        mEntry[13] = (Real)0;
        mEntry[14] = (Real)0;
        mEntry[15] = (Real)1;
    }

    void MakeDiagonal(Real m00, Real m11, Real m22)  // D
    {
        mEntry[0] = m00;
        mEntry[1] = (Real)0;
        mEntry[2] = (Real)0;
        mEntry[3] = (Real)0;
        mEntry[4] = (Real)0;
        mEntry[5] = m11;
        mEntry[6] = (Real)0;
        mEntry[7] = (Real)0;
        mEntry[8] = (Real)0;
        mEntry[9] = (Real)0;
        mEntry[10] = m22;
        mEntry[11] = (Real)0;
        mEntry[12] = (Real)0;
        mEntry[13] = (Real)0;
        mEntry[14] = (Real)0;
        mEntry[15] = (Real)1;
    }

    void MakeRotation(const TAVector<Real>& axis, Real angle)  // R
    {
        Real cs = cos(angle);
        Real sn = sin(angle);
        Real oneMinusCos = (Real)1 - cs;
        Real x2 = axis[0] * axis[0];
        Real y2 = axis[1] * axis[1];
        Real z2 = axis[2] * axis[2];
        Real xym = axis[0] * axis[1] * oneMinusCos;
        Real xzm = axis[0] * axis[2] * oneMinusCos;
        Real yzm = axis[1] * axis[2] * oneMinusCos;
        Real xSin = axis[0] * sn;
        Real ySin = axis[1] * sn;
        Real zSin = axis[2] * sn;

        mEntry[0] = x2 * oneMinusCos + cs;
        mEntry[1] = xym - zSin;
        mEntry[2] = xzm + ySin;
        mEntry[3] = (Real)0;
        mEntry[4] = xym + zSin;
        mEntry[5] = y2 * oneMinusCos + cs;
        mEntry[6] = yzm - xSin;
        mEntry[7] = (Real)0;
        mEntry[8] = xzm - ySin;
        mEntry[9] = yzm + xSin;
        mEntry[10] = z2 * oneMinusCos + cs;
        mEntry[11] = (Real)0;
        mEntry[12] = (Real)0;
        mEntry[13] = (Real)0;
        mEntry[14] = (Real)0;
        mEntry[15] = (Real)1;
    }

    THMatrix Transpose() const  // M^T
    {
        THMatrix result;
        result.mEntry[0] = mEntry[0];
        result.mEntry[1] = mEntry[4];
        result.mEntry[2] = mEntry[8];
        result.mEntry[3] = mEntry[12];
        result.mEntry[4] = mEntry[1];
        result.mEntry[5] = mEntry[5];
        result.mEntry[6] = mEntry[9];
        result.mEntry[7] = mEntry[13];
        result.mEntry[8] = mEntry[2];
        result.mEntry[9] = mEntry[6];
        result.mEntry[10] = mEntry[10];
        result.mEntry[11] = mEntry[14];
        result.mEntry[12] = mEntry[3];
        result.mEntry[13] = mEntry[7];
        result.mEntry[14] = mEntry[11];
        result.mEntry[15] = mEntry[15];
        return result;
    }

    THMatrix Inverse(const Real epsilon = (Real)0) const  // M^{-1}
    {
        Real a0 = mEntry[0] * mEntry[5] - mEntry[1] * mEntry[4];
        Real a1 = mEntry[0] * mEntry[6] - mEntry[2] * mEntry[4];
        Real a2 = mEntry[0] * mEntry[7] - mEntry[3] * mEntry[4];
        Real a3 = mEntry[1] * mEntry[6] - mEntry[2] * mEntry[5];
        Real a4 = mEntry[1] * mEntry[7] - mEntry[3] * mEntry[5];
        Real a5 = mEntry[2] * mEntry[7] - mEntry[3] * mEntry[6];
        Real b0 = mEntry[8] * mEntry[13] - mEntry[9] * mEntry[12];
        Real b1 = mEntry[8] * mEntry[14] - mEntry[10] * mEntry[12];
        Real b2 = mEntry[8] * mEntry[15] - mEntry[11] * mEntry[12];
        Real b3 = mEntry[9] * mEntry[14] - mEntry[10] * mEntry[13];
        Real b4 = mEntry[9] * mEntry[15] - mEntry[11] * mEntry[13];
        Real b5 = mEntry[10] * mEntry[15] - mEntry[11] * mEntry[14];

        Real det = a0*b5 - a1*b4 + a2*b3 + a3*b2 - a4*b1 + a5*b0;
        if (fabs(det) <= epsilon)
        {
            return ZERO();
        }

        THMatrix inverse;
        inverse.mEntry[0] = +mEntry[5] * b5 - mEntry[6] * b4 + mEntry[7] * b3;
        inverse.mEntry[4] = -mEntry[4] * b5 + mEntry[6] * b2 - mEntry[7] * b1;
        inverse.mEntry[8] = +mEntry[4] * b4 - mEntry[5] * b2 + mEntry[7] * b0;
        inverse.mEntry[12] = -mEntry[4] * b3 + mEntry[5] * b1 - mEntry[6] * b0;
        inverse.mEntry[1] = -mEntry[1] * b5 + mEntry[2] * b4 - mEntry[3] * b3;
        inverse.mEntry[5] = +mEntry[0] * b5 - mEntry[2] * b2 + mEntry[3] * b1;
        inverse.mEntry[9] = -mEntry[0] * b4 + mEntry[1] * b2 - mEntry[3] * b0;
        inverse.mEntry[13] = +mEntry[0] * b3 - mEntry[1] * b1 + mEntry[2] * b0;
        inverse.mEntry[2] = +mEntry[13] * a5 - mEntry[14] * a4 + mEntry[15] * a3;
        inverse.mEntry[6] = -mEntry[12] * a5 + mEntry[14] * a2 - mEntry[15] * a1;
        inverse.mEntry[10] = +mEntry[12] * a4 - mEntry[13] * a2 + mEntry[15] * a0;
        inverse.mEntry[14] = -mEntry[12] * a3 + mEntry[13] * a1 - mEntry[14] * a0;
        inverse.mEntry[3] = -mEntry[9] * a5 + mEntry[10] * a4 - mEntry[11] * a3;
        inverse.mEntry[7] = +mEntry[8] * a5 - mEntry[10] * a2 + mEntry[11] * a1;
        inverse.mEntry[11] = -mEntry[8] * a4 + mEntry[9] * a2 - mEntry[11] * a0;
        inverse.mEntry[15] = +mEntry[8] * a3 - mEntry[9] * a1 + mEntry[10] * a0;

        Real invDet = (Real)1 / det;
        inverse.mEntry[0] *= invDet;
        inverse.mEntry[1] *= invDet;
        inverse.mEntry[2] *= invDet;
        inverse.mEntry[3] *= invDet;
        inverse.mEntry[4] *= invDet;
        inverse.mEntry[5] *= invDet;
        inverse.mEntry[6] *= invDet;
        inverse.mEntry[7] *= invDet;
        inverse.mEntry[8] *= invDet;
        inverse.mEntry[9] *= invDet;
        inverse.mEntry[10] *= invDet;
        inverse.mEntry[11] *= invDet;
        inverse.mEntry[12] *= invDet;
        inverse.mEntry[13] *= invDet;
        inverse.mEntry[14] *= invDet;
        inverse.mEntry[15] *= invDet;

        return inverse;
    }

    THMatrix Adjoint() const  // M^{adj}
    {
        Real a0 = mEntry[0] * mEntry[5] - mEntry[1] * mEntry[4];
        Real a1 = mEntry[0] * mEntry[6] - mEntry[2] * mEntry[4];
        Real a2 = mEntry[0] * mEntry[7] - mEntry[3] * mEntry[4];
        Real a3 = mEntry[1] * mEntry[6] - mEntry[2] * mEntry[5];
        Real a4 = mEntry[1] * mEntry[7] - mEntry[3] * mEntry[5];
        Real a5 = mEntry[2] * mEntry[7] - mEntry[3] * mEntry[6];
        Real b0 = mEntry[8] * mEntry[13] - mEntry[9] * mEntry[12];
        Real b1 = mEntry[8] * mEntry[14] - mEntry[10] * mEntry[12];
        Real b2 = mEntry[8] * mEntry[15] - mEntry[11] * mEntry[12];
        Real b3 = mEntry[9] * mEntry[14] - mEntry[10] * mEntry[13];
        Real b4 = mEntry[9] * mEntry[15] - mEntry[11] * mEntry[13];
        Real b5 = mEntry[10] * mEntry[15] - mEntry[11] * mEntry[14];

        return THMatrix(
            +mEntry[5] * b5 - mEntry[6] * b4 + mEntry[7] * b3,
            -mEntry[1] * b5 + mEntry[2] * b4 - mEntry[3] * b3,
            +mEntry[13] * a5 - mEntry[14] * a4 + mEntry[15] * a3,
            -mEntry[9] * a5 + mEntry[10] * a4 - mEntry[11] * a3,
            -mEntry[4] * b5 + mEntry[6] * b2 - mEntry[7] * b1,
            +mEntry[0] * b5 - mEntry[2] * b2 + mEntry[3] * b1,
            -mEntry[12] * a5 + mEntry[14] * a2 - mEntry[15] * a1,
            +mEntry[8] * a5 - mEntry[10] * a2 + mEntry[11] * a1,
            +mEntry[4] * b4 - mEntry[5] * b2 + mEntry[7] * b0,
            -mEntry[0] * b4 + mEntry[1] * b2 - mEntry[3] * b0,
            +mEntry[12] * a4 - mEntry[13] * a2 + mEntry[15] * a0,
            -mEntry[8] * a4 + mEntry[9] * a2 - mEntry[11] * a0,
            -mEntry[4] * b3 + mEntry[5] * b1 - mEntry[6] * b0,
            +mEntry[0] * b3 - mEntry[1] * b1 + mEntry[2] * b0,
            -mEntry[12] * a3 + mEntry[13] * a1 - mEntry[14] * a0,
            +mEntry[8] * a3 - mEntry[9] * a1 + mEntry[10] * a0);
    }

    Real Determinant() const  // det(M)
    {
        Real a0 = mEntry[0] * mEntry[5] - mEntry[1] * mEntry[4];
        Real a1 = mEntry[0] * mEntry[6] - mEntry[2] * mEntry[4];
        Real a2 = mEntry[0] * mEntry[7] - mEntry[3] * mEntry[4];
        Real a3 = mEntry[1] * mEntry[6] - mEntry[2] * mEntry[5];
        Real a4 = mEntry[1] * mEntry[7] - mEntry[3] * mEntry[5];
        Real a5 = mEntry[2] * mEntry[7] - mEntry[3] * mEntry[6];
        Real b0 = mEntry[8] * mEntry[13] - mEntry[9] * mEntry[12];
        Real b1 = mEntry[8] * mEntry[14] - mEntry[10] * mEntry[12];
        Real b2 = mEntry[8] * mEntry[15] - mEntry[11] * mEntry[12];
        Real b3 = mEntry[9] * mEntry[14] - mEntry[10] * mEntry[13];
        Real b4 = mEntry[9] * mEntry[15] - mEntry[11] * mEntry[13];
        Real b5 = mEntry[10] * mEntry[15] - mEntry[11] * mEntry[14];
        Real det = a0 * b5 - a1 * b4 + a2 * b3 + a3 * b2 - a4 * b1 + a5 * b0;
        return det;
    }

    THMatrix operator* (const THMatrix& mat) const  // M*mat
    {
        return THMatrix(
            mEntry[0] * mat.mEntry[0] +
            mEntry[1] * mat.mEntry[4] +
            mEntry[2] * mat.mEntry[8] +
            mEntry[3] * mat.mEntry[12],

            mEntry[0] * mat.mEntry[1] +
            mEntry[1] * mat.mEntry[5] +
            mEntry[2] * mat.mEntry[9] +
            mEntry[3] * mat.mEntry[13],

            mEntry[0] * mat.mEntry[2] +
            mEntry[1] * mat.mEntry[6] +
            mEntry[2] * mat.mEntry[10] +
            mEntry[3] * mat.mEntry[14],

            mEntry[0] * mat.mEntry[3] +
            mEntry[1] * mat.mEntry[7] +
            mEntry[2] * mat.mEntry[11] +
            mEntry[3] * mat.mEntry[15],

            mEntry[4] * mat.mEntry[0] +
            mEntry[5] * mat.mEntry[4] +
            mEntry[6] * mat.mEntry[8] +
            mEntry[7] * mat.mEntry[12],

            mEntry[4] * mat.mEntry[1] +
            mEntry[5] * mat.mEntry[5] +
            mEntry[6] * mat.mEntry[9] +
            mEntry[7] * mat.mEntry[13],

            mEntry[4] * mat.mEntry[2] +
            mEntry[5] * mat.mEntry[6] +
            mEntry[6] * mat.mEntry[10] +
            mEntry[7] * mat.mEntry[14],

            mEntry[4] * mat.mEntry[3] +
            mEntry[5] * mat.mEntry[7] +
            mEntry[6] * mat.mEntry[11] +
            mEntry[7] * mat.mEntry[15],

            mEntry[8] * mat.mEntry[0] +
            mEntry[9] * mat.mEntry[4] +
            mEntry[10] * mat.mEntry[8] +
            mEntry[11] * mat.mEntry[12],

            mEntry[8] * mat.mEntry[1] +
            mEntry[9] * mat.mEntry[5] +
            mEntry[10] * mat.mEntry[9] +
            mEntry[11] * mat.mEntry[13],

            mEntry[8] * mat.mEntry[2] +
            mEntry[9] * mat.mEntry[6] +
            mEntry[10] * mat.mEntry[10] +
            mEntry[11] * mat.mEntry[14],

            mEntry[8] * mat.mEntry[3] +
            mEntry[9] * mat.mEntry[7] +
            mEntry[10] * mat.mEntry[11] +
            mEntry[11] * mat.mEntry[15],

            mEntry[12] * mat.mEntry[0] +
            mEntry[13] * mat.mEntry[4] +
            mEntry[14] * mat.mEntry[8] +
            mEntry[15] * mat.mEntry[12],

            mEntry[12] * mat.mEntry[1] +
            mEntry[13] * mat.mEntry[5] +
            mEntry[14] * mat.mEntry[9] +
            mEntry[15] * mat.mEntry[13],

            mEntry[12] * mat.mEntry[2] +
            mEntry[13] * mat.mEntry[6] +
            mEntry[14] * mat.mEntry[10] +
            mEntry[15] * mat.mEntry[14],

            mEntry[12] * mat.mEntry[3] +
            mEntry[13] * mat.mEntry[7] +
            mEntry[14] * mat.mEntry[11] +
            mEntry[15] * mat.mEntry[15]);
    }

    THMatrix TransposeTimes(const THMatrix& mat) const  // M^T*mat
    {
        // P = A^T*B
        return THMatrix(
            mEntry[0] * mat.mEntry[0] +
            mEntry[4] * mat.mEntry[4] +
            mEntry[8] * mat.mEntry[8] +
            mEntry[12] * mat.mEntry[12],

            mEntry[0] * mat.mEntry[1] +
            mEntry[4] * mat.mEntry[5] +
            mEntry[8] * mat.mEntry[9] +
            mEntry[12] * mat.mEntry[13],

            mEntry[0] * mat.mEntry[2] +
            mEntry[4] * mat.mEntry[6] +
            mEntry[8] * mat.mEntry[10] +
            mEntry[12] * mat.mEntry[14],

            mEntry[0] * mat.mEntry[3] +
            mEntry[4] * mat.mEntry[7] +
            mEntry[8] * mat.mEntry[11] +
            mEntry[12] * mat.mEntry[15],

            mEntry[1] * mat.mEntry[0] +
            mEntry[5] * mat.mEntry[4] +
            mEntry[9] * mat.mEntry[8] +
            mEntry[13] * mat.mEntry[12],

            mEntry[1] * mat.mEntry[1] +
            mEntry[5] * mat.mEntry[5] +
            mEntry[9] * mat.mEntry[9] +
            mEntry[13] * mat.mEntry[13],

            mEntry[1] * mat.mEntry[2] +
            mEntry[5] * mat.mEntry[6] +
            mEntry[9] * mat.mEntry[10] +
            mEntry[13] * mat.mEntry[14],

            mEntry[1] * mat.mEntry[3] +
            mEntry[5] * mat.mEntry[7] +
            mEntry[9] * mat.mEntry[11] +
            mEntry[13] * mat.mEntry[15],

            mEntry[2] * mat.mEntry[0] +
            mEntry[6] * mat.mEntry[4] +
            mEntry[10] * mat.mEntry[8] +
            mEntry[14] * mat.mEntry[12],

            mEntry[2] * mat.mEntry[1] +
            mEntry[6] * mat.mEntry[5] +
            mEntry[10] * mat.mEntry[9] +
            mEntry[14] * mat.mEntry[13],

            mEntry[2] * mat.mEntry[2] +
            mEntry[6] * mat.mEntry[6] +
            mEntry[10] * mat.mEntry[10] +
            mEntry[14] * mat.mEntry[14],

            mEntry[2] * mat.mEntry[3] +
            mEntry[6] * mat.mEntry[7] +
            mEntry[10] * mat.mEntry[11] +
            mEntry[14] * mat.mEntry[15],

            mEntry[3] * mat.mEntry[0] +
            mEntry[7] * mat.mEntry[4] +
            mEntry[11] * mat.mEntry[8] +
            mEntry[15] * mat.mEntry[12],

            mEntry[3] * mat.mEntry[1] +
            mEntry[7] * mat.mEntry[5] +
            mEntry[11] * mat.mEntry[9] +
            mEntry[15] * mat.mEntry[13],

            mEntry[3] * mat.mEntry[2] +
            mEntry[7] * mat.mEntry[6] +
            mEntry[11] * mat.mEntry[10] +
            mEntry[15] * mat.mEntry[14],

            mEntry[3] * mat.mEntry[3] +
            mEntry[7] * mat.mEntry[7] +
            mEntry[11] * mat.mEntry[11] +
            mEntry[15] * mat.mEntry[15]);
    }

    THMatrix TimesTranspose(const THMatrix& mat) const  // M*mat^T
    {
        // P = A*B^T
        return THMatrix(
            mEntry[0] * mat.mEntry[0] +
            mEntry[1] * mat.mEntry[1] +
            mEntry[2] * mat.mEntry[2] +
            mEntry[3] * mat.mEntry[3],

            mEntry[0] * mat.mEntry[4] +
            mEntry[1] * mat.mEntry[5] +
            mEntry[2] * mat.mEntry[6] +
            mEntry[3] * mat.mEntry[7],

            mEntry[0] * mat.mEntry[8] +
            mEntry[1] * mat.mEntry[9] +
            mEntry[2] * mat.mEntry[10] +
            mEntry[3] * mat.mEntry[11],

            mEntry[0] * mat.mEntry[12] +
            mEntry[1] * mat.mEntry[13] +
            mEntry[2] * mat.mEntry[14] +
            mEntry[3] * mat.mEntry[15],

            mEntry[4] * mat.mEntry[0] +
            mEntry[5] * mat.mEntry[1] +
            mEntry[6] * mat.mEntry[2] +
            mEntry[7] * mat.mEntry[3],

            mEntry[4] * mat.mEntry[4] +
            mEntry[5] * mat.mEntry[5] +
            mEntry[6] * mat.mEntry[6] +
            mEntry[7] * mat.mEntry[7],

            mEntry[4] * mat.mEntry[8] +
            mEntry[5] * mat.mEntry[9] +
            mEntry[6] * mat.mEntry[10] +
            mEntry[7] * mat.mEntry[11],

            mEntry[4] * mat.mEntry[12] +
            mEntry[5] * mat.mEntry[13] +
            mEntry[6] * mat.mEntry[14] +
            mEntry[7] * mat.mEntry[15],

            mEntry[8] * mat.mEntry[0] +
            mEntry[9] * mat.mEntry[1] +
            mEntry[10] * mat.mEntry[2] +
            mEntry[11] * mat.mEntry[3],

            mEntry[8] * mat.mEntry[4] +
            mEntry[9] * mat.mEntry[5] +
            mEntry[10] * mat.mEntry[6] +
            mEntry[11] * mat.mEntry[7],

            mEntry[8] * mat.mEntry[8] +
            mEntry[9] * mat.mEntry[9] +
            mEntry[10] * mat.mEntry[10] +
            mEntry[11] * mat.mEntry[11],

            mEntry[8] * mat.mEntry[12] +
            mEntry[9] * mat.mEntry[13] +
            mEntry[10] * mat.mEntry[14] +
            mEntry[11] * mat.mEntry[15],

            mEntry[12] * mat.mEntry[0] +
            mEntry[13] * mat.mEntry[1] +
            mEntry[14] * mat.mEntry[2] +
            mEntry[15] * mat.mEntry[3],

            mEntry[12] * mat.mEntry[4] +
            mEntry[13] * mat.mEntry[5] +
            mEntry[14] * mat.mEntry[6] +
            mEntry[15] * mat.mEntry[7],

            mEntry[12] * mat.mEntry[8] +
            mEntry[13] * mat.mEntry[9] +
            mEntry[14] * mat.mEntry[10] +
            mEntry[15] * mat.mEntry[11],

            mEntry[12] * mat.mEntry[12] +
            mEntry[13] * mat.mEntry[13] +
            mEntry[14] * mat.mEntry[14] +
            mEntry[15] * mat.mEntry[15]);
    }

    THMatrix TransposeTimesTranspose(const THMatrix& mat) const  // M^T*mat^T
    {
        // P = A^T*B^T
        return THMatrix(
            mEntry[0] * mat.mEntry[0] +
            mEntry[4] * mat.mEntry[1] +
            mEntry[8] * mat.mEntry[2] +
            mEntry[12] * mat.mEntry[3],

            mEntry[0] * mat.mEntry[4] +
            mEntry[4] * mat.mEntry[5] +
            mEntry[8] * mat.mEntry[6] +
            mEntry[12] * mat.mEntry[7],

            mEntry[0] * mat.mEntry[8] +
            mEntry[4] * mat.mEntry[9] +
            mEntry[8] * mat.mEntry[10] +
            mEntry[12] * mat.mEntry[11],

            mEntry[0] * mat.mEntry[12] +
            mEntry[4] * mat.mEntry[13] +
            mEntry[8] * mat.mEntry[14] +
            mEntry[12] * mat.mEntry[15],

            mEntry[1] * mat.mEntry[0] +
            mEntry[5] * mat.mEntry[1] +
            mEntry[9] * mat.mEntry[2] +
            mEntry[13] * mat.mEntry[3],

            mEntry[1] * mat.mEntry[4] +
            mEntry[5] * mat.mEntry[5] +
            mEntry[9] * mat.mEntry[6] +
            mEntry[13] * mat.mEntry[7],

            mEntry[1] * mat.mEntry[8] +
            mEntry[5] * mat.mEntry[9] +
            mEntry[9] * mat.mEntry[10] +
            mEntry[13] * mat.mEntry[11],

            mEntry[1] * mat.mEntry[12] +
            mEntry[5] * mat.mEntry[13] +
            mEntry[9] * mat.mEntry[14] +
            mEntry[13] * mat.mEntry[15],

            mEntry[2] * mat.mEntry[0] +
            mEntry[6] * mat.mEntry[1] +
            mEntry[10] * mat.mEntry[2] +
            mEntry[14] * mat.mEntry[3],

            mEntry[2] * mat.mEntry[4] +
            mEntry[6] * mat.mEntry[5] +
            mEntry[10] * mat.mEntry[6] +
            mEntry[14] * mat.mEntry[7],

            mEntry[2] * mat.mEntry[8] +
            mEntry[6] * mat.mEntry[9] +
            mEntry[10] * mat.mEntry[10] +
            mEntry[14] * mat.mEntry[11],

            mEntry[2] * mat.mEntry[12] +
            mEntry[6] * mat.mEntry[13] +
            mEntry[10] * mat.mEntry[14] +
            mEntry[14] * mat.mEntry[15],

            mEntry[3] * mat.mEntry[0] +
            mEntry[7] * mat.mEntry[1] +
            mEntry[11] * mat.mEntry[2] +
            mEntry[15] * mat.mEntry[3],

            mEntry[3] * mat.mEntry[4] +
            mEntry[7] * mat.mEntry[5] +
            mEntry[11] * mat.mEntry[6] +
            mEntry[15] * mat.mEntry[7],

            mEntry[3] * mat.mEntry[8] +
            mEntry[7] * mat.mEntry[9] +
            mEntry[11] * mat.mEntry[10] +
            mEntry[15] * mat.mEntry[11],

            mEntry[3] * mat.mEntry[12] +
            mEntry[7] * mat.mEntry[13] +
            mEntry[11] * mat.mEntry[14] +
            mEntry[15] * mat.mEntry[15]);
    }

    THMatrix TimesDiagonal(const TAPoint<Real>& diag) const  // M*D
    {
        return THMatrix(
            mEntry[0] * diag[0], mEntry[1] * diag[1], mEntry[2] * diag[2],
            mEntry[3],
            mEntry[4] * diag[0], mEntry[5] * diag[1], mEntry[6] * diag[2],
            mEntry[7],
            mEntry[8] * diag[0], mEntry[9] * diag[1], mEntry[10] * diag[2],
            mEntry[11],
            mEntry[12] * diag[0], mEntry[13] * diag[1], mEntry[14] * diag[2],
            mEntry[15]
        );
    }

    THMatrix DiagonalTimes(const TAPoint<Real>& diag) const  // D*M
    {
        return THMatrix(
            diag[0] * mEntry[0], diag[0] * mEntry[1], diag[0] * mEntry[2],
            mEntry[3],
            diag[1] * mEntry[4], diag[1] * mEntry[5], diag[1] * mEntry[6],
            mEntry[7],
            diag[2] * mEntry[8], diag[2] * mEntry[9], diag[2] * mEntry[10],
            mEntry[11],
            mEntry[12], mEntry[13], mEntry[14], mEntry[15]
        );
    }

    // This operation applies to 3x3 upper-left block of M.
    void Orthonormalize()
    {
        // Algorithm uses Gram-Schmidt orthogonalization.  If 'this' matrix has
        // upper-left 3x3 block M = [m0|m1|m2], then the orthonormal output matrix
        // is Q = [q0|q1|q2],
        //
        //   q0 = m0/|m0|
        //   q1 = (m1-(q0*m1)q0)/|m1-(q0*m1)q0|
        //   q2 = (m2-(q0*m2)q0-(q1*m2)q1)/|m2-(q0*m2)q0-(q1*m2)q1|
        //
        // where |V| indicates length of vector V and A*B indicates dot
        // product of vectors A and B.

        // Compute q0.
        Real invLength = (Real)1 / sqrt(mEntry[0] * mEntry[0] +
            mEntry[4] * mEntry[4] + mEntry[8] * mEntry[8]);

        mEntry[0] *= invLength;
        mEntry[4] *= invLength;
        mEntry[8] *= invLength;

        // Compute q1.
        Real dot0 = mEntry[0] * mEntry[1] + mEntry[4] * mEntry[5] +
            mEntry[8] * mEntry[9];

        mEntry[1] -= dot0*mEntry[0];
        mEntry[5] -= dot0*mEntry[4];
        mEntry[9] -= dot0*mEntry[8];

        invLength = (Real)1 / sqrt(mEntry[1] * mEntry[1] +
            mEntry[5] * mEntry[5] + mEntry[9] * mEntry[9]);

        mEntry[1] *= invLength;
        mEntry[5] *= invLength;
        mEntry[9] *= invLength;

        // Compute q2.
        Real dot1 = mEntry[1] * mEntry[2] + mEntry[5] * mEntry[6] +
            mEntry[9] * mEntry[10];

        dot0 = mEntry[0] * mEntry[2] + mEntry[4] * mEntry[6] +
            mEntry[8] * mEntry[10];

        mEntry[2] -= dot0*mEntry[0] + dot1*mEntry[1];
        mEntry[6] -= dot0*mEntry[4] + dot1*mEntry[5];
        mEntry[10] -= dot0*mEntry[8] + dot1*mEntry[9];

        invLength = (Real)1 / sqrt(mEntry[2] * mEntry[2] +
            mEntry[6] * mEntry[6] + mEntry[10] * mEntry[10]);

        mEntry[2] *= invLength;
        mEntry[6] *= invLength;
        mEntry[10] *= invLength;
    }

    // Operations between matrices and homogeneous points.  Both M and p
    // are general homogeneous objects (M not required to be affine, p
    // not required to have w = 1).
    THPoint<Real> operator* (const THPoint<Real>& p) const  // M*p
    {
        return THPoint<Real>(
            mEntry[0] * p[0] +
            mEntry[1] * p[1] +
            mEntry[2] * p[2] +
            mEntry[3] * p[3],

            mEntry[4] * p[0] +
            mEntry[5] * p[1] +
            mEntry[6] * p[2] +
            mEntry[7] * p[3],

            mEntry[8] * p[0] +
            mEntry[9] * p[1] +
            mEntry[10] * p[2] +
            mEntry[11] * p[3],

            mEntry[12] * p[0] +
            mEntry[13] * p[1] +
            mEntry[14] * p[2] +
            mEntry[15] * p[3]);
    }

    void BatchMultiply(int numPoints, const THPoint<Real>* input,
        THPoint<Real>* output) const  // M*p[0], ..., M*p[n-1]
    {
        const THPoint<Real>* src = input;
        THPoint<Real>* trg = output;
        for (int i = 0; i < numPoints; ++i, src++, trg++)
        {
            *trg = (*this)*(*src);
        }
    }

    // Operations between affine matrices and affine points.
    TAPoint<Real> operator* (const TAPoint<Real>& p) const  // M*p
    {
        return TAPoint<Real>(
            mEntry[0] * p[0] +
            mEntry[1] * p[1] +
            mEntry[2] * p[2] +
            mEntry[3],

            mEntry[4] * p[0] +
            mEntry[5] * p[1] +
            mEntry[6] * p[2] +
            mEntry[7],

            mEntry[8] * p[0] +
            mEntry[9] * p[1] +
            mEntry[10] * p[2] +
            mEntry[11]);
    }

    void BatchMultiply(int numPoints, const TAPoint<Real>* input,
        TAPoint<Real>* output) const  // M*p[0], ..., M*p[n-1]
    {
        const TAPoint<Real>* src = input;
        TAPoint<Real>* trg = output;
        for (int i = 0; i < numPoints; ++i, src++, trg++)
        {
            *trg = (*this)*(*src);
        }
    }

    // Operations between affine matrices and affine vectors.
    TAVector<Real> operator* (const TAVector<Real>& p) const  // M*v
    {
        return TAVector<Real>(
            mEntry[0] * p[0] +
            mEntry[1] * p[1] +
            mEntry[2] * p[2],

            mEntry[4] * p[0] +
            mEntry[5] * p[1] +
            mEntry[6] * p[2],

            mEntry[8] * p[0] +
            mEntry[9] * p[1] +
            mEntry[10] * p[2]);
    }

    void BatchMultiply(int numPoints, const TAVector<Real>* input,
        TAVector<Real>* output) const  // M*v[0], ..., M*v[n-1]
    {
        const TAVector<Real>* src = input;
        TAVector<Real>* trg = output;
        for (int i = 0; i < numPoints; ++i, src++, trg++)
        {
            *trg = (*this)*(*src);
        }
    }

    // Compute a quadratic forms.
    Real QForm(const THPoint<Real>& p0, const THPoint<Real>& p1) const  // p0^T*M*p1
    {
        THPoint<Real> Mp1(
            mEntry[0] * p1[0] +
            mEntry[1] * p1[1] +
            mEntry[2] * p1[2] +
            mEntry[3] * p1[3],

            mEntry[4] * p1[0] +
            mEntry[5] * p1[1] +
            mEntry[6] * p1[2] +
            mEntry[7] * p1[3],

            mEntry[8] * p1[0] +
            mEntry[9] * p1[1] +
            mEntry[10] * p1[2] +
            mEntry[11] * p1[3],

            mEntry[12] * p1[0] +
            mEntry[13] * p1[1] +
            mEntry[14] * p1[2] +
            mEntry[15] * p1[3]);

        Real dot = p0[0] * Mp1[0] + p0[1] * Mp1[1] + p0[2] * Mp1[2] + p0[3] * Mp1[3];
        return dot;
    }

    // Set the transformation to an oblique projection matrix onto a
    // specified plane.  The plane has an 'origin' point and a unit-length
    // 'normal'.
    void MakeObliqueProjection(const TAPoint<Real>& origin,
        const TAVector<Real>& normal, const TAVector<Real>& direction)
    {
        // The projection plane is Dot(N,X-P) = 0 where N is a 3-by-1 unit-length
        // normal vector and P is a 3-by-1 point on the plane.  The projection
        // is oblique to the plane, in the direction of the 3-by-1 vector D.
        // Necessarily Dot(N,D) is not zero for this projection to make sense.
        // Given a 3-by-1 point U, compute the intersection of the line U+t*D
        // with the plane to obtain t = -Dot(N,U-P)/Dot(N,D).  Then
        //
        //   projection(U) = P + [I - D*N^T/Dot(N,D)]*(U-P)
        //
        // A 4-by-4 homogeneous transformation representing the projection is
        //
        //       +-                               -+
        //   M = | D*N^T - Dot(N,D)*I   -Dot(N,P)D |
        //       |          0^T          -Dot(N,D) |
        //       +-                               -+
        //
        // where M applies to [U^T 1]^T by M*[U^T 1]^T.  The matrix is chosen so
        // that M[3][3] > 0 whenever Dot(N,D) < 0 (projection is onto the
        // "positive side" of the plane).

        Real dotND = normal.Dot(direction);
        Real dotNO = origin.Dot(normal);

        mEntry[0] = direction[0] * normal[0] - dotND;
        mEntry[1] = direction[0] * normal[1];
        mEntry[2] = direction[0] * normal[2];
        mEntry[3] = -dotNO*direction[0];
        mEntry[4] = direction[1] * normal[0];
        mEntry[5] = direction[1] * normal[1] - dotND;
        mEntry[6] = direction[1] * normal[2];
        mEntry[7] = -dotNO*direction[1];
        mEntry[8] = direction[2] * normal[0];
        mEntry[9] = direction[2] * normal[1];
        mEntry[10] = direction[2] * normal[2] - dotND;
        mEntry[11] = -dotNO*direction[2];
        mEntry[12] = (Real)0;
        mEntry[13] = (Real)0;
        mEntry[14] = (Real)0;
        mEntry[15] = -dotND;
    }

    // Set the transformation to a perspective projection matrix onto a
    // specified plane.  The plane has an 'origin' point and a unit-length
    // 'normal'.  The 'eye' is the origin of projection.
    void MakePerspectiveProjection(const TAPoint<Real>& origin,
        const TAVector<Real>& normal, const TAPoint<Real>& eye)
    {
        //     +-                                                 -+
        // M = | Dot(N,E-P)*I - E*N^T    -(Dot(N,E-P)*I - E*N^T)*E |
        //     |        -N^t                      Dot(N,E)         |
        //     +-                                                 -+
        //
        // where E is the eye point, P is a point on the plane, and N is a
        // unit-length plane normal.

        Real dotND = normal.Dot(eye - origin);

        mEntry[0] = dotND - eye[0] * normal[0];
        mEntry[1] = -eye[0] * normal[1];
        mEntry[2] = -eye[0] * normal[2];
        mEntry[3] = -(mEntry[0] * eye[0] + mEntry[1] * eye[1] + mEntry[2] * eye[2]);
        mEntry[4] = -eye[1] * normal[0];
        mEntry[5] = dotND - eye[1] * normal[1];
        mEntry[6] = -eye[1] * normal[2];
        mEntry[7] = -(mEntry[4] * eye[0] + mEntry[5] * eye[1] + mEntry[6] * eye[2]);
        mEntry[8] = -eye[2] * normal[0];
        mEntry[9] = -eye[2] * normal[1];
        mEntry[10] = dotND - eye[2] * normal[2];
        mEntry[11] = -(mEntry[8] * eye[0] + mEntry[9] * eye[1] + mEntry[10] * eye[2]);
        mEntry[12] = -normal[0];
        mEntry[13] = -normal[1];
        mEntry[14] = -normal[2];
        mEntry[15] = eye.Dot(normal);
    }

    // Set the transformation to a reflection matrix through a specified
    // plane.  The plane has an 'origin' point and a unit-length 'normal'.
    void MakeReflection(const TAPoint<Real>& origin,
        const TAVector<Real>& normal)
    {
        //     +-                         -+
        // M = | I-2*N*N^T    2*Dot(N,P)*N |
        //     |     0^T            1      |
        //     +-                         -+
        //
        // where P is a point on the plane and N is a unit-length plane normal.

        Real twoDotNO = ((Real)2) * origin.Dot(normal);

        mEntry[0] = (Real)1 - ((Real)2) * normal[0] * normal[0];
        mEntry[1] = ((Real)-2) * normal[0] * normal[1];
        mEntry[2] = ((Real)-2) * normal[0] * normal[2];
        mEntry[3] = twoDotNO * normal[0];
        mEntry[4] = ((Real)-2) * normal[1] * normal[0];
        mEntry[5] = (Real)1 - 2.0f*normal[1] * normal[1];
        mEntry[6] = ((Real)-2) * normal[1] * normal[2];
        mEntry[7] = twoDotNO * normal[1];
        mEntry[8] = ((Real)-2) * normal[2] * normal[0];
        mEntry[9] = ((Real)-2) * normal[2] * normal[1];
        mEntry[10] = (Real)1 - ((Real)2) * normal[2] * normal[2];
        mEntry[11] = twoDotNO * normal[2];
        mEntry[12] = (Real)0;
        mEntry[13] = (Real)0;
        mEntry[14] = (Real)0;
        mEntry[15] = (Real)1;
    }

    // Special matrices.
    static const THMatrix ZERO()
    {
        return THMatrix(true);
    }

    static const THMatrix IDENTITY()
    {
        return THMatrix(false);
    }

private:
    // The matrix is stored in row-major order.
    Real mEntry[16];
};

template <typename Real>
THMatrix<Real> operator* (Real scalar, const THMatrix<Real>& mat)
{
    return mat*scalar;
}

template <typename Real>
THPoint<Real> operator* (const THPoint<Real>& p, const THMatrix<Real>& mat)  // p*M
{
    return THPoint<Real>(
        p[0] * mat[0][0] + p[1] * mat[1][0] + p[2] * mat[2][0] + p[3] * mat[3][0],
        p[0] * mat[0][1] + p[1] * mat[1][1] + p[2] * mat[2][1] + p[3] * mat[3][1],
        p[0] * mat[0][2] + p[1] * mat[1][2] + p[2] * mat[2][2] + p[3] * mat[3][2],
        p[0] * mat[0][3] + p[1] * mat[1][3] + p[2] * mat[2][3] + p[3] * mat[3][3]);
}

}

#endif
