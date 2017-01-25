// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <GTEngineDEF.h>
#include <cmath>
#include <cstdint>
#include <xmmintrin.h>
#include <emmintrin.h>

namespace gte
{
// Support for Intel's Streaming SIMD Extensions (SSE) using 128-bit registers
// that store four 32-bit floating-point numbers.
class GTE_IMPEXP SIMD
{
public:
    // The representation of the SIMD 4-tuple.
    class GTE_IMPEXP Vector
    {
    public:
        // Information about vectors.
        enum { NUM_ELEMENTS = 4 };
        typedef float ElementType;

        // Construction.
        Vector ();
        Vector (Vector const& vec);
        Vector (__m128 const vec);
        Vector (__m128i const vec);
        Vector (float number);
        Vector (float n0, float n1, float n2, float n3);
        Vector (uint32_t encoding);
        Vector (uint32_t e0, uint32_t e1, uint32_t e2, uint32_t e3);

        // Assignment.
        Vector& operator= (Vector const& vec);
        Vector& operator= (__m128 const vec);
        Vector& operator= (__m128i const vec);

        // Implicit conversions.
        operator __m128 ();
        operator __m128 () const;
        operator __m128i ();
        operator __m128i () const;

    protected:
        __m128 mTuple;
    };

    // The representation of the SIMD 4x4-table.
    class GTE_IMPEXP Matrix
    {
    public:
        // Information about matrices.
        enum
        {
            NUM_ROWS = 4,
            NUM_COLS = 4,
            NUM_ELEMENTS = 16,
#if defined(GTE_USE_ROW_MAJOR)
            STORAGE_ROW_MAJOR = 1,
#else
            STORAGE_ROW_MAJOR = 0,
#endif
        };
        typedef float ElementType;

        // Construction.
        Matrix ();
        Matrix (Matrix const& mat);
        Matrix (__m128 const* mat);
        Matrix (
            float m00, float m01, float m02, float m03,
            float m10, float m11, float m12, float m13,
            float m20, float m21, float m22, float m23,
            float m30, float m31, float m32, float m33);

        // Assignment.
        Matrix& operator= (Matrix const& mat);
        Matrix& operator= (__m128 const* mat);

        // Implicit conversions.
        operator __m128* ();
        operator __m128 const* () const;

        // Access to the slices (rows or columns) of the matrix.
        __m128 const& operator[] (int i) const;
        __m128& operator[] (int i);

    protected:
        // mTable[i] is row i for row-major storage but is column i for
        // column-major order.
        __m128 mTable[4];
    };

public:
    // Logical operations.
    inline static __m128 Not (__m128 const v);                          // ~v
    inline static __m128 And (__m128 const v0, __m128 const v1);        // v0 & v1
    inline static __m128 AndNot (__m128 const v0, __m128 const v1);     // ~v0 & v1
    inline static __m128 Or (__m128 const v0, __m128 const v1);         // v0 | v1
    inline static __m128 Xor (__m128 const v0, __m128 const v1);        // v0 ^ v1
    inline static __m128 Select (__m128 const c, __m128 const v0, __m128 const v1); // (c & v0) | (~c & v1)

    // Comparisons.
    inline static __m128 Equal (__m128 const v0, __m128 const v1);          // v0 == v1
    inline static __m128 NotEqual (__m128 const v0, __m128 const v1);       // v0 != v1
    inline static __m128 Less (__m128 const v0, __m128 const v1);           // v0 < v1
    inline static __m128 LessEqual (__m128 const v0, __m128 const v1);      // v0 <= v1
    inline static __m128 Greater (__m128 const v0, __m128 const v1);        // v0 > v1
    inline static __m128 GreaterEqual (__m128 const v0, __m128 const v1);   // v0 >= v1

    // Vector arithmetic operations.
    inline static __m128 Negate (__m128 const v);                       // -v
    inline static __m128 Add (__m128 const v0, __m128 const v1);        // v0 + v1
    inline static __m128 Subtract (__m128 const v0, __m128 const v1);   // v0 - v1
    inline static __m128 Multiply (__m128 const v0, __m128 const v1);   // v0 * v1
    inline static __m128 Divide (__m128 const v0, __m128 const v1);     // v0 / v1
    inline static __m128 Round (__m128 const v);
    inline static __m128 MaximumAbsoluteComponent (__m128 const v);

    // Vector algebraic operations.
    inline static __m128 Dot (__m128 const v0, __m128 const v1);
    inline static __m128 Length (__m128 const v);
    inline static __m128 LengthRobust (__m128 const v);
    inline static __m128 Normalize (__m128 const v);
    inline static __m128 NormalizeGetLength (__m128 const v, __m128& length);
    inline static __m128 NormalizeRobust (__m128 const v);
    inline static __m128 NormalizeRobustGetLength (__m128 const v, __m128& length);
    inline static __m128 Cross (__m128 const v0, __m128 const v1);

    // Matrix arithmetic operations.
    inline static void Negate (__m128 const* M, __m128* result);
    inline static void Add (__m128 const* A, __m128 const*B, __m128* result);
    inline static void Subtract (__m128 const* A, __m128 const* B, __m128* result);
    inline static void Multiply (__m128 const* M, __m128 const c, __m128* result);
    inline static void Divide (__m128 const* M, __m128 const c, __m128* result);

    // Matrix geometric operations.
    inline static void Transpose (__m128 const* mat, __m128* trn);
    inline static void Inverse (__m128 const* mat, __m128* inv);
    inline static void Adjoint (__m128 const* mat, __m128* adj);
    inline static __m128 Determinant (__m128 const* mat);
    inline static __m128 L1Norm (__m128 const* mat);
    inline static __m128 L2Norm (__m128 const* mat);
    inline static __m128 LInfinityNorm (__m128 const* mat);

    // Matrix-matrix products.
    inline static void MultiplyAB (__m128 const* A, __m128 const* B, __m128* AB);
    inline static void MultiplyATB (__m128 const* A, __m128 const* B, __m128* ATB);
    inline static void MultiplyABT (__m128 const* A, __m128 const* B, __m128* ABT);
    inline static void MultiplyATBT (__m128 const* A, __m128 const* B, __m128* ATBT);
    inline static void MultiplyDM (__m128 const D, __m128 const* M, __m128* DM);
    inline static void MultiplyMD (__m128 const* M, __m128 const D, __m128* MD);

    // Matrix-vector products.
    inline static __m128 MultiplyMV (__m128 const* M, __m128 const V);
    inline static __m128 MultiplyVM (__m128 const V, __m128 const* M);

    // Quaternion support.  In QSlerp, the 't' component must be the splat of
    // a floating-point scalar in [0,1], and q0 and q1 must be unit-length
    // quaternions.
    inline static __m128 QMultiply (__m128 const q0, __m128 const q1);
    inline static __m128 QConjugate (__m128 const q);
    inline static __m128 QInverse (__m128 const q);
    inline static __m128 QSlerp (__m128 const t, __m128 const q0, __m128 const q1);

    // Function evaluations (generally slow, CPU function call per component).
    inline static __m128 Sin (__m128 const v);
    inline static __m128 Cos (__m128 const v);
    inline static __m128 Tan (__m128 const v);
    inline static __m128 ASin (__m128 const v);
    inline static __m128 ACos (__m128 const v);
    inline static __m128 ATan (__m128 const v);

    // Fast function approximations.

    // The SinAppr* functions require |x| <= pi/2.  When x is in this domain,
    // just call SinAppr*(x).  When x is not in the domain, call
    // ReduceAnglesSin(x,y) to obtain y in the domain with sin(x) = sin(y),
    // and then call SinAppr*(y).
    inline static void ReduceAnglesSin (__m128 const x, __m128& y);
    inline static __m128 SinApprDeg11 (__m128 const x);
    inline static __m128 SinApprDeg7 (__m128 const x);

    // The CosAppr* functions require |x| <= pi/2.  When x is in this domain,
    // just call CosAppr*(x, Float4SEE::PPP).  When x is not in the domain,
    // call ReduceAnglesCos(x,y,sign) to obtain y in the domain with cos(x) =
    // sign*cos(y), and then call CosAppr*(y,sign).
    inline static void ReduceAnglesCos (__m128 const x, __m128& y, __m128& sign);
    inline static __m128 CosApprDeg10 (__m128 const x, __m128 const sign);
    inline static __m128 CosApprDeg6 (__m128 const x, __m128 const sign);

    // Integer masks.
    static Vector const ZZZZ;  // (0x00000000, 0x00000000, 0x00000000, 0x00000000)
    static Vector const ZZZF;  // (0x00000000, 0x00000000, 0x00000000, 0xFFFFFFFF)
    static Vector const ZZFZ;  // (0x00000000, 0x00000000, 0xFFFFFFFF, 0x00000000)
    static Vector const ZZFF;  // (0x00000000, 0x00000000, 0xFFFFFFFF, 0xFFFFFFFF)
    static Vector const ZFZZ;  // (0x00000000, 0xFFFFFFFF, 0x00000000, 0x00000000)
    static Vector const ZFZF;  // (0x00000000, 0xFFFFFFFF, 0x00000000, 0xFFFFFFFF)
    static Vector const ZFFZ;  // (0x00000000, 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000)
    static Vector const ZFFF;  // (0x00000000, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF)
    static Vector const FZZZ;  // (0xFFFFFFFF, 0x00000000, 0x00000000, 0x00000000)
    static Vector const FZZF;  // (0xFFFFFFFF, 0x00000000, 0x00000000, 0xFFFFFFFF)
    static Vector const FZFZ;  // (0xFFFFFFFF, 0x00000000, 0xFFFFFFFF, 0x00000000)
    static Vector const FZFF;  // (0xFFFFFFFF, 0x00000000, 0xFFFFFFFF, 0xFFFFFFFF)
    static Vector const FFZZ;  // (0xFFFFFFFF, 0xFFFFFFFF, 0x00000000, 0x00000000)
    static Vector const FFZF;  // (0xFFFFFFFF, 0xFFFFFFFF, 0x00000000, 0xFFFFFFFF)
    static Vector const FFFZ;  // (0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000)
    static Vector const FFFF;  // (0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF)
    static Vector const SIGN;  // (0x80000000, 0x80000000, 0x80000000, 0x80000000)
    static Vector const NSIGN; // (0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF)
    static Vector const NOFRC; // (0x00800000, 0x00800000, 0x00800000, 0x00800000)

    // Numeric constants.
    static Vector const PZZZ;  // (+1.0f,  0.0f,  0.0f,  0.0f)
    static Vector const ZPZZ;  // ( 0.0f, +1.0f,  0.0f,  0.0f)
    static Vector const ZZPZ;  // ( 0.0f,  0.0f, +1.0f,  0.0f)
    static Vector const ZZZP;  // ( 0.0f,  0.0f,  0.0f, +1.0f)
    static Vector const MZZZ;  // (-1.0f,  0.0f,  0.0f,  0.0f)
    static Vector const ZMZZ;  // ( 0.0f, -1.0f,  0.0f,  0.0f)
    static Vector const ZZMZ;  // ( 0.0f,  0.0f, -1.0f,  0.0f)
    static Vector const ZZZM;  // ( 0.0f,  0.0f,  0.0f, -1.0f)
    static Vector const MMMM;  // (-1.0f, -1.0f, -1.0f, -1.0f)
    static Vector const MMMP;  // (-1.0f, -1.0f, -1.0f, +1.0f)
    static Vector const MMPM;  // (-1.0f, -1.0f, +1.0f, -1.0f)
    static Vector const MMPP;  // (-1.0f, -1.0f, +1.0f, +1.0f)
    static Vector const MPMM;  // (-1.0f, +1.0f, -1.0f, -1.0f)
    static Vector const MPMP;  // (-1.0f, +1.0f, -1.0f, +1.0f)
    static Vector const MPPM;  // (-1.0f, +1.0f, +1.0f, -1.0f)
    static Vector const MPPP;  // (-1.0f, +1.0f, +1.0f, +1.0f)
    static Vector const PMMM;  // (+1.0f, -1.0f, -1.0f, -1.0f)
    static Vector const PMMP;  // (+1.0f, -1.0f, -1.0f, +1.0f)
    static Vector const PMPM;  // (+1.0f, -1.0f, +1.0f, -1.0f)
    static Vector const PMPP;  // (+1.0f, -1.0f, +1.0f, +1.0f)
    static Vector const PPMM;  // (+1.0f, +1.0f, -1.0f, -1.0f)
    static Vector const PPMP;  // (+1.0f, +1.0f, -1.0f, +1.0f)
    static Vector const PPPM;  // (+1.0f, +1.0f, +1.0f, -1.0f)
    static Vector const PPPP;  // (+1.0f, +1.0f, +1.0f, +1.0f)
    static Vector const UNIT[4];  // = {PZZZ, ZPZZ, ZZPZ, ZZZP};

    // Constants involving pi.
    static Vector const PI;
    static Vector const HALF_PI;
    static Vector const TWO_PI;
    static Vector const INV_PI;
    static Vector const INV_TWO_PI;

private:
    // Support for computing the adjoint, determinanat, and inverse.
    inline static void GetAdjDet (__m128 const* mat, __m128* adj, __m128* det);

    // Constants to support approximations of sin(x).
    static Vector const C_SIN_APPR_DEG11_0;
    static Vector const C_SIN_APPR_DEG11_1;
    static Vector const C_SIN_APPR_DEG11_2;
    static Vector const C_SIN_APPR_DEG11_3;
    static Vector const C_SIN_APPR_DEG11_4;
    static Vector const C_SIN_APPR_DEG11_5;
    static Vector const C_SIN_APPR_DEG7_0;
    static Vector const C_SIN_APPR_DEG7_1;
    static Vector const C_SIN_APPR_DEG7_2;
    static Vector const C_SIN_APPR_DEG7_3;

    // Constants to support approximations of cos(x).
    static Vector const C_COS_APPR_DEG10_0;
    static Vector const C_COS_APPR_DEG10_1;
    static Vector const C_COS_APPR_DEG10_2;
    static Vector const C_COS_APPR_DEG10_3;
    static Vector const C_COS_APPR_DEG10_4;
    static Vector const C_COS_APPR_DEG10_5;
    static Vector const C_COS_APPR_DEG6_0;
    static Vector const C_COS_APPR_DEG6_1;
    static Vector const C_COS_APPR_DEG6_2;
    static Vector const C_COS_APPR_DEG6_3;
};


// SIMD::Vector

inline SIMD::Vector::Vector()
{
    // Uninitialized.
}

inline SIMD::Vector::Vector(Vector const& vec)
    :
    mTuple(vec.mTuple)
{
}

inline SIMD::Vector::Vector(__m128 const vec)
    :
    mTuple(vec)
{
}

inline SIMD::Vector::Vector(__m128i const vec)
    :
    mTuple(_mm_castsi128_ps(vec))
{
}

inline SIMD::Vector::Vector(float number)
{
    mTuple = _mm_set1_ps(number);
}

inline SIMD::Vector::Vector(float n0, float n1, float n2, float n3)
{
    mTuple = _mm_set_ps(n3, n2, n1, n0);
}

inline SIMD::Vector::Vector(uint32_t encoding)
{
    mTuple = _mm_castsi128_ps(_mm_set1_epi32(encoding));
}

inline SIMD::Vector::Vector(uint32_t e0, uint32_t e1, uint32_t e2, uint32_t e3)
{
    mTuple = _mm_castsi128_ps(_mm_set_epi32(e3, e2, e1, e0));
}

inline SIMD::Vector& SIMD::Vector::operator= (Vector const& vec)
{
    mTuple = vec.mTuple;
    return *this;
}

inline SIMD::Vector& SIMD::Vector::operator= (__m128 const vec)
{
    mTuple = vec;
    return *this;
}

inline SIMD::Vector& SIMD::Vector::operator= (__m128i const vec)
{
    mTuple = _mm_castsi128_ps(vec);
    return *this;
}

inline SIMD::Vector::operator __m128 ()
{
    return mTuple;
}

inline SIMD::Vector::operator __m128 () const
{
    return mTuple;
}

inline SIMD::Vector::operator __m128i ()
{
    return _mm_castps_si128(mTuple);
}

inline SIMD::Vector::operator __m128i () const
{
    return _mm_castps_si128(mTuple);
}



// SIMD::Matrix

inline SIMD::Matrix::Matrix()
{
    // Uninitialized.
}

inline SIMD::Matrix::Matrix(Matrix const& mat)
{
    mTable[0] = mat.mTable[0];
    mTable[1] = mat.mTable[1];
    mTable[2] = mat.mTable[2];
    mTable[3] = mat.mTable[3];
}

inline SIMD::Matrix::Matrix(__m128 const* mat)
{
    mTable[0] = mat[0];
    mTable[1] = mat[1];
    mTable[2] = mat[2];
    mTable[3] = mat[3];
}

inline SIMD::Matrix::Matrix(
    float m00, float m01, float m02, float m03,
    float m10, float m11, float m12, float m13,
    float m20, float m21, float m22, float m23,
    float m30, float m31, float m32, float m33)
{
#if defined(GTE_USE_ROW_MAJOR)
    mTable[0] = _mm_setr_ps(m00, m01, m02, m03);
    mTable[1] = _mm_setr_ps(m10, m11, m12, m13);
    mTable[2] = _mm_setr_ps(m20, m21, m22, m23);
    mTable[3] = _mm_setr_ps(m30, m31, m32, m33);
#else
    mTable[0] = _mm_setr_ps(m00, m10, m20, m30);
    mTable[1] = _mm_setr_ps(m01, m11, m21, m31);
    mTable[2] = _mm_setr_ps(m02, m12, m22, m32);
    mTable[3] = _mm_setr_ps(m03, m13, m23, m33);
#endif
}

inline SIMD::Matrix& SIMD::Matrix::operator= (Matrix const& mat)
{
    mTable[0] = mat.mTable[0];
    mTable[1] = mat.mTable[1];
    mTable[2] = mat.mTable[2];
    mTable[3] = mat.mTable[3];
    return *this;
}

inline SIMD::Matrix& SIMD::Matrix::operator= (__m128 const* mat)
{
    mTable[0] = mat[0];
    mTable[1] = mat[1];
    mTable[2] = mat[2];
    mTable[3] = mat[3];
    return *this;
}

inline SIMD::Matrix::operator __m128* ()
{
    return mTable;
}

inline SIMD::Matrix::operator __m128 const* () const
{
    return mTable;
}

inline __m128 const& SIMD::Matrix::operator[] (int i) const
{
    return mTable[i];
}

inline __m128& SIMD::Matrix::operator[] (int i)
{
    return mTable[i];
}



// Logical operations.

inline __m128 SIMD::Not(__m128 const v)
{
    return _mm_xor_ps(v, FFFF);
}

inline __m128 SIMD::And(__m128 const v0, __m128 const v1)
{
    return _mm_and_ps(v0, v1);
}

inline __m128 SIMD::AndNot(__m128 const v0, __m128 const v1)
{
    return _mm_andnot_ps(v0, v1);
}

inline __m128 SIMD::Or(__m128 const v0, __m128 const v1)
{
    return _mm_or_ps(v0, v1);
}

inline __m128 SIMD::Xor(__m128 const v0, __m128 const v1)
{
    return _mm_xor_ps(v0, v1);
}

inline __m128 SIMD::Select(__m128 const c, __m128 const v0, __m128 const v1)
{
    return _mm_or_ps(_mm_and_ps(c, v0), _mm_andnot_ps(c, v1));
}



// Comparisons.

inline __m128 SIMD::Equal(__m128 const v0, __m128 const v1)
{
    return _mm_cmpeq_ps(v0, v1);
}

inline __m128 SIMD::NotEqual(__m128 const v0, __m128 const v1)
{
    return _mm_cmpneq_ps(v0, v1);
}

inline __m128 SIMD::Less(__m128 const v0, __m128 const v1)
{
    return _mm_cmplt_ps(v0, v1);
}

inline __m128 SIMD::LessEqual(__m128 const v0, __m128 const v1)
{
    return _mm_cmple_ps(v0, v1);
}

inline __m128 SIMD::Greater(__m128 const v0, __m128 const v1)
{
    return _mm_cmpgt_ps(v0, v1);
}

inline __m128 SIMD::GreaterEqual(__m128 const v0, __m128 const v1)
{
    return _mm_cmpge_ps(v0, v1);
}



// Vector arithmetic operations.

inline __m128 SIMD::Negate(__m128 const v)
{
    return _mm_xor_ps(v, SIGN);
}

inline __m128 SIMD::Add(__m128 const v0, __m128 const v1)
{
    return _mm_add_ps(v0, v1);
}

inline __m128 SIMD::Subtract(__m128 const v0, __m128 const v1)
{
    return _mm_sub_ps(v0, v1);
}

inline __m128 SIMD::Multiply(__m128 const v0, __m128 const v1)
{
    return _mm_mul_ps(v0, v1);
}

inline __m128 SIMD::Divide(__m128 const v0, __m128 const v1)
{
    return _mm_div_ps(v0, v1);
}

inline __m128 SIMD::Round(__m128 const v)
{
    __m128 t0 = _mm_and_ps(NSIGN, v);
    t0 = _mm_castsi128_ps(_mm_cmplt_epi32(_mm_castps_si128(t0), NOFRC));
    __m128i t1 = _mm_cvtps_epi32(v);  // float-to-int
    __m128 t2 = _mm_cvtepi32_ps(t1);  // int-to-float
    t2 = _mm_and_ps(t2, t0);
    t0 = _mm_andnot_ps(t0, v);
    t2 = _mm_or_ps(t2, t0);
    return t2;
}

inline __m128 SIMD::MaximumAbsoluteComponent(__m128 const v)
{
    __m128 vAbs = _mm_andnot_ps(SIGN, v);
    __m128 max0 = _mm_shuffle_ps(vAbs, vAbs, _MM_SHUFFLE(0, 0, 0, 0));
    __m128 max1 = _mm_shuffle_ps(vAbs, vAbs, _MM_SHUFFLE(1, 1, 1, 1));
    __m128 max2 = _mm_shuffle_ps(vAbs, vAbs, _MM_SHUFFLE(2, 2, 2, 2));
    __m128 max3 = _mm_shuffle_ps(vAbs, vAbs, _MM_SHUFFLE(3, 3, 3, 3));
    max0 = _mm_max_ps(max0, max1);
    max2 = _mm_max_ps(max2, max3);
    max0 = _mm_max_ps(max0, max2);
    return max0;
}



// Vector algebraic operations.

inline __m128 SIMD::Dot(__m128 const v0, __m128 const v1)
{
    // (x0*x1, y0*y1, z0*z1, w0*w1)
    __m128 t0 = _mm_mul_ps(v0, v1);

    // (y0*y1, x0*x1, w0*w1, z0*z1)
    __m128 t1 = _mm_shuffle_ps(t0, t0, _MM_SHUFFLE(2, 3, 0, 1));

    // (x0*x1 + y0*y1, x0*x1 + y0*y1, z0*z1 + w0*w1, z0*z1 + w0*w1)
    __m128 t2 = _mm_add_ps(t0, t1);

    // (z0*z1 + w0*w1, z0*z1 + w0*w1, x0*x1 + y0*y1, x0*x1 + y0*y1)
    __m128 t3 = _mm_shuffle_ps(t2, t2, _MM_SHUFFLE(0, 0, 2, 2));

    // (dot, dot, dot, dot)
    __m128 dotSplat = _mm_add_ps(t2, t3);
    return dotSplat;
}

inline __m128 SIMD::Length(__m128 const v)
{
    __m128 sqrLength = Dot(v, v);
    return _mm_sqrt_ps(sqrLength);
}

inline __m128 SIMD::LengthRobust(__m128 const v)
{
    // Compute the maximum absolute value component.
    __m128 maxComponent = MaximumAbsoluteComponent(v);

    // Divide by the maximum absolute component.  This is potentially a
    // divide by zero.
    __m128 normalized = _mm_div_ps(v, maxComponent);

    // Set to zero when the original length is zero.
    __m128 mask = _mm_cmpneq_ps(ZZZZ, maxComponent);
    normalized = _mm_and_ps(mask, normalized);

    // (sqrLength, sqrLength, sqrLength, sqrLength)
    __m128 sqrLength = Dot(normalized, normalized);
    __m128 length = _mm_sqrt_ps(sqrLength);
    return length;
}

inline __m128 SIMD::Normalize(__m128 const v)
{
    // (sqrLength, sqrLength, sqrLength, sqrLength)
    __m128 sqrLength = Dot(v, v);

    // (length, length, length, length)
    __m128 length = _mm_sqrt_ps(sqrLength);

    // Divide by the length to normalize.  This is potentially a divide by
    // zero or a divide by infinity.
    __m128 normalized = _mm_div_ps(v, length);

    // Set to zero when the original length is zero.
    __m128 mask = _mm_cmpneq_ps(ZZZZ, length);
    normalized = _mm_and_ps(mask, normalized);
    return normalized;
}

inline __m128 SIMD::NormalizeGetLength(__m128 const v, __m128& length)
{
    // (sqrLength, sqrLength, sqrLength, sqrLength)
    __m128 sqrLength = Dot(v, v);

    // (length, length, length, length)
    length = _mm_sqrt_ps(sqrLength);

    // Divide by the length to normalize.  This is potentially a divide by
    // zero or a divide by infinity.
    __m128 normalized = _mm_div_ps(v, length);

    // Set to zero when the original length is zero.
    __m128 mask = _mm_cmpneq_ps(ZZZZ, length);
    normalized = _mm_and_ps(mask, normalized);
    length = _mm_and_ps(mask, length);
    return normalized;
}

inline __m128 SIMD::NormalizeRobust(__m128 const v)
{
    // Compute the maximum absolute value component.
    __m128 maxComponent = MaximumAbsoluteComponent(v);

    // Divide by the maximum absolute component.  This is potentially a
    // divide by zero.
    __m128 normalized = _mm_div_ps(v, maxComponent);

    // Set to zero when the original length is zero.
    __m128 mask = _mm_cmpneq_ps(ZZZZ, maxComponent);
    normalized = _mm_and_ps(mask, normalized);

    // (sqrLength, sqrLength, sqrLength, sqrLength)
    __m128 sqrLength = Dot(normalized, normalized);

    // (length, length, length, length)
    __m128 length = _mm_sqrt_ps(sqrLength);

    // Divide by the length to normalize.  This is potentially a divide by
    // zero.
    normalized = _mm_div_ps(normalized, length);

    // Set to zero when the original length is zero or infinity.  In the
    // latter case, this is considered to be an unexpected condition.
    normalized = _mm_and_ps(mask, normalized);
    return normalized;
}

inline __m128 SIMD::NormalizeRobustGetLength(__m128 const v, __m128& length)
{
    // Compute the maximum absolute value component.
    __m128 maxComponent = MaximumAbsoluteComponent(v);

    // Divide by the maximum absolute component.  This is potentially a
    // divide by zero.
    __m128 normalized = _mm_div_ps(v, maxComponent);

    // Set to zero when the original length is zero.
    __m128 mask = _mm_cmpneq_ps(ZZZZ, maxComponent);
    normalized = _mm_and_ps(mask, normalized);

    // (sqrLength, sqrLength, sqrLength, sqrLength)
    __m128 sqrLength = Dot(normalized, normalized);

    // (length, length, length, length)
    length = _mm_sqrt_ps(sqrLength);

    // Divide by the length to normalize.  This is potentially a divide by
    // zero.
    normalized = _mm_div_ps(normalized, length);
    length = _mm_mul_ps(length, maxComponent);  // true length

    // Set to zero when the original length is zero or infinity.  In the
    // latter case, this is considered to be an unexpected condition.
    normalized = _mm_and_ps(mask, normalized);
    length = _mm_and_ps(mask, length);
    return normalized;
}

inline __m128 SIMD::Cross(__m128 const v0, __m128 const v1)
{
    // v0 = (x0, y0, z0, 0), v1 = (x1, y1, z1, 0)
    // cross = (y0*z1 - z0*y1, z0*x1 - x0*z1, x0*y1 - y0*x1, 0)

    // (y0, z0, x0, 0)
    __m128 t0 = _mm_shuffle_ps(v0, v0, _MM_SHUFFLE(3, 0, 2, 1));

    // (z1, x1, y1, 0)
    __m128 t1 = _mm_shuffle_ps(v1, v1, _MM_SHUFFLE(3, 1, 0, 2));

    // (y0*z1, z0*x1, x0*y1, 0)
    __m128 cross = _mm_mul_ps(t0, t1);

    // (z0, x0, y0, 0), computed from t0, not v0
    t0 = _mm_shuffle_ps(t0, t0, _MM_SHUFFLE(3, 0, 2, 1));

    // (y1, z1, x1, 0), computed from t1, not v1
    t1 = _mm_shuffle_ps(t1, t1, _MM_SHUFFLE(3, 1, 0, 2));

    // (z0*y1, x0*z1, y0*x1, 0)
    t0 = _mm_mul_ps(t0, t1);

    // (y0*z1 - z0*y1, z0*x1 - x0*z1, x0*y1 - y0*x1, 0)
    cross = _mm_sub_ps(cross, t0);
    return cross;
}



// Matrix arithmetic operations.

inline void SIMD::Negate(__m128 const* M, __m128* result)
{
    result[0] = _mm_sub_ps(SIMD::ZZZZ, M[0]);
    result[1] = _mm_sub_ps(SIMD::ZZZZ, M[1]);
    result[2] = _mm_sub_ps(SIMD::ZZZZ, M[2]);
    result[3] = _mm_sub_ps(SIMD::ZZZZ, M[3]);
}

inline void SIMD::Add(__m128 const* A, __m128 const* B, __m128* result)
{
    result[0] = _mm_add_ps(A[0], B[0]);
    result[1] = _mm_add_ps(A[1], B[1]);
    result[2] = _mm_add_ps(A[2], B[2]);
    result[3] = _mm_add_ps(A[3], B[3]);
}

inline void SIMD::Subtract(__m128 const* A, __m128 const*B, __m128* result)
{
    result[0] = _mm_sub_ps(A[0], B[0]);
    result[1] = _mm_sub_ps(A[1], B[1]);
    result[2] = _mm_sub_ps(A[2], B[2]);
    result[3] = _mm_sub_ps(A[3], B[3]);
}

inline void SIMD::Multiply(__m128 const* M, __m128 const c, __m128* result)
{
    result[0] = _mm_mul_ps(M[0], c);
    result[1] = _mm_mul_ps(M[1], c);
    result[2] = _mm_mul_ps(M[2], c);
    result[3] = _mm_mul_ps(M[3], c);
}

inline void SIMD::Divide(__m128 const* M, __m128 const c, __m128* result)
{
    result[0] = _mm_div_ps(M[0], c);
    result[1] = _mm_div_ps(M[1], c);
    result[2] = _mm_div_ps(M[2], c);
    result[3] = _mm_div_ps(M[3], c);
}



// Matrix geometric operations.

inline void SIMD::Transpose(__m128 const* mat, __m128* trn)
{
    // VM:(m00, m01, m10, m11), MV:(m00, m10, m01, m11)
    __m128 s0 = _mm_shuffle_ps(mat[0], mat[1], _MM_SHUFFLE(1, 0, 1, 0));
    // VM:(m20, m21, m30, m31), MV:(m02, m12, m03, m13)
    __m128 s1 = _mm_shuffle_ps(mat[2], mat[3], _MM_SHUFFLE(1, 0, 1, 0));
    // VM:(m02, m03, m12, m13), MV:(m20, m30, m21, m31)
    __m128 s2 = _mm_shuffle_ps(mat[0], mat[1], _MM_SHUFFLE(3, 2, 3, 2));
    // VM:(m22, m23, m32, m33), MV:(m22, m32, m23, m33)
    __m128 s3 = _mm_shuffle_ps(mat[2], mat[3], _MM_SHUFFLE(3, 2, 3, 2));

    // VM:(m00, m10, m20, m30), MV:(m00, m01, m02, m03)
    trn[0] = _mm_shuffle_ps(s0, s1, _MM_SHUFFLE(2, 0, 2, 0));
    // VM:(m01, m11, m21, m31), MV:(m10, m11, m12, m13)
    trn[1] = _mm_shuffle_ps(s0, s1, _MM_SHUFFLE(3, 1, 3, 1));
    // VM:(m02, m12, m22, m32), MV:(m20, m21, m22, m23)
    trn[2] = _mm_shuffle_ps(s2, s3, _MM_SHUFFLE(2, 0, 2, 0));
    // VM:(m03, m13, m23, m33), MV:(m30, m31, m32, m33)
    trn[3] = _mm_shuffle_ps(s2, s3, _MM_SHUFFLE(3, 1, 3, 1));
}

inline void SIMD::Inverse(__m128 const* mat, __m128* inv)
{
    __m128 det;
    GetAdjDet(mat, inv, &det);

    // Compute the reciprocal of the determinant.  Guard against a division
    // by zero.  When the determinant is zero, the zero matrix is returned.
    __m128 invDet = _mm_div_ps(PPPP, det);
    __m128 neqZero = _mm_cmpneq_ps(det, ZZZZ);
    invDet = _mm_and_ps(neqZero, invDet);

    inv[0] = _mm_mul_ps(inv[0], invDet);
    inv[1] = _mm_mul_ps(inv[1], invDet);
    inv[2] = _mm_mul_ps(inv[2], invDet);
    inv[3] = _mm_mul_ps(inv[3], invDet);
}

inline void SIMD::Adjoint(__m128 const* mat, __m128* adj)
{
    GetAdjDet(mat, adj, 0);
}

inline __m128 SIMD::Determinant(__m128 const* mat)
{
    __m128 det;
    GetAdjDet(mat, 0, &det);
    return det;
}

inline __m128 SIMD::L1Norm(__m128 const* mat)
{
    __m128 sum = _mm_andnot_ps(SIMD::SIGN, mat[0]);
    __m128 tmp = _mm_andnot_ps(SIMD::SIGN, mat[1]);
    sum = _mm_add_ps(sum, tmp);
    tmp = _mm_andnot_ps(SIMD::SIGN, mat[2]);
    sum = _mm_add_ps(sum, tmp);
    tmp = _mm_andnot_ps(SIMD::SIGN, mat[3]);
    sum = _mm_add_ps(sum, tmp);
    __m128 norm = SIMD::Dot(sum, SIMD::PPPP);
    return norm;
}

inline __m128 SIMD::L2Norm(__m128 const* mat)
{
    __m128 sum = _mm_mul_ps(mat[0], mat[0]);
    __m128 tmp = _mm_mul_ps(mat[1], mat[1]);
    sum = _mm_add_ps(sum, tmp);
    tmp = _mm_mul_ps(mat[2], mat[2]);
    sum = _mm_add_ps(sum, tmp);
    tmp = _mm_mul_ps(mat[3], mat[3]);
    sum = _mm_add_ps(sum, tmp);
    __m128 norm = SIMD::Dot(sum, SIMD::PPPP);
    return norm;
}

inline __m128 SIMD::LInfinityNorm(__m128 const* mat)
{
    __m128 max = SIMD::MaximumAbsoluteComponent(mat[0]);
    __m128 tmp = SIMD::MaximumAbsoluteComponent(mat[1]);
    max = _mm_max_ps(max, tmp);
    tmp = SIMD::MaximumAbsoluteComponent(mat[2]);
    max = _mm_max_ps(max, tmp);
    tmp = SIMD::MaximumAbsoluteComponent(mat[3]);
    max = _mm_max_ps(max, tmp);
    return max;
}



// Matrix-matrix products.

inline void SIMD::MultiplyAB(__m128 const* A, __m128 const* B, __m128* AB)
{
    __m128 t0, t1, t2, t3;

#if defined(GTE_USE_ROW_MAJOR)
    t0 = _mm_shuffle_ps(A[0], A[0], _MM_SHUFFLE(0, 0, 0, 0));
    t1 = _mm_shuffle_ps(A[0], A[0], _MM_SHUFFLE(1, 1, 1, 1));
    t2 = _mm_shuffle_ps(A[0], A[0], _MM_SHUFFLE(2, 2, 2, 2));
    t3 = _mm_shuffle_ps(A[0], A[0], _MM_SHUFFLE(3, 3, 3, 3));
    t0 = _mm_mul_ps(t0, B[0]);
    t1 = _mm_mul_ps(t1, B[1]);
    t2 = _mm_mul_ps(t2, B[2]);
    t3 = _mm_mul_ps(t3, B[3]);
    t0 = _mm_add_ps(t0, t1);
    t2 = _mm_add_ps(t2, t3);
    AB[0] = _mm_add_ps(t0, t2);

    t0 = _mm_shuffle_ps(A[1], A[1], _MM_SHUFFLE(0, 0, 0, 0));
    t1 = _mm_shuffle_ps(A[1], A[1], _MM_SHUFFLE(1, 1, 1, 1));
    t2 = _mm_shuffle_ps(A[1], A[1], _MM_SHUFFLE(2, 2, 2, 2));
    t3 = _mm_shuffle_ps(A[1], A[1], _MM_SHUFFLE(3, 3, 3, 3));
    t0 = _mm_mul_ps(t0, B[0]);
    t1 = _mm_mul_ps(t1, B[1]);
    t2 = _mm_mul_ps(t2, B[2]);
    t3 = _mm_mul_ps(t3, B[3]);
    t0 = _mm_add_ps(t0, t1);
    t2 = _mm_add_ps(t2, t3);
    AB[1] = _mm_add_ps(t0, t2);

    t0 = _mm_shuffle_ps(A[2], A[2], _MM_SHUFFLE(0, 0, 0, 0));
    t1 = _mm_shuffle_ps(A[2], A[2], _MM_SHUFFLE(1, 1, 1, 1));
    t2 = _mm_shuffle_ps(A[2], A[2], _MM_SHUFFLE(2, 2, 2, 2));
    t3 = _mm_shuffle_ps(A[2], A[2], _MM_SHUFFLE(3, 3, 3, 3));
    t0 = _mm_mul_ps(t0, B[0]);
    t1 = _mm_mul_ps(t1, B[1]);
    t2 = _mm_mul_ps(t2, B[2]);
    t3 = _mm_mul_ps(t3, B[3]);
    t0 = _mm_add_ps(t0, t1);
    t2 = _mm_add_ps(t2, t3);
    AB[2] = _mm_add_ps(t0, t2);

    t0 = _mm_shuffle_ps(A[3], A[3], _MM_SHUFFLE(0, 0, 0, 0));
    t1 = _mm_shuffle_ps(A[3], A[3], _MM_SHUFFLE(1, 1, 1, 1));
    t2 = _mm_shuffle_ps(A[3], A[3], _MM_SHUFFLE(2, 2, 2, 2));
    t3 = _mm_shuffle_ps(A[3], A[3], _MM_SHUFFLE(3, 3, 3, 3));
    t0 = _mm_mul_ps(t0, B[0]);
    t1 = _mm_mul_ps(t1, B[1]);
    t2 = _mm_mul_ps(t2, B[2]);
    t3 = _mm_mul_ps(t3, B[3]);
    t0 = _mm_add_ps(t0, t1);
    t2 = _mm_add_ps(t2, t3);
    AB[3] = _mm_add_ps(t0, t2);
#else
    t0 = _mm_shuffle_ps(B[0], B[0], _MM_SHUFFLE(0, 0, 0, 0));
    t1 = _mm_shuffle_ps(B[0], B[0], _MM_SHUFFLE(1, 1, 1, 1));
    t2 = _mm_shuffle_ps(B[0], B[0], _MM_SHUFFLE(2, 2, 2, 2));
    t3 = _mm_shuffle_ps(B[0], B[0], _MM_SHUFFLE(3, 3, 3, 3));
    t0 = _mm_mul_ps(t0, A[0]);
    t1 = _mm_mul_ps(t1, A[1]);
    t2 = _mm_mul_ps(t2, A[2]);
    t3 = _mm_mul_ps(t3, A[3]);
    t0 = _mm_add_ps(t0, t1);
    t2 = _mm_add_ps(t2, t3);
    AB[0] = _mm_add_ps(t0, t2);

    t0 = _mm_shuffle_ps(B[1], B[1], _MM_SHUFFLE(0, 0, 0, 0));
    t1 = _mm_shuffle_ps(B[1], B[1], _MM_SHUFFLE(1, 1, 1, 1));
    t2 = _mm_shuffle_ps(B[1], B[1], _MM_SHUFFLE(2, 2, 2, 2));
    t3 = _mm_shuffle_ps(B[1], B[1], _MM_SHUFFLE(3, 3, 3, 3));
    t0 = _mm_mul_ps(t0, A[0]);
    t1 = _mm_mul_ps(t1, A[1]);
    t2 = _mm_mul_ps(t2, A[2]);
    t3 = _mm_mul_ps(t3, A[3]);
    t0 = _mm_add_ps(t0, t1);
    t2 = _mm_add_ps(t2, t3);
    AB[1] = _mm_add_ps(t0, t2);

    t0 = _mm_shuffle_ps(B[2], B[2], _MM_SHUFFLE(0, 0, 0, 0));
    t1 = _mm_shuffle_ps(B[2], B[2], _MM_SHUFFLE(1, 1, 1, 1));
    t2 = _mm_shuffle_ps(B[2], B[2], _MM_SHUFFLE(2, 2, 2, 2));
    t3 = _mm_shuffle_ps(B[2], B[2], _MM_SHUFFLE(3, 3, 3, 3));
    t0 = _mm_mul_ps(t0, A[0]);
    t1 = _mm_mul_ps(t1, A[1]);
    t2 = _mm_mul_ps(t2, A[2]);
    t3 = _mm_mul_ps(t3, A[3]);
    t0 = _mm_add_ps(t0, t1);
    t2 = _mm_add_ps(t2, t3);
    AB[2] = _mm_add_ps(t0, t2);

    t0 = _mm_shuffle_ps(B[3], B[3], _MM_SHUFFLE(0, 0, 0, 0));
    t1 = _mm_shuffle_ps(B[3], B[3], _MM_SHUFFLE(1, 1, 1, 1));
    t2 = _mm_shuffle_ps(B[3], B[3], _MM_SHUFFLE(2, 2, 2, 2));
    t3 = _mm_shuffle_ps(B[3], B[3], _MM_SHUFFLE(3, 3, 3, 3));
    t0 = _mm_mul_ps(t0, A[0]);
    t1 = _mm_mul_ps(t1, A[1]);
    t2 = _mm_mul_ps(t2, A[2]);
    t3 = _mm_mul_ps(t3, A[3]);
    t0 = _mm_add_ps(t0, t1);
    t2 = _mm_add_ps(t2, t3);
    AB[3] = _mm_add_ps(t0, t2);
#endif
}

inline void SIMD::MultiplyATB(__m128 const* A, __m128 const* B, __m128* ATB)
{
    __m128 ATrn[4];
    Transpose(A, ATrn);
    MultiplyAB(ATrn, B, ATB);
}

inline void SIMD::MultiplyABT(__m128 const* A, __m128 const* B, __m128* ABT)
{
    __m128 BTrn[4];
    Transpose(B, BTrn);
    MultiplyAB(A, BTrn, ABT);
}

inline void SIMD::MultiplyATBT(__m128 const* A, __m128 const* B, __m128* ATBT)
{
    __m128 BA[4];
    MultiplyAB(B, A, BA);
    Transpose(BA, ATBT);
}

inline void SIMD::MultiplyDM(__m128 const D, __m128 const* M, __m128* DM)
{
#if defined(GTE_USE_ROW_MAJOR)
    DM[0] = _mm_mul_ps(D, M[0]);
    DM[1] = _mm_mul_ps(D, M[1]);
    DM[2] = _mm_mul_ps(D, M[2]);
    DM[3] = _mm_mul_ps(D, M[3]);
#else
    __m128 d0 = _mm_shuffle_ps(D, D, _MM_SHUFFLE(0, 0, 0, 0));
    __m128 d1 = _mm_shuffle_ps(D, D, _MM_SHUFFLE(1, 1, 1, 1));
    __m128 d2 = _mm_shuffle_ps(D, D, _MM_SHUFFLE(2, 2, 2, 2));
    __m128 d3 = _mm_shuffle_ps(D, D, _MM_SHUFFLE(3, 3, 3, 3));
    DM[0] = _mm_mul_ps(d0, M[0]);
    DM[1] = _mm_mul_ps(d1, M[1]);
    DM[2] = _mm_mul_ps(d2, M[2]);
    DM[3] = _mm_mul_ps(d3, M[3]);
#endif
}

inline void SIMD::MultiplyMD(__m128 const* M, __m128 const D, __m128* MD)
{
#if defined(GTE_USE_ROW_MAJOR)
    __m128 d0 = _mm_shuffle_ps(D, D, _MM_SHUFFLE(0, 0, 0, 0));
    __m128 d1 = _mm_shuffle_ps(D, D, _MM_SHUFFLE(1, 1, 1, 1));
    __m128 d2 = _mm_shuffle_ps(D, D, _MM_SHUFFLE(2, 2, 2, 2));
    __m128 d3 = _mm_shuffle_ps(D, D, _MM_SHUFFLE(3, 3, 3, 3));
    MD[0] = _mm_mul_ps(M[0], d0);
    MD[1] = _mm_mul_ps(M[1], d1);
    MD[2] = _mm_mul_ps(M[2], d2);
    MD[3] = _mm_mul_ps(M[3], d3);
#else
    MD[0] = _mm_mul_ps(M[0], D);
    MD[1] = _mm_mul_ps(M[1], D);
    MD[2] = _mm_mul_ps(M[2], D);
    MD[3] = _mm_mul_ps(M[3], D);
#endif
}



// Matrix-vector products.

inline __m128 SIMD::MultiplyMV(__m128 const* M, __m128 const V)
{
#if defined(GTE_USE_ROW_MAJOR)
    __m128 MTrn[4];
    Transpose(M, MTrn);
    return MultiplyVM(V, MTrn);
#else
    // u0 = m00*v0 + m01*v1 + m02*v2 + m03*v3
    // u1 = m10*v0 + m11*v1 + m12*v2 + m13*v3
    // u2 = m20*v0 + m21*v1 + m22*v2 + m23*v3
    // u3 = m30*v0 + m31*v1 + m32*v2 + m33*v3

    // (v0, v0, v0, v0)
    __m128 t0 = _mm_shuffle_ps(V, V, _MM_SHUFFLE(0, 0, 0, 0));
    // (v1, v1, v1, v1)
    __m128 t1 = _mm_shuffle_ps(V, V, _MM_SHUFFLE(1, 1, 1, 1));
    // (v2, v2, v2, v2)
    __m128 t2 = _mm_shuffle_ps(V, V, _MM_SHUFFLE(2, 2, 2, 2));
    // (v3, v3, v3, v3)
    __m128 t3 = _mm_shuffle_ps(V, V, _MM_SHUFFLE(3, 3, 3, 3));

    // (m00*v0, m10*v0, m20*v0, m30*v0)
    t0 = _mm_mul_ps(M[0], t0);
    // (m01*v1, m11*v1, m21*v1, m31*v1)
    t1 = _mm_mul_ps(M[1], t1);
    // (m02*v2, m12*v2, m22*v2, m32*v2)
    t2 = _mm_mul_ps(M[2], t2);
    // (m03*v3, m13*v3, m23*v3, m33*v3)
    t3 = _mm_mul_ps(M[3], t3);

    // (m00*v0+m01*v1, m10*v0+m11*v1, m20*v0+m21*v1, m30*v0+m31*v1)
    t0 = _mm_add_ps(t0, t1);
    // (m02*v2+m03*v3, m12*v2+m13*v3, m22*v2+m23*v3, m32*v2+m33*v3)
    t2 = _mm_add_ps(t2, t3);

    // M*V
    t0 = _mm_add_ps(t0, t2);
    return t0;
#endif
}

inline __m128 SIMD::MultiplyVM(__m128 const V, __m128 const* M)
{
#if defined(GTE_USE_ROW_MAJOR)
    // u0 = v0*m00 + v1*m10 + v2*m20 + v3*m30
    // u1 = v0*m01 + v1*m11 + v2*m21 + v3*m31
    // u2 = v0*m02 + v1*m12 + v2*m22 + v3*m32
    // u3 = v0*m03 + v1*m13 + v2*m23 + v3*m33

    // (v0, v0, v0, v0)
    __m128 t0 = _mm_shuffle_ps(V, V, _MM_SHUFFLE(0, 0, 0, 0));
    // (v1, v1, v1, v1)
    __m128 t1 = _mm_shuffle_ps(V, V, _MM_SHUFFLE(1, 1, 1, 1));
    // (v2, v2, v2, v2)
    __m128 t2 = _mm_shuffle_ps(V, V, _MM_SHUFFLE(2, 2, 2, 2));
    // (v3, v3, v3, v3)
    __m128 t3 = _mm_shuffle_ps(V, V, _MM_SHUFFLE(3, 3, 3, 3));

    // (v0*m00, v0*m01, v0*m02, v0*m03)
    t0 = _mm_mul_ps(t0, M[0]);
    // (v1*m10, v1*m11, v1*m12, v1*m13)
    t1 = _mm_mul_ps(t1, M[1]);
    // (v2*m20, v2*m21, v2*m22, v2*m23)
    t2 = _mm_mul_ps(t2, M[2]);
    // (v3*m30, v3*m31, v3*m32, v3*m33)
    t3 = _mm_mul_ps(t3, M[3]);

    // (v0*m00+v1*m10, v0*m01+v1*m11, v0*m02+v1*m12, v0*m03+v1*m13)
    t0 = _mm_add_ps(t0, t1);
    // (v2*m20+v3*m30, v2*m21+v3*m31, v2*m22+v3*m32, v2*m23+v3*m33)
    t2 = _mm_add_ps(t2, t3);

    // V*M
    t0 = _mm_add_ps(t0, t2);
    return t0;
#else
    __m128 MTrn[4];
    Transpose(M, MTrn);
    return MultiplyMV(MTrn, V);
#endif
}



// Quaternion support.

inline __m128 SIMD::QMultiply(__m128 const q0, __m128 const q1)
{
    // (x0*i + y0*j + z0*k + w0)*(x1*i + y1*j + z1*k + w1)
    // =
    // i*(+x0*w1 + y0*z1 - z0*y1 + w0*x1) +
    // j*(-x0*z1 + y0*w1 + z0*x1 + w0*y1) +
    // k*(+x0*y1 - y0*x1 + z0*w1 + w0*z1) +
    // 1*(-x0*x1 - y0*y1 - z0*z1 + w0*w1)

    __m128 product;
    {
        // (x0, x0, x0, x0)
        __m128 t0 = _mm_shuffle_ps(q0, q0, _MM_SHUFFLE(0, 0, 0, 0));
        // (w1, z1, y1, x1)
        __m128 t1 = _mm_shuffle_ps(q1, q1, _MM_SHUFFLE(0, 1, 2, 3));
        // (+w1, -z1, +y1, -x1)
        t1 = _mm_mul_ps(SIMD::PMPM, t1);
        // (+x0*w1, -x0*z1, +x0*y1, -x0*x1)
        product = _mm_mul_ps(t0, t1);
    }
    {
        // (y0, y0, y0, y0)
        __m128 t0 = _mm_shuffle_ps(q0, q0, _MM_SHUFFLE(1, 1, 1, 1));
        // (z1, w1, x1, y1)
        __m128 t1 = _mm_shuffle_ps(q1, q1, _MM_SHUFFLE(1, 0, 3, 2));
        // (+z1, +w1, -x1, -y1)
        t1 = _mm_mul_ps(SIMD::PPMM, t1);
        // product += (+y0*z1, +y0*w1, -y0*x1, -y0*y1)
        t1 = _mm_mul_ps(t0, t1);
        product = _mm_add_ps(product, t1);
    }
    {
        // (z0, z0, z0, z0)
        __m128 t0 = _mm_shuffle_ps(q0, q0, _MM_SHUFFLE(2, 2, 2, 2));
        // (y1, x1, w1, z1)
        __m128 t1 = _mm_shuffle_ps(q1, q1, _MM_SHUFFLE(2, 3, 0, 1));
        // (-y1, +x1, +w1, -z1)
        t1 = _mm_mul_ps(SIMD::MPPM, t1);
        // product += (-z0*y1, +z0*x1, +z0*w1, -z0*z1)
        t1 = _mm_mul_ps(t0, t1);
        product = _mm_add_ps(product, t1);
    }
    {
        // (w0, w0, w0, w0)
        __m128 t0 = _mm_shuffle_ps(q0, q0, _MM_SHUFFLE(3, 3, 3, 3));
        // (+w0*x1, +w0*y1, +w0*z1, +w0*w1)
        __m128 t1 = _mm_mul_ps(t0, q1);
        // product += (+w0*x1, +w0*y1, +w0*z1, +w0*w1)
        product = _mm_add_ps(product, t1);
    }
    return product;
}

inline __m128 SIMD::QConjugate(__m128 const q)
{
    __m128 conjugate = _mm_mul_ps(SIMD::MMMP, q);
    return conjugate;
}

inline __m128 SIMD::QInverse(__m128 const q)
{
    // (-x, -y, -z, +w)
    __m128 conjugate = _mm_mul_ps(SIMD::MMMP, q);

    // (sqrlen, sqrlen, sqrlen, sqrlen)
    __m128 sqrlen = SIMD::Dot(conjugate, conjugate);

    // Divide by the squared length.  This is potentially a divide by
    // zero or a divide by infinity.
    __m128 inverse = _mm_div_ps(conjugate, sqrlen);

    // Set to zero when the squared length is zero.
    __m128 mask = _mm_cmpneq_ps(SIMD::ZZZZ, sqrlen);
    inverse = _mm_and_ps(mask, inverse);
    return inverse;
}

inline __m128 SIMD::QSlerp(__m128 const t, __m128 const q0, __m128 const q1)
{
    float const onePlusMuFPU = 1.90110745351730037f;

    __m128 cs = Dot(q0, q1);
    __m128 negative = _mm_cmplt_ps(cs, ZZZZ);
    __m128 term0 = _mm_and_ps(negative, MMMM);
    __m128 term1 = _mm_andnot_ps(negative, PPPP);
    __m128 sign = _mm_or_ps(term0, term1);
    cs = _mm_mul_ps(cs, sign);

    __m128 csm1 = _mm_sub_ps(cs, PPPP);
    __m128 omt = _mm_sub_ps(PPPP, t);
    // (1-t,1-t,t,t)
    __m128 temp = _mm_shuffle_ps(omt, t, _MM_SHUFFLE(0, 0, 0, 0));
    // (1-t,t,0,0)
    __m128 coeff = _mm_shuffle_ps(temp, ZZZZ, _MM_SHUFFLE(0, 0, 2, 0));
    __m128 u = coeff;
    __m128 sqr = _mm_mul_ps(coeff, u);

    __m128 avalue = _mm_set1_ps(1.0f / (1.0f*3.0f));
    __m128 bvalue = _mm_set1_ps(1.0f / 3.0f);
    temp = _mm_mul_ps(avalue, sqr);
    temp = _mm_sub_ps(temp, bvalue);
    temp = _mm_mul_ps(temp, csm1);
    coeff = _mm_mul_ps(coeff, temp);
    u = _mm_add_ps(u, coeff);

    avalue = _mm_set1_ps(1.0f / (2.0f*5.0f));
    bvalue = _mm_set1_ps(2.0f / 5.0f);
    temp = _mm_mul_ps(avalue, sqr);
    temp = _mm_sub_ps(temp, bvalue);
    temp = _mm_mul_ps(temp, csm1);
    coeff = _mm_mul_ps(coeff, temp);
    u = _mm_add_ps(u, coeff);

    avalue = _mm_set1_ps(1.0f / (3.0f*7.0f));
    bvalue = _mm_set1_ps(3.0f / 7.0f);
    temp = _mm_mul_ps(avalue, sqr);
    temp = _mm_sub_ps(temp, bvalue);
    temp = _mm_mul_ps(temp, csm1);
    coeff = _mm_mul_ps(coeff, temp);
    u = _mm_add_ps(u, coeff);

    avalue = _mm_set1_ps(1.0f / (4.0f*9.0f));
    bvalue = _mm_set1_ps(4.0f / 9.0f);
    temp = _mm_mul_ps(avalue, sqr);
    temp = _mm_sub_ps(temp, bvalue);
    temp = _mm_mul_ps(temp, csm1);
    coeff = _mm_mul_ps(coeff, temp);
    u = _mm_add_ps(u, coeff);

    avalue = _mm_set1_ps(1.0f / (5.0f*11.0f));
    bvalue = _mm_set1_ps(5.0f / 11.0f);
    temp = _mm_mul_ps(avalue, sqr);
    temp = _mm_sub_ps(temp, bvalue);
    temp = _mm_mul_ps(temp, csm1);
    coeff = _mm_mul_ps(coeff, temp);
    u = _mm_add_ps(u, coeff);

    avalue = _mm_set1_ps(1.0f / (6.0f*13.0f));
    bvalue = _mm_set1_ps(6.0f / 13.0f);
    temp = _mm_mul_ps(avalue, sqr);
    temp = _mm_sub_ps(temp, bvalue);
    temp = _mm_mul_ps(temp, csm1);
    coeff = _mm_mul_ps(coeff, temp);
    u = _mm_add_ps(u, coeff);

    avalue = _mm_set1_ps(1.0f / (7.0f*15.0f));
    bvalue = _mm_set1_ps(7.0f / 15.0f);
    temp = _mm_mul_ps(avalue, sqr);
    temp = _mm_sub_ps(temp, bvalue);
    temp = _mm_mul_ps(temp, csm1);
    coeff = _mm_mul_ps(coeff, temp);
    u = _mm_add_ps(u, coeff);

    avalue = _mm_set1_ps(1.0f / (8.0f*17.0f));
    bvalue = _mm_set1_ps(8.0f / 17.0f);
    temp = _mm_mul_ps(avalue, sqr);
    temp = _mm_sub_ps(temp, bvalue);
    temp = _mm_mul_ps(temp, csm1);
    coeff = _mm_mul_ps(coeff, temp);
    u = _mm_add_ps(u, coeff);

    avalue = _mm_set1_ps(onePlusMuFPU*1.0f / (9.0f*19.0f));
    bvalue = _mm_set1_ps(onePlusMuFPU*9.0f / 19.0f);
    temp = _mm_mul_ps(avalue, sqr);
    temp = _mm_sub_ps(temp, bvalue);
    temp = _mm_mul_ps(temp, csm1);
    coeff = _mm_mul_ps(coeff, temp);
    u = _mm_add_ps(u, coeff);

    term0 = _mm_shuffle_ps(u, u, _MM_SHUFFLE(0, 0, 0, 0));
    term1 = _mm_shuffle_ps(u, u, _MM_SHUFFLE(1, 1, 1, 1));
    term0 = _mm_mul_ps(term0, q0);
    term1 = _mm_mul_ps(term1, q1);
    term1 = _mm_mul_ps(term1, sign);
    __m128 slerp = _mm_add_ps(term0, term1);
    return slerp;
}



// Function evaluations (generally slow, function call per component).

inline __m128 SIMD::Sin(__m128 const v)
{
    __m128 result;
    result.m128_f32[0] = sin(v.m128_f32[0]);
    result.m128_f32[1] = sin(v.m128_f32[1]);
    result.m128_f32[2] = sin(v.m128_f32[2]);
    result.m128_f32[3] = sin(v.m128_f32[3]);
    return result;
}

inline __m128 SIMD::Cos(__m128 const v)
{
    __m128 result;
    result.m128_f32[0] = cos(v.m128_f32[0]);
    result.m128_f32[1] = cos(v.m128_f32[1]);
    result.m128_f32[2] = cos(v.m128_f32[2]);
    result.m128_f32[3] = cos(v.m128_f32[3]);
    return result;
}

inline __m128 SIMD::Tan(__m128 const v)
{
    __m128 result;
    result.m128_f32[0] = tan(v.m128_f32[0]);
    result.m128_f32[1] = tan(v.m128_f32[1]);
    result.m128_f32[2] = tan(v.m128_f32[2]);
    result.m128_f32[3] = tan(v.m128_f32[3]);
    return result;
}

inline __m128 SIMD::ASin(__m128 const v)
{
    __m128 result;
    result.m128_f32[0] = asin(v.m128_f32[0]);
    result.m128_f32[1] = asin(v.m128_f32[1]);
    result.m128_f32[2] = asin(v.m128_f32[2]);
    result.m128_f32[3] = asin(v.m128_f32[3]);
    return result;
}

inline __m128 SIMD::ACos(__m128 const v)
{
    __m128 result;
    result.m128_f32[0] = acos(v.m128_f32[0]);
    result.m128_f32[1] = acos(v.m128_f32[1]);
    result.m128_f32[2] = acos(v.m128_f32[2]);
    result.m128_f32[3] = acos(v.m128_f32[3]);
    return result;
}

inline __m128 SIMD::ATan(__m128 const v)
{
    __m128 result;
    result.m128_f32[0] = atan(v.m128_f32[0]);
    result.m128_f32[1] = atan(v.m128_f32[1]);
    result.m128_f32[2] = atan(v.m128_f32[2]);
    result.m128_f32[3] = atan(v.m128_f32[3]);
    return result;
}



// Fast function approximations.

inline void SIMD::ReduceAnglesSin(__m128 const x, __m128& y)
{
    // Map x to y in [-pi,pi], x = 2*pi*quotient + remainder.
    __m128 quotient = _mm_mul_ps(x, INV_TWO_PI);
    quotient = Round(quotient);
    y = _mm_mul_ps(quotient, TWO_PI);
    y = _mm_sub_ps(x, y);

    // Map y to [-pi/2,pi/2] with sin(y) = sin(x).
    __m128 sign = _mm_and_ps(x, SIGN);
    __m128 c = _mm_or_ps(PI, sign);  // pi when x >= 0, -pi when x < 0
    __m128 absx = _mm_andnot_ps(sign, x);  // |x|
    __m128 rflx = _mm_sub_ps(c, x);
    __m128 comp = _mm_cmple_ps(absx, HALF_PI);
    __m128 select0 = _mm_and_ps(comp, x);
    __m128 select1 = _mm_andnot_ps(comp, rflx);
    y = _mm_or_ps(select0, select1);
}

inline __m128 SIMD::SinApprDeg11(__m128 const x)
{
    __m128 xsqr = _mm_mul_ps(x, x);
    __m128 poly = _mm_mul_ps(C_SIN_APPR_DEG11_5, xsqr);
    poly = _mm_add_ps(poly, C_SIN_APPR_DEG11_4);
    poly = _mm_mul_ps(poly, xsqr);
    poly = _mm_add_ps(poly, C_SIN_APPR_DEG11_3);
    poly = _mm_mul_ps(poly, xsqr);
    poly = _mm_add_ps(poly, C_SIN_APPR_DEG11_2);
    poly = _mm_mul_ps(poly, xsqr);
    poly = _mm_add_ps(poly, C_SIN_APPR_DEG11_1);
    poly = _mm_mul_ps(poly, xsqr);
    poly = _mm_add_ps(poly, C_SIN_APPR_DEG11_0);
    poly = _mm_mul_ps(poly, x);
    return poly;
}

inline __m128 SIMD::SinApprDeg7(__m128 const x)
{
    __m128 xsqr = _mm_mul_ps(x, x);
    __m128 poly = _mm_mul_ps(C_SIN_APPR_DEG7_3, xsqr);
    poly = _mm_add_ps(poly, C_SIN_APPR_DEG7_2);
    poly = _mm_mul_ps(poly, xsqr);
    poly = _mm_add_ps(poly, C_SIN_APPR_DEG7_1);
    poly = _mm_mul_ps(poly, xsqr);
    poly = _mm_add_ps(poly, C_SIN_APPR_DEG7_0);
    poly = _mm_mul_ps(poly, x);
    return poly;
}

inline void SIMD::ReduceAnglesCos(__m128 const x, __m128& y, __m128& sign)
{
    // Map x to y in [-pi,pi], x = 2*pi*quotient + remainder.
    __m128 quotient = _mm_mul_ps(x, INV_TWO_PI);
    quotient = Round(quotient);
    y = _mm_mul_ps(quotient, TWO_PI);
    y = _mm_sub_ps(x, y);

    // Map x to y in [-pi/2,pi/2] with cos(y) = sign*cos(x).
    sign = _mm_and_ps(x, SIGN);
    __m128 c = _mm_or_ps(PI, sign);  // pi when x >= 0, -pi when x < 0
    __m128 absx = _mm_andnot_ps(sign, x);  // |x|
    __m128 rflx = _mm_sub_ps(c, x);
    __m128 comp = _mm_cmple_ps(absx, HALF_PI);
    __m128 select0 = _mm_and_ps(comp, x);
    __m128 select1 = _mm_andnot_ps(comp, rflx);
    y = _mm_or_ps(select0, select1);
    select0 = _mm_and_ps(comp, PPPP);
    select1 = _mm_andnot_ps(comp, MMMM);
    sign = _mm_or_ps(select0, select1);
}

inline __m128 SIMD::CosApprDeg10(__m128 const x, __m128 const sign)
{
    __m128 xsqr = _mm_mul_ps(x, x);
    __m128 poly = _mm_mul_ps(C_COS_APPR_DEG10_5, xsqr);
    poly = _mm_add_ps(poly, C_COS_APPR_DEG10_4);
    poly = _mm_mul_ps(poly, xsqr);
    poly = _mm_add_ps(poly, C_COS_APPR_DEG10_3);
    poly = _mm_mul_ps(poly, xsqr);
    poly = _mm_add_ps(poly, C_COS_APPR_DEG10_2);
    poly = _mm_mul_ps(poly, xsqr);
    poly = _mm_add_ps(poly, C_COS_APPR_DEG10_1);
    poly = _mm_mul_ps(poly, xsqr);
    poly = _mm_add_ps(poly, C_COS_APPR_DEG10_0);
    poly = _mm_mul_ps(poly, sign);
    return poly;
}

inline __m128 SIMD::CosApprDeg6(__m128 const x, __m128 const sign)
{
    __m128 xsqr = _mm_mul_ps(x, x);
    __m128 poly = _mm_mul_ps(C_COS_APPR_DEG6_3, xsqr);
    poly = _mm_add_ps(poly, C_COS_APPR_DEG6_2);
    poly = _mm_mul_ps(poly, xsqr);
    poly = _mm_add_ps(poly, C_COS_APPR_DEG6_1);
    poly = _mm_mul_ps(poly, xsqr);
    poly = _mm_add_ps(poly, C_COS_APPR_DEG6_0);
    poly = _mm_mul_ps(poly, sign);
    return poly;
}

inline void SIMD::GetAdjDet(__m128 const* mat, __m128* adj, __m128* det)
{
    // [GTE_USE_MAT_VEC]
    // a0 = m00*m11 - m01*m10, b0 = m20*m31 - m21*m30
    // a1 - m00*m12 - m02*m10, b1 = m20*m32 - m22*m30
    // a2 = m00*m13 - m03*m10, b2 = m20*m33 - m23*m30
    // a3 = m01*m12 - m02*m11, b3 = m21*m32 - m22*m31
    // a4 = m01*m13 - m03*m11, b4 = m21*m33 - m23*m31
    // a5 = m02*m13 - m03*m12, b5 = m22*m33 - m23*m32
    // +c00 = m11*b5 - m12*b4 + m13*b3
    // -c10 = m10*b5 - m12*b2 + m13*b1
    // +c20 = m10*b4 - m11*b2 + m13*b0
    // -c30 = m10*b3 - m11*b1 + m12*b0
    // -c01 = m01*b5 - m02*b4 + m03*b3
    // +c11 = m00*b5 - m02*b2 + m03*b1
    // -c21 = m00*b4 - m01*b2 + m03*b0
    // +c31 = m00*b3 - m01*b1 + m02*b0
    // +c02 = m31*a5 - m32*a4 + m33*a3
    // -c12 = m30*a5 - m32*a2 + m33*a1
    // +c22 = m30*a4 - m31*a2 + m33*a0
    // -c32 = m30*a3 - m31*a1 + m32*a0
    // -c03 = m21*a5 - m22*a4 + m23*a3
    // +c13 = m20*a5 - m22*a2 + m23*a1
    // -c23 = m20*a4 - m21*a2 + m23*a0
    // +c33 = m20*a3 - m21*a1 + m22*a0
    //
    // [GTE_USE_VEC_MAT]
    // a0 = m00*m11 - m01*m10,  b0 = m02*m13 - m03*m12
    // a1 = m00*m21 - m01*m20,  b1 = m02*m23 - m03*m22
    // a2 = m00*m31 - m01*m30,  b2 = m02*m33 - m03*m32
    // a3 = m10*m21 - m11*m20,  b3 = m12*m23 - m13*m22
    // a4 = m10*m31 - m11*m30,  b4 = m12*m33 - m13*m32
    // a5 = m20*m31 - m21*m30,  b5 = m22*m33 - m23*m32
    // +c00 = m11*b5 - m21*b4 + m31*b3
    // -c01 = m01*b5 - m21*b2 + m31*b1
    // +c02 = m01*b4 - m11*b2 + m31*b0
    // -c03 = m01*b3 - m11*b1 + m21*b0
    // -c10 = m10*b5 - m20*b4 + m30*b3
    // +c11 = m00*b5 - m20*b2 + m30*b1
    // -c12 = m00*b4 - m10*b2 + m30*b0
    // +c13 = m00*b3 - m10*b1 + m20*b0
    // +c20 = m13*a5 - m23*a4 + m33*a3
    // -c21 = m03*a5 - m23*a2 + m33*a1
    // +c22 = m03*a4 - m13*a2 + m33*a0
    // -c23 = m03*a3 - m13*a1 + m23*a0
    // -c30 = m12*a5 - m22*a4 + m32*a3
    // +c31 = m02*a5 - m22*a2 + m32*a1
    // -c32 = m02*a4 - m12*a2 + m32*a0
    // +c33 = m02*a3 - m12*a1 + m22*a0
    //
    // det = a0*b5 - a1*b4 + a2*b3 + a3*b2 - a4*b1 + a5*b0
    // inverse[i][j] = c[i][j]/det
    __m128 t0, t1;

    // Compute a1, a2, a3, a4.
    __m128 a1a2a3a4;
    {
        // MV:(m00, m00, m01, m01), VM:(m00, m00, m10, m10)
        t0 = _mm_shuffle_ps(mat[0], mat[1], _MM_SHUFFLE(0, 0, 0, 0));
        // MV:(m12, m12, m13, m13), VM:(m21, m21, m31, m31)
        t1 = _mm_shuffle_ps(mat[2], mat[3], _MM_SHUFFLE(1, 1, 1, 1));
        // MV:(m12, m13, m12, m13), VM:(m21, m31, m21, m31)
        t1 = _mm_shuffle_ps(t1, t1, _MM_SHUFFLE(2, 0, 2, 0));
        // MV:(m00*m12, m00*m13, m01*m12, m01*m13)
        // VM:(m00*m21, m00*m31, m10*m21, m10*m31)
        a1a2a3a4 = _mm_mul_ps(t0, t1);
        // MV:(m10, m10, m11, m11), VM:(m01, m01, m11, m11)
        t0 = _mm_shuffle_ps(mat[0], mat[1], _MM_SHUFFLE(1, 1, 1, 1));
        // MV:(m02, m02, m03, m03), VM:(m20, m20, m30, m30)
        t1 = _mm_shuffle_ps(mat[2], mat[3], _MM_SHUFFLE(0, 0, 0, 0));
        // MV:(m02, m03, m02, m03), VM:(m20, m30, m20, m30)
        t1 = _mm_shuffle_ps(t1, t1, _MM_SHUFFLE(2, 0, 2, 0));
        // MV:(m10*m02, m10*m03, m11*m02, m11*m03)
        // VM:(m01*m20, m01*m30, m11*m20, m11*m30)
        t0 = _mm_mul_ps(t0, t1);
        // MV:(m00*m12-m10*m02,m00*m13-m10*m03,m01*m12-m11*m02,m01*m13-m11*m03)
        // VM:(m00*m21-m01*m20,m00*m31-m01*m30,m10*m21-m11*m20,m10*m31-m11*m30)
        a1a2a3a4 = _mm_sub_ps(a1a2a3a4, t0);
    }

    // Compute b1, b2, b3, b4.
    __m128 b1b2b3b4;
    {
        // MV:(m20, m20, m21, m21), VM:(m02, m02, m12, m12)
        t0 = _mm_shuffle_ps(mat[0], mat[1], _MM_SHUFFLE(2, 2, 2, 2));
        // MV:(m32, m32, m33, m33), VM:(m23, m23, m33, m33)
        t1 = _mm_shuffle_ps(mat[2], mat[3], _MM_SHUFFLE(3, 3, 3, 3));
        // MV:(m32, m33, m32, m33), VM:(m23, m33, m23, m33)
        t1 = _mm_shuffle_ps(t1, t1, _MM_SHUFFLE(2, 0, 2, 0));
        // MV:(m20*m32, m20*m33, m21*m32, m21*m33)
        // VM:(m02*m23, m02*m33, m12*m23, m12*m33)
        b1b2b3b4 = _mm_mul_ps(t0, t1);
        // MV:(m30, m30, m31, m31), VM:(m03, m03, m13, m13)
        t0 = _mm_shuffle_ps(mat[0], mat[1], _MM_SHUFFLE(3, 3, 3, 3));
        // MV:(m22, m22, m23, m23), VM:(m22, m22, m32, m32)
        t1 = _mm_shuffle_ps(mat[2], mat[3], _MM_SHUFFLE(2, 2, 2, 2));
        // MV:(m22, m23, m22, m23), VM:(m22, m32, m22, m32)
        t1 = _mm_shuffle_ps(t1, t1, _MM_SHUFFLE(2, 0, 2, 0));
        // MV:(m30*m22, m30*m23, m31*m22, m31*m23)
        // VM:(m03*m22, m03*m32, m13*m22, m13*m32)
        t0 = _mm_mul_ps(t0, t1);
        // MV:(m20*m32-m30*m22,m20*m33-m30*m23,m21*m32-m31*m22,m21*m33-m31*m22)
        // VM:(m02*m23-m03*m22,m02*m33-m03*m32,m12*m23-m13*m22,m12*m33-m13*m32)
        b1b2b3b4 = _mm_sub_ps(b1b2b3b4, t0);
    }

    // Compute a0, b0, a5, b5.
    __m128 a0b0a5b5;
    {
        // MV:(m00, m20, m02, m22), VM:(m00, m02, m20, m22)
        t0 = _mm_shuffle_ps(mat[0], mat[2], _MM_SHUFFLE(2, 0, 2, 0));
        // MV:(m11, m31, m13, m33), VM:(m11, m13, m31, m33)
        t1 = _mm_shuffle_ps(mat[1], mat[3], _MM_SHUFFLE(3, 1, 3, 1));
        // MV:(m00*m11, m20*m31, m02*m13, m22*m33)
        // VM:(m00*m11, m02*m13, m20*m31, m22*m33)
        a0b0a5b5 = _mm_mul_ps(t0, t1);
        // MV:(m10, m30, m12, m32), VM:(m01, m03, m21, m23)
        t0 = _mm_shuffle_ps(mat[0], mat[2], _MM_SHUFFLE(3, 1, 3, 1));
        // MV:(m01, m21, m03, m23), VM:(m10, m12, m30, m32)
        t1 = _mm_shuffle_ps(mat[1], mat[3], _MM_SHUFFLE(2, 0, 2, 0));
        // MV:(m10*m01, m30*m21, m12*m03, m32*m23)
        // VM:(m01*m10, m03*m12, m21*m30, m23*m32)
        t0 = _mm_mul_ps(t0, t1);
        // MV:(m00*m11-m10*m01,m20*m31-m30*m21,m02*m13-m12*m03,m22*m33-m32*m23)
        // VM:(m00*m11-m01*m10,m02*m13-m03*m12,m20*m31-m21*m30,m22*m33-m23*m32)
        a0b0a5b5 = _mm_sub_ps(a0b0a5b5, t0);
    }

    if (adj)
    {
        // Compute slices 0 and 1 of the adjoint matrix.  An MV slice is a
        // column and a VM slice is a row.
        __m128 slice0, slice1;
        {
            __m128 b5b5b4b3 = _mm_shuffle_ps(a0b0a5b5, b1b2b3b4, _MM_SHUFFLE(2, 3, 3, 3));
            __m128 b4b2b2b1 = _mm_shuffle_ps(b1b2b3b4, b1b2b3b4, _MM_SHUFFLE(0, 1, 1, 3));
            __m128 b3b1b0b0 = _mm_shuffle_ps(b1b2b3b4, a0b0a5b5, _MM_SHUFFLE(1, 1, 0, 2));

            // Compute slice 0 of the adjoint matrix.
            {
                // MV:(m11, m11, m10, m10), VM:(m11, m11, m01, m01)
                t0 = _mm_shuffle_ps(mat[1], mat[0], _MM_SHUFFLE(1, 1, 1, 1));
                // MV:(m11, m10, m10, m10), VM:(m11, m01, m01, m01)
                t1 = _mm_shuffle_ps(t0, t0, _MM_SHUFFLE(2, 2, 2, 0));
                // MV:(m11*b5, m10*b5, m10*b4, m10*b3)
                // VM:(m11*b5, m01*b5, m01*b4, m01*b3)
                slice0 = _mm_mul_ps(t1, b5b5b4b3);

                // MV:(m12, m12, m11, m11), VM:(m21, m21, m11, m11)
                t0 = _mm_shuffle_ps(mat[2], mat[1], _MM_SHUFFLE(1, 1, 1, 1));
                // MV:(m12*b4, m12*b2, m11*b2, m11*b1)
                // VM:(m21*b4, m21*b2, m11*b2, m11*b1)
                t1 = _mm_mul_ps(t0, b4b2b2b1);
                slice0 = _mm_sub_ps(slice0, t1);

                // MV:(m13, m13, m12, m12), VM:(m31, m31, m21, m21)
                t0 = _mm_shuffle_ps(mat[3], mat[2], _MM_SHUFFLE(1, 1, 1, 1));
                // MV:(m13, m13, m13, m12), VM:(m31, m31, m31, m21)
                t1 = _mm_shuffle_ps(t0, t0, _MM_SHUFFLE(2, 0, 0, 0));
                // MV:(m13*b3, m13*b1, m13*b0, m12*b0)
                // VM:(m31*b3, m31*b1, m31*b0, m21*b0)
                t0 = _mm_mul_ps(t1, b3b1b0b0);
                slice0 = _mm_add_ps(slice0, t0);

                // MV:(c00, c10, c20, c30), VM:(c00, c01, c02, c03)
                slice0 = _mm_mul_ps(slice0, PMPM);
            }

            // Compute slice 1 of the adjoint matrix.
            {
                // MV:(m01, m01, m00, m00), VM:(m10, m10, m00, m00)
                t0 = _mm_shuffle_ps(mat[1], mat[0], _MM_SHUFFLE(0, 0, 0, 0));
                // MV:(m01, m00, m00, m00), VM:(m10, m00, m00, m00)
                t1 = _mm_shuffle_ps(t0, t0, _MM_SHUFFLE(2, 2, 2, 0));
                // MV:(m01*b5, m00*b5, m00*b4, m00*b3)
                // VM:(m10*b5, m00*b5, m00*b4, m00*b3)
                slice1 = _mm_mul_ps(t1, b5b5b4b3);

                // MV:(m02, m02, m01, m01), VM:(m20, m20, m10, m10)
                t0 = _mm_shuffle_ps(mat[2], mat[1], _MM_SHUFFLE(0, 0, 0, 0));
                // MV:(m02*b4, m02*b2, m01*b2, m01*b1)
                // VM:(m20*b4, m20*b2, m10*b2, m10*b1)
                t1 = _mm_mul_ps(t0, b4b2b2b1);
                slice1 = _mm_sub_ps(slice1, t1);

                // MV:(m03, m03, m02, m02), VM:(m30, m30, m20, m20)
                t0 = _mm_shuffle_ps(mat[3], mat[2], _MM_SHUFFLE(0, 0, 0, 0));
                // MV:(m03, m03, m03, m02), VM:(m30, m30, m30, m20)
                t1 = _mm_shuffle_ps(t0, t0, _MM_SHUFFLE(2, 0, 0, 0));
                // MV:(m03*b3, m03*b1, m03*b0, m02*b0)
                // VM:(m30*b3, m30*b1, m30*b0, m20*b0)
                t0 = _mm_mul_ps(t1, b3b1b0b0);
                slice1 = _mm_add_ps(slice1, t0);

                // MV:(c01, c11, c21, c31), VM:(c10, c11, c12, c13)
                slice1 = _mm_mul_ps(slice1, MPMP);
            }
        }

        // Compute slices 2 and 3 of the adjoint matrix.  An MV slice is a
        // column and a VM slice is a row.
        __m128 slice2, slice3;
        {
            __m128 a5a5a4a3 = _mm_shuffle_ps(a0b0a5b5, a1a2a3a4, _MM_SHUFFLE(2, 3, 2, 2));
            __m128 a4a2a2a1 = _mm_shuffle_ps(a1a2a3a4, a1a2a3a4, _MM_SHUFFLE(0, 1, 1, 3));
            __m128 a3a1a0a0 = _mm_shuffle_ps(a1a2a3a4, a0b0a5b5, _MM_SHUFFLE(0, 0, 0, 2));

            // Compute slice 2 of the adjoint matrix.
            {
                // MV:(m31, m31, m30, m30), VM:(m13, m13, m03, m03)
                t0 = _mm_shuffle_ps(mat[1], mat[0], _MM_SHUFFLE(3, 3, 3, 3));
                // MV:(m13, m03, m03, m03), VM:(m13, m03, m03, m03)
                t1 = _mm_shuffle_ps(t0, t0, _MM_SHUFFLE(2, 2, 2, 0));
                // MV:(m31*a5, m30*a5, m30*a4, m30*a3)
                // VM:(m13*a5, m03*a5, m03*a4, m03*a3)
                slice2 = _mm_mul_ps(t1, a5a5a4a3);

                // MV:(m32, m32, m31, m31), VM:(m23, m23, m13, m13)
                t0 = _mm_shuffle_ps(mat[2], mat[1], _MM_SHUFFLE(3, 3, 3, 3));
                // MV:(m32*a4, m32*a2, m31*a2, m31*a1)
                // VM:(m23*a4, m23*a2, m13*a2, m13*a1)
                t1 = _mm_mul_ps(t0, a4a2a2a1);
                slice2 = _mm_sub_ps(slice2, t1);

                // MV:(m33, m33, m32, m32), VM:(m33, m33, m23, m23)
                t0 = _mm_shuffle_ps(mat[3], mat[2], _MM_SHUFFLE(3, 3, 3, 3));
                // MV:(m33, m33, m33, m32), VM:(m33, m33, m33, m23)
                t1 = _mm_shuffle_ps(t0, t0, _MM_SHUFFLE(2, 0, 0, 0));
                // MV:(m33*a3, m33*a1, m33*a0, m32*a0)
                // VM:(m33*a3, m33*a1, m33*a0, m23*a0)
                t0 = _mm_mul_ps(t1, a3a1a0a0);
                slice2 = _mm_add_ps(slice2, t0);

                // MV:(c02, c12, c22, c32), VM:(c20, c21, c22, c23)
                slice2 = _mm_mul_ps(slice2, PMPM);
            }

            // Compute slice 3 of the adjoint matrix.
            {
                // MV:(m21, m21, m20, m20), VM:(m12, m12, m02, m02)
                t0 = _mm_shuffle_ps(mat[1], mat[0], _MM_SHUFFLE(2, 2, 2, 2));
                // MV:(m21, m20, m20, m20), VM:(m12, m02, m02, m02)
                t1 = _mm_shuffle_ps(t0, t0, _MM_SHUFFLE(2, 2, 2, 0));
                // MV:(m21*a5, m20*a5, m20*a4, m20*a3)
                // VM:(m12*a5, m02*a5, m02*a4, m02*a3)
                slice3 = _mm_mul_ps(t1, a5a5a4a3);

                // MV:(m22, m22, m21, m21), VM:(m22, m22, m12, m12)
                t0 = _mm_shuffle_ps(mat[2], mat[1], _MM_SHUFFLE(2, 2, 2, 2));
                // MV:(m22*a4, m22*a2, m21*a2, m21*a1)
                // VM:(m22*a4, m22*a2, m12*a2, m12*a1)
                t1 = _mm_mul_ps(t0, a4a2a2a1);
                slice3 = _mm_sub_ps(slice3, t1);

                // MV:(m23, m23, m22, m22), VM:(m32, m32, m22, m22)
                t0 = _mm_shuffle_ps(mat[3], mat[2], _MM_SHUFFLE(2, 2, 2, 2));
                // MV:(m23, m23, m23, m22), VM:(m32, m32, m32, m22)
                t1 = _mm_shuffle_ps(t0, t0, _MM_SHUFFLE(2, 0, 0, 0));
                // MV:(m23*a3, m23*a1, m23*a0, m22*a0)
                // VM:(m32*a3, m32*a1, m32*a0, m22*a0)
                t0 = _mm_mul_ps(t1, a3a1a0a0);
                slice3 = _mm_add_ps(slice3, t0);

                // MV:(c03, c13, c23, c33), VM:(c30, c31, c32, c33)
                slice3 = _mm_mul_ps(slice3, MPMP);
            }
        }

        adj[0] = slice0;
        adj[1] = slice1;
        adj[2] = slice2;
        adj[3] = slice3;

        if (det)
        {
            // Compute the determinant using the cofactors.
            // det = m00*c00 + m01*c10 + m02*c20 + m03*c30
            {
                // MV:(c00, c00, c01, c01), VM:(c00, c00, c10, c10)
                t0 = _mm_shuffle_ps(slice0, slice1, _MM_SHUFFLE(0, 0, 0, 0));
                // MV:(c02, c02, c03, c03), VM:(c20, c20, c30, c30)
                t1 = _mm_shuffle_ps(slice2, slice3, _MM_SHUFFLE(0, 0, 0, 0));
                // MV:(c00, c01, c02, c03), VM:(c00, c10, c20, c30)
                t1 = _mm_shuffle_ps(t0, t1, _MM_SHUFFLE(2, 0, 2, 0));
                // MV:(m00*c00, m10*c01, m20*c02, m30*c03)
                // VM:(m00*c00, m01*c10, m02*c20, m03*c30)
                t0 = _mm_mul_ps(mat[0], t1);
                // MV:(m10*c01, m00*c00, m30*c03, m20*c02)
                // VM:(m01*c10, m00*c00, m03*c30, m02*c20)
                t1 = _mm_shuffle_ps(t0, t0, _MM_SHUFFLE(2, 3, 0, 1));
                // MV:(m00*c00+m10*c01,m00*c00+m01*c01,m20*c02+m30*c03,m20*c02+m30*c03)
                // VM:(m00*c00+m01*c10,m00*c00+m10*c10,m02*c20+m03*c30,m02*c20+m03*c30)
                t0 = _mm_add_ps(t0, t1);
                // MV:(m20*c02+m30*c03,m20*c02+m30*c03,m00*c00+m10*c01,m00*c00+m10*c01)
                // VM:(m02*c20+m03*c30,m02*c20+m03*c30,m00*c00+m01*c10,m00*c00+m01*c10)
                t1 = _mm_shuffle_ps(t0, t0, _MM_SHUFFLE(0, 0, 2, 2));
                // (det, det, det, det)
                *det = _mm_add_ps(t0, t1);
            }
        }
    }
    else if (det)
    {
        // Compute the determinant using the a- and b-coefficients.
        // det = (a0*b5+a5*b0+a2*b3+a3*b2) - (a1*b4+a4*b1) = dot0 - dot1

        // (a0, a5, a2, a3)
        t0 = _mm_shuffle_ps(a0b0a5b5, a1a2a3a4, _MM_SHUFFLE(2, 1, 2, 0));
        // (b5, b0, b3, b2)
        t1 = _mm_shuffle_ps(a0b0a5b5, b1b2b3b4, _MM_SHUFFLE(1, 2, 1, 3));
        // (a0*b5, a5*b0, a2*b3, a3*b2)
        t0 = _mm_mul_ps(t0, t1);
        // (a5*b0, a0*b5, a3*b2, a2*b3)
        t1 = _mm_shuffle_ps(t0, t0, _MM_SHUFFLE(2, 3, 0, 1));
        // (a0*b5+a5*b0, a0*b5+a5*b0, a2*b3+a3*b2, a2*b3+a3*b2)
        t0 = _mm_add_ps(t0, t1);
        // (a2*b3+a3*b2, a2*b3+a3*b2, a0*b5+a5*b0, a0*b5+a5*b0)
        t1 = _mm_shuffle_ps(t0, t0, _MM_SHUFFLE(0, 0, 2, 2));
        // (dot0, dot0, dot0, dot0)
        __m128 dot0 = _mm_add_ps(t0, t1);

        // (a1, a4, a1, a4)
        t0 = _mm_shuffle_ps(a1a2a3a4, a1a2a3a4, _MM_SHUFFLE(3, 0, 3, 0));
        // (b4, b1, b4, b1)
        t1 = _mm_shuffle_ps(b1b2b3b4, b1b2b3b4, _MM_SHUFFLE(0, 3, 0, 3));
        // (a1*b4, a4*b1, a1*b4, a4*b1)
        t0 = _mm_mul_ps(t0, t1);
        // (a4*b1, a1*b4, a4*b1, a1*b4)
        t1 = _mm_shuffle_ps(t0, t0, _MM_SHUFFLE(0, 1, 0, 1));
        // (a1*b4+a4*b1, a1*b4+a4*b1, a1*b4+a4*b1, a1*b4+a4*b1)
        __m128 dot1 = _mm_add_ps(t0, t1);
        *det = _mm_sub_ps(dot0, dot1);
    }
}


}
