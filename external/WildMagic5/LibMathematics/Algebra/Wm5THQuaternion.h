// Geometric Tools, LLC
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.16.0 (2017/08/24)

#ifndef WM5THQUATERNION_H
#define WM5THQUATERNION_H

#include "Wm5MathematicsLIB.h"
#include "Wm5THMatrix.h"

namespace Wm5
{

template <typename Real>
class THQuaternion
{
public:
    // A quaternion is q = w + x*i + y*j + z*k where (w,x,y,z) is not
    // necessarily a unit-length vector in 4D.

    // Construction.
    THQuaternion()
    {
        // uninitialized members
    }

    THQuaternion(Real w, Real x, Real y, Real z)
    {
        mTuple[0] = w;
        mTuple[1] = x;
        mTuple[2] = y;
        mTuple[3] = z;
    }

    THQuaternion(const THQuaternion& q)
    {
        mTuple[0] = q.mTuple[0];
        mTuple[1] = q.mTuple[1];
        mTuple[2] = q.mTuple[2];
        mTuple[3] = q.mTuple[3];
    }

    // THQuaternion for the input rotation matrix.
    THQuaternion(const THMatrix<Real>& rot)
    {
        FromRotationMatrix(rot);
    }

    // THQuaternion for the rotation of the axis-angle pair.
    THQuaternion(const TAVector<Real>& axis, Real angle)
    {
        FromAxisAngle(axis, angle);
    }

    // Coordinate access as an array:  0 = w, 1 = x, 2 = y, 3 = z.
    inline operator const Real* () const
    {
        return mTuple;
    }

    inline operator Real* ()
    {
        return mTuple;
    }

    inline const Real& operator[] (int i) const
    {
        return mTuple[i];
    }

    inline Real& operator[] (int i)
    {
        return mTuple[i];
    }

    inline Real W() const
    {
        return mTuple[0];
    }

    inline Real& W()
    {
        return mTuple[0];
    }

    inline Real X() const
    {
        return mTuple[1];
    }

    inline Real& X()
    {
        return mTuple[1];
    }

    inline Real Y() const
    {
        return mTuple[2];
    }

    inline Real& Y()
    {
        return mTuple[2];
    }

    inline Real Z() const
    {
        return mTuple[3];
    }

    inline Real& Z()
    {
        return mTuple[3];
    }

    // Assignment.
    THQuaternion& operator= (const THQuaternion& q)
    {
        mTuple[0] = q.mTuple[0];
        mTuple[1] = q.mTuple[1];
        mTuple[2] = q.mTuple[2];
        mTuple[3] = q.mTuple[3];
        return *this;
    }

    // Comparison (for use by STL containers).
    bool operator== (const THQuaternion& q) const
    {
        for (int i = 0; i < 4; ++i)
        {
            if (mTuple[i] != q.mTuple[i])
            {
                return false;
            }
        }
        return true;
    }

    bool operator!= (const THQuaternion& q) const
    {
        return !operator==(q);
    }

    bool operator< (const THQuaternion& q) const
    {
        // lexicographical ordering
        for (int i = 0; i < 4; ++i)
        {
            if (mTuple[i] < q.mTuple[i])
            {
                return true;
            }
            if (mTuple[i] > q.mTuple[i])
            {
                return false;
            }
        }
        return false;
    }

    bool operator<= (const THQuaternion& q) const
    {
        // (x <= y) <=> !(y < x)
        return !(q.operator<(*this));
    }

    bool operator> (const THQuaternion& q) const
    {
        // (x > y) <=> (y < x)
        return q.operator<(*this);
    }

    bool operator>= (const THQuaternion& q) const
    {
        // (x >= y) <=> !(x < y)
        return !operator<(q);
    }

    // Arithmetic operations.
    THQuaternion operator+ (const THQuaternion& q) const
    {
        THQuaternion result;
        for (int i = 0; i < 4; ++i)
        {
            result.mTuple[i] = mTuple[i] + q.mTuple[i];
        }
        return result;
    }

    THQuaternion operator- (const THQuaternion& q) const
    {
        THQuaternion result;
        for (int i = 0; i < 4; ++i)
        {
            result.mTuple[i] = mTuple[i] - q.mTuple[i];
        }
        return result;
    }

    THQuaternion operator* (const THQuaternion& q) const
    {
        // NOTE:  Multiplication is not generally commutative, so in most
        // cases p*q != q*p.

        THQuaternion result;

        result.mTuple[0] =
            mTuple[0] * q.mTuple[0] -
            mTuple[1] * q.mTuple[1] -
            mTuple[2] * q.mTuple[2] -
            mTuple[3] * q.mTuple[3];

        result.mTuple[1] =
            mTuple[0] * q.mTuple[1] +
            mTuple[1] * q.mTuple[0] +
            mTuple[2] * q.mTuple[3] -
            mTuple[3] * q.mTuple[2];

        result.mTuple[2] =
            mTuple[0] * q.mTuple[2] +
            mTuple[2] * q.mTuple[0] +
            mTuple[3] * q.mTuple[1] -
            mTuple[1] * q.mTuple[3];

        result.mTuple[3] =
            mTuple[0] * q.mTuple[3] +
            mTuple[3] * q.mTuple[0] +
            mTuple[1] * q.mTuple[2] -
            mTuple[2] * q.mTuple[1];

        return result;
    }

    THQuaternion operator* (Real scalar) const
    {
        THQuaternion result;
        for (int i = 0; i < 4; ++i)
        {
            result.mTuple[i] = scalar  * mTuple[i];
        }
        return result;
    }

    THQuaternion operator/ (Real scalar) const
    {
        THQuaternion result;
        for (int i = 0; i < 4; ++i)
        {
            result.mTuple[i] = mTuple[i] / scalar;
        }
        return result;
    }

    THQuaternion operator- () const
    {
        THQuaternion result;
        for (int i = 0; i < 4; ++i)
        {
            result.mTuple[i] = -mTuple[i];
        }
        return result;
    }

    // Arithmetic updates.
    THQuaternion& operator+= (const THQuaternion& q)
    {
        for (int i = 0; i < 4; ++i)
        {
            mTuple[i] += q.mTuple[i];
        }
        return *this;
    }

    THQuaternion& operator-= (const THQuaternion& q)
    {
        for (int i = 0; i < 4; ++i)
        {
            mTuple[i] -= q.mTuple[i];
        }
        return *this;
    }

    THQuaternion& operator*= (Real scalar)
    {
        for (int i = 0; i < 4; ++i)
        {
            mTuple[i] *= scalar;
        }
        return *this;
    }

    THQuaternion& operator/= (Real scalar)
    {
        for (int i = 0; i < 4; ++i)
        {
            mTuple[i] /= scalar;
        }
        return *this;
    }

    // Conversion between quaternions, matrices, and axis-angle.
    void FromRotationMatrix(const THMatrix<Real>& rot)
    {
        // Algorithm in Ken Shoemake's article in 1987 SIGGRAPH course notes
        // article "HQuaternion Calculus and Fast Animation".

        const int next[3] = { 1, 2, 0 };

        Real trace = rot(0, 0) + rot(1, 1) + rot(2, 2);
        Real root;

        if (trace > (Real)0)
        {
            // |w| > 1/2, may as well choose w > 1/2
            root = sqrt(trace + (Real)1);  // 2w
            mTuple[0] = ((Real)0.5) * root;
            root = ((Real)0.5) / root;  // 1/(4w)
            mTuple[1] = (rot(2, 1) - rot(1, 2)) * root;
            mTuple[2] = (rot(0, 2) - rot(2, 0)) * root;
            mTuple[3] = (rot(1, 0) - rot(0, 1)) * root;
        }
        else
        {
            // |w| <= 1/2
            int i = 0;
            if (rot(1, 1) > rot(0, 0))
            {
                i = 1;
            }
            if (rot(2, 2) > rot(i, i))
            {
                i = 2;
            }
            int j = next[i];
            int k = next[j];

            root = sqrt(rot(i, i) - rot(j, j) - rot(k, k) + (Real)1);
            Real* quat[3] = { &mTuple[1], &mTuple[2], &mTuple[3] };
            *quat[i] = ((Real)0.5) * root;
            root = ((Real)0.5) / root;
            mTuple[0] = (rot(k, j) - rot(j, k)) * root;
            *quat[j] = (rot(j, i) + rot(i, j)) * root;
            *quat[k] = (rot(k, i) + rot(i, k)) * root;
        }
    }

    void ToRotationMatrix(THMatrix<Real>& rot) const
    {
        Real twoX = ((Real)2) * mTuple[1];
        Real twoY = ((Real)2) * mTuple[2];
        Real twoZ = ((Real)2) * mTuple[3];
        Real twoWX = twoX * mTuple[0];
        Real twoWY = twoY * mTuple[0];
        Real twoWZ = twoZ * mTuple[0];
        Real twoXX = twoX * mTuple[1];
        Real twoXY = twoY * mTuple[1];
        Real twoXZ = twoZ * mTuple[1];
        Real twoYY = twoY * mTuple[2];
        Real twoYZ = twoZ * mTuple[2];
        Real twoZZ = twoZ * mTuple[3];

        rot(0, 0) = (Real)1 - (twoYY + twoZZ);
        rot(0, 1) = twoXY - twoWZ;
        rot(0, 2) = twoXZ + twoWY;
        rot(0, 3) = (Real)0;
        rot(1, 0) = twoXY + twoWZ;
        rot(1, 1) = (Real)1 - (twoXX + twoZZ);
        rot(1, 2) = twoYZ - twoWX;
        rot(1, 3) = (Real)0;
        rot(2, 0) = twoXZ - twoWY;
        rot(2, 1) = twoYZ + twoWX;
        rot(2, 2) = (Real)1 - (twoXX + twoYY);
        rot(2, 3) = (Real)0;
        rot(3, 0) = (Real)0;
        rot(3, 1) = (Real)0;
        rot(3, 2) = (Real)0;
        rot(3, 3) = (Real)1;
    }

    void FromAxisAngle(const TAVector<Real>& axis, Real angle)
    {
        // assert:  axis[] is unit length
        //
        // The quaternion representing the rotation is
        //   q = cos(A/2)+sin(A/2)*(x*i+y*j+z*k)

        Real halfAngle = ((Real)0.5) * angle;
        Real sn = sin(halfAngle);
        mTuple[0] = cos(halfAngle);
        mTuple[1] = sn * axis[0];
        mTuple[2] = sn * axis[1];
        mTuple[3] = sn * axis[2];
    }

    void ToAxisAngle(TAVector<Real>& axis, Real& angle) const
    {
        // The quaternion representing the rotation is
        //   q = cos(A/2)+sin(A/2)*(x*i+y*j+z*k)

        Real sqrLength = mTuple[1] * mTuple[1] + mTuple[2] * mTuple[2]
            + mTuple[3] * mTuple[3];

        if (sqrLength > (Real)0)
        {
            angle = ((Real)2) * acos(mTuple[0]);
            Real invLength = ((Real)1) / sqrt(sqrLength);
            axis[0] = mTuple[1] * invLength;
            axis[1] = mTuple[2] * invLength;
            axis[2] = mTuple[3] * invLength;
        }
        else
        {
            // Angle is 0 (mod 2*pi), so any axis will do.
            angle = (Real)0;
            axis[0] = (Real)1;
            axis[1] = (Real)0;
            axis[2] = (Real)0;
        }
    }

    // Functions of a quaternion.
    Real Length() const  // length of 4-tuple
    {
        return sqrt(mTuple[0] * mTuple[0] + mTuple[1] * mTuple[1] +
            mTuple[2] * mTuple[2] + mTuple[3] * mTuple[3]);
    }

    Real SquaredLength() const  // squared length of 4-tuple
    {
        return mTuple[0] * mTuple[0] + mTuple[1] * mTuple[1] +
            mTuple[2] * mTuple[2] + mTuple[3] * mTuple[3];
    }

    Real Dot(const THQuaternion& q) const  // dot product of 4-tuples
    {
        return mTuple[0] * q.mTuple[0] + mTuple[1] * q.mTuple[1] +
            mTuple[2] * q.mTuple[2] + mTuple[3] * q.mTuple[3];
    }

    Real Normalize(const Real epsilon = (Real)0)
    {
        Real length = Length();

        if (length > epsilon)
        {
            Real invLength = ((Real)1) / length;
            mTuple[0] *= invLength;
            mTuple[1] *= invLength;
            mTuple[2] *= invLength;
            mTuple[3] *= invLength;
        }
        else
        {
            length = (Real)0;
            mTuple[0] = (Real)0;
            mTuple[1] = (Real)0;
            mTuple[2] = (Real)0;
            mTuple[3] = (Real)0;
        }

        return length;
    }

    THQuaternion Inverse() const  // apply to non-zero quaternion
    {
        THQuaternion inverse;

        Real norm = SquaredLength();
        if (norm > (Real)0)
        {
            Real invNorm = ((Real)1) / norm;
            inverse.mTuple[0] = mTuple[0] * invNorm;
            inverse.mTuple[1] = -mTuple[1] * invNorm;
            inverse.mTuple[2] = -mTuple[2] * invNorm;
            inverse.mTuple[3] = -mTuple[3] * invNorm;
        }
        else
        {
            // Return an invalid result to flag the error.
            for (int i = 0; i < 4; ++i)
            {
                inverse.mTuple[i] = (Real)0;
            }
        }

        return inverse;
    }

    THQuaternion Conjugate() const  // negate x, y, and z terms
    {
        return THQuaternion(mTuple[0], -mTuple[1], -mTuple[2], -mTuple[3]);
    }

    THQuaternion Exp() const  // apply to quaternion with w = 0
    {
        // If q = A*(x*i+y*j+z*k) where (x,y,z) is unit length, then
        // exp(q) = cos(A)+sin(A)*(x*i+y*j+z*k).  If sin(A) is near zero,
        // use exp(q) = cos(A)+A*(x*i+y*j+z*k) since A/sin(A) has limit 1.

        THQuaternion result;

        Real angle = sqrt(mTuple[1] * mTuple[1] +
            mTuple[2] * mTuple[2] + mTuple[3] * mTuple[3]);

        Real sn = sin(angle);
        result.mTuple[0] = cos(angle);

        int i;

        if (fabs(sn) > (Real)0)
        {
            Real coeff = sn / angle;
            for (i = 1; i < 4; ++i)
            {
                result.mTuple[i] = coeff * mTuple[i];
            }
        }
        else
        {
            for (i = 1; i < 4; ++i)
            {
                result.mTuple[i] = mTuple[i];
            }
        }

        return result;
    }

    THQuaternion Log() const  // apply to unit-length quaternion
    {
        // If q = cos(A)+sin(A)*(x*i+y*j+z*k) where (x,y,z) is unit length, then
        // log(q) = A*(x*i+y*j+z*k).  If sin(A) is near zero, use log(q) =
        // sin(A)*(x*i+y*j+z*k) since sin(A)/A has limit 1.

        THQuaternion result;
        result.mTuple[0] = (Real)0;

        int i;

        if (fabs(mTuple[0]) < (Real)1)
        {
            Real angle = acos(mTuple[0]);
            Real sn = sin(angle);
            if (fabs(sn) > (Real)0)
            {
                Real coeff = angle / sn;
                for (i = 1; i < 4; ++i)
                {
                    result.mTuple[i] = coeff * mTuple[i];
                }
                return result;
            }
        }

        for (i = 1; i < 4; ++i)
        {
            result.mTuple[i] = mTuple[i];
        }
        return result;
    }

    // Rotation of a vector by a quaternion.
    TAVector<Real> Rotate (const TAVector<Real>& vec) const
    {
        // Given a vector u = (x0,y0,z0) and a unit length quaternion
        // q = <w,x,y,z>, the vector v = (x1,y1,z1) which represents the
        // rotation of u by q is v = q*u*q^{-1} where * indicates quaternion
        // multiplication and where u is treated as the quaternion <0,x0,y0,z0>.
        // Note that q^{-1} = <w,-x,-y,-z>, so no real work is required to
        // invert q.  Now
        //
        //   q*u*q^{-1} = q*<0,x0,y0,z0>*q^{-1}
        //     = q*(x0*i+y0*j+z0*k)*q^{-1}
        //     = x0*(q*i*q^{-1})+y0*(q*j*q^{-1})+z0*(q*k*q^{-1})
        //
        // As 3-vectors, q*i*q^{-1}, q*j*q^{-1}, and 2*k*q^{-1} are the columns
        // of the rotation matrix computed in HQuaternion::ToRotationMatrix.
        // The vector v is obtained as the product of that rotation matrix with
        // vector u.  As such, the quaternion representation of a rotation
        // matrix requires less space than the matrix and more time to compute
        // the rotated vector.  Typical space-time tradeoff...

        THMatrix<Real> rot;
        ToRotationMatrix(rot);
        return rot * vec;
    }

    // Spherical linear interpolation.
    THQuaternion& Slerp(Real t, const THQuaternion& p, const THQuaternion& q)
    {
        Real cs = p.Dot(q);
        Real angle = acos(cs);

        if (fabs(angle) > (Real)0)
        {
            Real sn = sin(angle);
            Real invSn = ((Real)1) / sn;
            Real tAngle = t * angle;
            Real coeff0 = sin(angle - tAngle) * invSn;
            Real coeff1 = sin(tAngle) * invSn;

            mTuple[0] = coeff0 * p.mTuple[0] + coeff1 * q.mTuple[0];
            mTuple[1] = coeff0 * p.mTuple[1] + coeff1 * q.mTuple[1];
            mTuple[2] = coeff0 * p.mTuple[2] + coeff1 * q.mTuple[2];
            mTuple[3] = coeff0 * p.mTuple[3] + coeff1 * q.mTuple[3];
        }
        else
        {
            mTuple[0] = p.mTuple[0];
            mTuple[1] = p.mTuple[1];
            mTuple[2] = p.mTuple[2];
            mTuple[3] = p.mTuple[3];
        }

        return *this;
    }

    // Intermediate terms for spherical quadratic interpolation.
    THQuaternion& Intermediate(const THQuaternion& q0, const THQuaternion& q1,
        const THQuaternion& q2)
    {
        // assert:  Q0, Q1, Q2 all unit-length
        THQuaternion q1Inv = q1.Conjugate();
        THQuaternion p0 = q1Inv * q0;
        THQuaternion p2 = q1Inv * q2;
        THQuaternion arg = ((Real)-0.25) * (p0.Log() + p2.Log());
        THQuaternion a = q1 * arg.Exp();
        *this = a;
        return *this;
    }

    // Spherical quadratic interpolation.
    THQuaternion& Squad(Real t, const THQuaternion& q0, const THQuaternion& a0,
        const THQuaternion& a1, const THQuaternion& q1)
    {
        Real slerpT = ((Real)2) * t * ((Real)1 - t);
        THQuaternion slerpP = Slerp(t, q0, q1);
        THQuaternion slerpQ = Slerp(t, a0, a1);
        return Slerp(slerpT, slerpP, slerpQ);
    }

    // Special quaternions.
    static const THQuaternion ZERO()
    {
        return THQuaternion((Real)0, (Real)0, (Real)0, (Real)0);
    }

    static const THQuaternion IDENTITY()
    {
        return THQuaternion((Real)1, (Real)0, (Real)0, (Real)0);
    }

private:
    // Order of storage is (w,x,y,z).
    Real mTuple[4];
};

template <typename Real>
THQuaternion<Real> operator* (Real scalar, const THQuaternion<Real>& q)
{
    return q * scalar;
}

}

#endif
