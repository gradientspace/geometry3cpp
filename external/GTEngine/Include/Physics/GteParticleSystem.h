// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <Mathematics/GteVector.h>
#include <limits>
#include <vector>

namespace gte
{

template <int N, typename Real>
class ParticleSystem
{
public:
    // Construction and destruction.  If a particle is to be immovable, set
    // its mass to std::numeric_limits<Real>::max().
    virtual ~ParticleSystem();
    ParticleSystem(int numParticles, Real step);

    // Member access.
    inline int GetNumParticles() const;
    void SetMass(int i, Real mass);
    inline void SetPosition(int i, Vector<N, Real> const& position);
    inline void SetVelocity(int i, Vector<N, Real> const& velocity);
    void SetStep(Real step);
    inline Real const& GetMass(int i) const;
    inline Vector<N, Real> const& GetPosition(int i) const;
    inline Vector<N, Real> const& GetVelocity(int i) const;
    inline Real GetStep() const;

    // Update the particle positions based on current time and particle state.
    // The Acceleration(...) function is called in this update for each
    // particle.  This function is virtual so that derived classes can perform
    // pre-update and/or post-update semantics.
    virtual void Update(Real time);

protected:
    // Callback for acceleration (ODE solver uses x" = F/m) applied to
    // particle i.  The positions and velocities are not necessarily
    // mPosition and mVelocity, because the ODE solver evaluates the
    // impulse function at intermediate positions.
    virtual Vector<N, Real> Acceleration(int i, Real time,
        std::vector<Vector<N, Real>> const& position,
        std::vector<Vector<N, Real>> const& velocity) = 0;

    int mNumParticles;
    std::vector<Real> mMass, mInvMass;
    std::vector<Vector<N, Real>> mPosition, mVelocity;
    Real mStep, mHalfStep, mSixthStep;

    // Temporary storage for the Runge-Kutta differential equation solver.
    struct Temporary
    {
        Vector<N, Real> d1, d2, d3, d4;
    };
    std::vector<Vector<N, Real>> mPTmp, mVTmp;
    std::vector<Temporary> mPAllTmp, mVAllTmp;
};


template <int N, typename Real>
ParticleSystem<N, Real>::~ParticleSystem()
{
}

template <int N, typename Real> inline
ParticleSystem<N, Real>::ParticleSystem(int numParticles, Real step)
:
mNumParticles(numParticles),
mMass(numParticles),
mInvMass(numParticles),
mPosition(numParticles),
mVelocity(numParticles),
mStep(step),
mHalfStep(step / (Real)2),
mSixthStep(step / (Real)6),
mPTmp(numParticles),
mVTmp(numParticles),
mPAllTmp(numParticles),
mVAllTmp(numParticles)
{
    std::fill(mMass.begin(), mMass.end(), (Real)0);
    std::fill(mInvMass.begin(), mInvMass.end(), (Real)0);
    std::fill(mPosition.begin(), mPosition.end(), Vector<N, Real>::Zero());
    std::fill(mVelocity.begin(), mVelocity.end(), Vector<N, Real>::Zero());
}

template <int N, typename Real> inline
int ParticleSystem<N, Real>::GetNumParticles() const
{
    return mNumParticles;
}

template <int N, typename Real>
void ParticleSystem<N, Real>::SetMass(int i, Real mass)
{
    if ((Real)0 < mass && mass < std::numeric_limits<Real>::max())
    {
        mMass[i] = mass;
        mInvMass[i] = ((Real)1) / mass;
    }
    else
    {
        mMass[i] = std::numeric_limits<Real>::max();
        mInvMass[i] = (Real)0;
    }
}

template <int N, typename Real> inline
void ParticleSystem<N, Real>::SetPosition(int i,
Vector<N, Real> const& position)
{
    mPosition[i] = position;
}

template <int N, typename Real> inline
void ParticleSystem<N, Real>::SetVelocity(int i,
Vector<N, Real> const& velocity)
{
    mVelocity[i] = velocity;
}

template <int N, typename Real>
void ParticleSystem<N, Real>::SetStep(Real step)
{
    mStep = step;
    mHalfStep = mStep / (Real)2;
    mSixthStep = mStep / (Real)6;
}

template <int N, typename Real> inline
Real const& ParticleSystem<N, Real>::GetMass(int i) const
{
    return mMass[i];
}

template <int N, typename Real> inline
Vector<N, Real> const& ParticleSystem<N, Real>::GetPosition(int i) const
{
    return mPosition[i];
}

template <int N, typename Real> inline
Vector<N, Real> const& ParticleSystem<N, Real>::GetVelocity(int i) const
{
    return mVelocity[i];
}

template <int N, typename Real> inline
Real ParticleSystem<N, Real>::GetStep() const
{
    return mStep;
}

template <int N, typename Real>
void ParticleSystem<N, Real>::Update(Real time)
{
    // Runge-Kutta fourth-order solver.
    Real halfTime = time + mHalfStep;
    Real fullTime = time + mStep;

    // Compute the first step.
    int i;
    for (i = 0; i < mNumParticles; ++i)
    {
        if (mInvMass[i] >(Real)0)
        {
            mPAllTmp[i].d1 = mVelocity[i];
            mVAllTmp[i].d1 = Acceleration(i, time, mPosition, mVelocity);
        }
    }
    for (i = 0; i < mNumParticles; ++i)
    {
        if (mInvMass[i] >(Real)0)
        {
            mPTmp[i] = mPosition[i] + mHalfStep * mPAllTmp[i].d1;
            mVTmp[i] = mVelocity[i] + mHalfStep * mVAllTmp[i].d1;
        }
        else
        {
            mPTmp[i] = mPosition[i];
            mVTmp[i].MakeZero();
        }
    }

    // Compute the second step.
    for (i = 0; i < mNumParticles; ++i)
    {
        if (mInvMass[i] >(Real)0)
        {
            mPAllTmp[i].d2 = mVTmp[i];
            mVAllTmp[i].d2 = Acceleration(i, halfTime, mPTmp, mVTmp);
        }
    }
    for (i = 0; i < mNumParticles; ++i)
    {
        if (mInvMass[i] >(Real)0)
        {
            mPTmp[i] = mPosition[i] + mHalfStep * mPAllTmp[i].d2;
            mVTmp[i] = mVelocity[i] + mHalfStep * mVAllTmp[i].d2;
        }
        else
        {
            mPTmp[i] = mPosition[i];
            mVTmp[i].MakeZero();
        }
    }

    // Compute the third step.
    for (i = 0; i < mNumParticles; ++i)
    {
        if (mInvMass[i] >(Real)0)
        {
            mPAllTmp[i].d3 = mVTmp[i];
            mVAllTmp[i].d3 = Acceleration(i, halfTime, mPTmp, mVTmp);
        }
    }
    for (i = 0; i < mNumParticles; ++i)
    {
        if (mInvMass[i] >(Real)0)
        {
            mPTmp[i] = mPosition[i] + mStep * mPAllTmp[i].d3;
            mVTmp[i] = mVelocity[i] + mStep * mVAllTmp[i].d3;
        }
        else
        {
            mPTmp[i] = mPosition[i];
            mVTmp[i].MakeZero();
        }
    }

    // Compute the fourth step.
    for (i = 0; i < mNumParticles; ++i)
    {
        if (mInvMass[i] >(Real)0)
        {
            mPAllTmp[i].d4 = mVTmp[i];
            mVAllTmp[i].d4 = Acceleration(i, fullTime, mPTmp, mVTmp);
        }
    }
    for (i = 0; i < mNumParticles; ++i)
    {
        if (mInvMass[i] >(Real)0)
        {
            mPosition[i] += mSixthStep * (mPAllTmp[i].d1 +
                ((Real)2) * (mPAllTmp[i].d2 + mPAllTmp[i].d3) +
                mPAllTmp[i].d4);

            mVelocity[i] += mSixthStep * (mVAllTmp[i].d1 +
                ((Real)2) * (mVAllTmp[i].d2 + mVAllTmp[i].d3) +
                mVAllTmp[i].d4);
        }
    }
}


}
