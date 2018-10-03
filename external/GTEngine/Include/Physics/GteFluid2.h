// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <Physics/GteFluid2AdjustVelocity.h>
#include <Physics/GteFluid2ComputeDivergence.h>
#include <Physics/GteFluid2EnforceStateBoundary.h>
#include <Physics/GteFluid2InitializeSource.h>
#include <Physics/GteFluid2InitializeState.h>
#include <Physics/GteFluid2SolvePoisson.h>
#include <Physics/GteFluid2UpdateState.h>

namespace gte
{

class GraphicsEngine;
class ProgramFactory;

class GTE_IMPEXP Fluid2
{
public:
    // Construction.  The (x,y) grid covers [0,1]^2.
    Fluid2(std::shared_ptr<GraphicsEngine> const& engine,
        std::shared_ptr<ProgramFactory> const& factory,
        int xSize, int ySize, float dt, float densityViscosity, float velocityViscosity);

    void Initialize();
    void DoSimulationStep();
    inline std::shared_ptr<Texture2> const& GetState() const;

private:
    // Constructor inputs.
    std::shared_ptr<GraphicsEngine> mEngine;
    int mXSize, mYSize;
    float mDt;

    // Current simulation time.
    float mTime;

    std::shared_ptr<ConstantBuffer> mParameters;
    std::shared_ptr<Fluid2InitializeSource> mInitializeSource;
    std::shared_ptr<Fluid2InitializeState> mInitializeState;
    std::shared_ptr<Fluid2EnforceStateBoundary> mEnforceStateBoundary;
    std::shared_ptr<Fluid2UpdateState> mUpdateState;
    std::shared_ptr<Fluid2ComputeDivergence> mComputeDivergence;
    std::shared_ptr<Fluid2SolvePoisson> mSolvePoisson;
    std::shared_ptr<Fluid2AdjustVelocity> mAdjustVelocity;

    std::shared_ptr<Texture2> mSourceTexture;
    std::shared_ptr<Texture2> mStateTm1Texture;
    std::shared_ptr<Texture2> mStateTTexture;
    std::shared_ptr<Texture2> mStateTp1Texture;
    std::shared_ptr<Texture2> mDivergenceTexture;
    std::shared_ptr<Texture2> mPoissonTexture;
};

inline std::shared_ptr<Texture2> const& Fluid2::GetState() const
{
    return mStateTTexture;
}

}
