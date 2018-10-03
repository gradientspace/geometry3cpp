// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <Physics/GteFluid3AdjustVelocity.h>
#include <Physics/GteFluid3ComputeDivergence.h>
#include <Physics/GteFluid3EnforceStateBoundary.h>
#include <Physics/GteFluid3InitializeSource.h>
#include <Physics/GteFluid3InitializeState.h>
#include <Physics/GteFluid3SolvePoisson.h>
#include <Physics/GteFluid3UpdateState.h>

namespace gte
{

class GraphicsEngine;
class ProgramFactory;

class GTE_IMPEXP Fluid3
{
public:
    // Construction.  The (x,y,z) grid covers [0,1]^3.
    Fluid3(std::shared_ptr<GraphicsEngine> const& engine,
        std::shared_ptr<ProgramFactory> const& factory,
        int xSize, int ySize, int zSize, float dt);

    void Initialize();
    void DoSimulationStep();
    inline std::shared_ptr<Texture3> const& GetState() const;

private:
    // Constructor inputs.
    std::shared_ptr<GraphicsEngine> mEngine;
    int mXSize, mYSize, mZSize;
    float mDt;

    // Current simulation time.
    float mTime;

    std::shared_ptr<ConstantBuffer> mParameters;
    std::shared_ptr<Fluid3InitializeSource> mInitializeSource;
    std::shared_ptr<Fluid3InitializeState> mInitializeState;
    std::shared_ptr<Fluid3EnforceStateBoundary> mEnforceStateBoundary;
    std::shared_ptr<Fluid3UpdateState> mUpdateState;
    std::shared_ptr<Fluid3ComputeDivergence> mComputeDivergence;
    std::shared_ptr<Fluid3SolvePoisson> mSolvePoisson;
    std::shared_ptr<Fluid3AdjustVelocity> mAdjustVelocity;

    std::shared_ptr<Texture3> mSourceTexture;
    std::shared_ptr<Texture3> mStateTm1Texture;
    std::shared_ptr<Texture3> mStateTTexture;
    std::shared_ptr<Texture3> mStateTp1Texture;
    std::shared_ptr<Texture3> mDivergenceTexture;
    std::shared_ptr<Texture3> mPoissonTexture;
};

inline std::shared_ptr<Texture3> const& Fluid3::GetState() const
{
    return mStateTTexture;
}

}
