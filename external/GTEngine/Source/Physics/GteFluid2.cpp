// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#include <GTEnginePCH.h>
#include <Physics/GteFluid2.h>
using namespace gte;

Fluid2::Fluid2(std::shared_ptr<GraphicsEngine> const& engine,
    std::shared_ptr<ProgramFactory> const& factory,
    int xSize, int ySize, float dt, float densityViscosity, float velocityViscosity)
    :
    mEngine(engine),
    mXSize(xSize),
    mYSize(ySize),
    mDt(dt),
    mTime(0.0f)
{
    // Create the shared parameters for many of the simulation shaders.
    float dx = 1.0f/static_cast<float>(mXSize);
    float dy = 1.0f/static_cast<float>(mYSize);
    float dtDivDxDx = (dt/dx)/dx;
    float dtDivDyDy = (dt/dy)/dy;
    float ratio = dx/dy;
    float ratioSqr = ratio*ratio;
    float factor = 0.5f/(1.0f + ratioSqr);
    float epsilonX = factor;
    float epsilonY = ratioSqr*factor;
    float epsilon0 = dx*dx*factor;
    float denVX = densityViscosity*dtDivDxDx;
    float denVY = densityViscosity*dtDivDyDy;
    float velVX = velocityViscosity*dtDivDxDx;
    float velVY = velocityViscosity*dtDivDyDy;

    mParameters = std::make_shared<ConstantBuffer>(sizeof(Fluid2Parameters), false);
    Fluid2Parameters& p = *mParameters->Get<Fluid2Parameters>();
    p.spaceDelta = { dx, dy, 0.0f, 0.0f };
    p.halfDivDelta = { 0.5f / dx, 0.5f / dy, 0.0f, 0.0f };
    p.timeDelta = { dt / dx, dt / dy, 0.0f, dt };
    p.viscosityX = { velVX, velVX, 0.0f, denVX };
    p.viscosityY = { velVY, velVY, 0.0f, denVY };
    p.epsilon = { epsilonX, epsilonY, 0.0f, epsilon0 };

    // Create the compute shaders and textures for the simulation.
    mInitializeSource = std::make_shared<Fluid2InitializeSource>(factory,
        mXSize, mYSize, 16, 16, mParameters);
    mSourceTexture = mInitializeSource->GetSource();

    mInitializeState = std::make_shared<Fluid2InitializeState>(factory,
        mXSize, mYSize, 16, 16);
    mStateTm1Texture = mInitializeState->GetStateTm1();
    mStateTTexture = mInitializeState->GetStateT();

    mEnforceStateBoundary = std::make_shared<Fluid2EnforceStateBoundary>(
        factory, mXSize, mYSize, 16, 16);

    mUpdateState = std::make_shared<Fluid2UpdateState>(factory, mXSize,
        mYSize, 16, 16, mParameters);
    mStateTp1Texture = mUpdateState->GetUpdateState();

    mComputeDivergence = std::make_shared<Fluid2ComputeDivergence>(factory,
        mXSize, mYSize, 16, 16, mParameters);
    mDivergenceTexture = mComputeDivergence->GetDivergence();

    mSolvePoisson = std::make_shared<Fluid2SolvePoisson>(factory, mXSize,
        mYSize, 16, 16, mParameters, 32);
    mPoissonTexture = mSolvePoisson->GetPoisson();

    mAdjustVelocity = std::make_shared<Fluid2AdjustVelocity>(factory, mXSize,
        mYSize, 16, 16, mParameters);
}

void Fluid2::Initialize()
{
    mInitializeSource->Execute(mEngine);
    mInitializeState->Execute(mEngine);
    mEnforceStateBoundary->Execute(mEngine, mStateTm1Texture);
    mEnforceStateBoundary->Execute(mEngine, mStateTTexture);
}

void Fluid2::DoSimulationStep()
{
    mUpdateState->Execute(mEngine, mSourceTexture, mStateTm1Texture, mStateTTexture);
    mEnforceStateBoundary->Execute(mEngine, mStateTp1Texture);
    mComputeDivergence->Execute(mEngine, mStateTp1Texture);
    mSolvePoisson->Execute(mEngine, mDivergenceTexture);
    mAdjustVelocity->Execute(mEngine, mStateTp1Texture, mPoissonTexture, mStateTm1Texture);
    mEnforceStateBoundary->Execute(mEngine, mStateTm1Texture);
    std::swap(mStateTm1Texture, mStateTTexture);

    mTime += mDt;
}

