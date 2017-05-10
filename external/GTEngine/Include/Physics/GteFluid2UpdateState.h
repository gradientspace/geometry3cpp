// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <Physics/GteFluid2Parameters.h>
#include <Graphics/GteConstantBuffer.h>
#include <Graphics/GteProgramFactory.h>
#include <Graphics/GteSamplerState.h>
#include <Graphics/GteTexture2.h>

namespace gte
{

class GraphicsEngine;

class GTE_IMPEXP Fluid2UpdateState
{
public:
    // Construction.
    Fluid2UpdateState(std::shared_ptr<ProgramFactory> const& factory,
        int xSize, int ySize, int numXThreads, int numYThreads,
        std::shared_ptr<ConstantBuffer> const& parameters);

    // Member access.  The texels are (velocity.x, velocity.y, 0, density).
    // The third component is unused in the simulation (a 3D simulation will
    // store velocity.z in this component).
    inline std::shared_ptr<Texture2> const& GetUpdateState() const;

    // Update the state for the fluid simulation.
    void Execute(std::shared_ptr<GraphicsEngine> const& engine,
        std::shared_ptr<Texture2> const& source,
        std::shared_ptr<Texture2> const& stateTm1,
        std::shared_ptr<Texture2> const& stateT);

private:
    int mNumXGroups, mNumYGroups;
    std::shared_ptr<ComputeProgram> mComputeUpdateState;
    std::shared_ptr<SamplerState> mAdvectionSampler;
    std::shared_ptr<Texture2> mUpdateState;

    // Shader source code as strings.
    static std::string const msGLSLSource;
    static std::string const msHLSLSource;
    static std::string const* msSource[ProgramFactory::PF_NUM_API];
};

inline std::shared_ptr<Texture2> const& Fluid2UpdateState::GetUpdateState()  const
{
    return mUpdateState;
}

}
