// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <Graphics/GteConstantBuffer.h>
#include <Graphics/GteProgramFactory.h>
#include <Graphics/GteTexture3.h>

namespace gte
{

class GraphicsEngine;

class GTE_IMPEXP Fluid3InitializeState
{
public:
    // Construction.  The initial velocity is zero and the initial density
    // is randomly generated with values in [0,1].
    Fluid3InitializeState(std::shared_ptr<ProgramFactory> const& factory,
        int xSize, int ySize, int zSize, int numXThreads, int numYThreads, int numZThreads);

    // Member access.  The texels are (velocity.xyz, density).
    inline std::shared_ptr<Texture3> const& GetStateTm1() const;
    inline std::shared_ptr<Texture3> const& GetStateT() const;

    // Compute the initial density and initial velocity for the fluid
    // simulation.
    void Execute(std::shared_ptr<GraphicsEngine> const& engine);

private:
    int mNumXGroups, mNumYGroups, mNumZGroups;
    std::shared_ptr<ComputeProgram> mInitializeState;
    std::shared_ptr<Texture3> mDensity;
    std::shared_ptr<Texture3> mVelocity;
    std::shared_ptr<Texture3> mStateTm1;
    std::shared_ptr<Texture3> mStateT;

    // Shader source code as strings.
    static std::string const msGLSLSource;
    static std::string const msHLSLSource;
    static std::string const* msSource[ProgramFactory::PF_NUM_API];
};

inline std::shared_ptr<Texture3> const& Fluid3InitializeState::GetStateTm1() const
{
    return mStateTm1;
}

inline std::shared_ptr<Texture3> const& Fluid3InitializeState::GetStateT() const
{
    return mStateT;
}

}
