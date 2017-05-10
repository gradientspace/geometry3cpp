// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <Physics/GteFluid3Parameters.h>
#include <Graphics/GteConstantBuffer.h>
#include <Graphics/GteProgramFactory.h>
#include <Graphics/GteTexture3.h>

namespace gte
{

class GraphicsEngine;

class GTE_IMPEXP Fluid3AdjustVelocity
{
public:
    // Construction.  Adjust the velocities using the solution to the
    // Poisson equation.
    Fluid3AdjustVelocity(std::shared_ptr<ProgramFactory> const& factory,
        int xSize, int ySize, int zSize, int numXThreads, int numYThreads, int numZThreads,
        std::shared_ptr<ConstantBuffer> const& parameters);

    // Update the state for the fluid simulation.
    void Execute(std::shared_ptr<GraphicsEngine> const& engine,
        std::shared_ptr<Texture3> const& inState,
        std::shared_ptr<Texture3> const& poisson,
        std::shared_ptr<Texture3> const& outState);

private:
    int mNumXGroups, mNumYGroups, mNumZGroups;
    std::shared_ptr<ComputeProgram> mAdjustVelocity;

    // Shader source code as strings.
    static std::string const msGLSLSource;
    static std::string const msHLSLSource;
    static std::string const* msSource[ProgramFactory::PF_NUM_API];
};

}
