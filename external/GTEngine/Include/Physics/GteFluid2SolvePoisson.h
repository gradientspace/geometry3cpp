// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <Physics/GteFluid2Parameters.h>
#include <Graphics/GteConstantBuffer.h>
#include <Graphics/GteProgramFactory.h>
#include <Graphics/GteTexture2.h>

namespace gte
{

class GraphicsEngine;

class GTE_IMPEXP Fluid2SolvePoisson
{
public:
    // Construction.  Solve the Poisson equation where numIterations is the
    // number of Gauss-Seidel steps to use in Execute.
    Fluid2SolvePoisson(std::shared_ptr<ProgramFactory> const& factory,
        int xSize, int ySize, int numXThreads, int numYThreads,
        std::shared_ptr<ConstantBuffer> const& parameters, int numIterations);

    // Member access.  The texels are (velocity.x, velocity.y, 0, density).
    // The third component is unused in the simulation (a 3D simulation will
    // store velocity.z in this component).
    inline std::shared_ptr<Texture2> const& GetPoisson() const;

    // Update the state for the fluid simulation.
    void Execute(std::shared_ptr<GraphicsEngine> const& engine,
        std::shared_ptr<Texture2> const& divergence);

private:
    int mNumXGroups, mNumYGroups;
    std::shared_ptr<ComputeProgram> mZeroPoisson;
    std::shared_ptr<ComputeProgram> mSolvePoisson;
    std::shared_ptr<ComputeProgram> mWriteXEdge;
    std::shared_ptr<ComputeProgram> mWriteYEdge;
    std::shared_ptr<Texture2> mPoisson0;
    std::shared_ptr<Texture2> mPoisson1;
    int mNumIterations;

    // Shader source code as strings.
    static std::string const msGLSLZeroSource;
    static std::string const msGLSLSolveSource;
    static std::string const msGLSLEnforceSource;
    static std::string const msHLSLZeroSource;
    static std::string const msHLSLSolveSource;
    static std::string const msHLSLEnforceSource;
    static std::string const* msZeroSource[ProgramFactory::PF_NUM_API];
    static std::string const* msSolveSource[ProgramFactory::PF_NUM_API];
    static std::string const* msEnforceSource[ProgramFactory::PF_NUM_API];
};

inline std::shared_ptr<Texture2> const& Fluid2SolvePoisson::GetPoisson() const
{
    return mPoisson0;
}

}
