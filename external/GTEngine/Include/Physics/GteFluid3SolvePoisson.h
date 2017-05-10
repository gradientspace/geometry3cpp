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

class GTE_IMPEXP Fluid3SolvePoisson
{
public:
    // Construction.  Solve the Poisson equation where numIterations is the
    // number of Gauss-Seidel steps to use in Execute.
    Fluid3SolvePoisson(std::shared_ptr<ProgramFactory> const& factory,
        int xSize, int ySize, int zSize, int numXThreads, int numYThreads, int numZThreads,
        std::shared_ptr<ConstantBuffer> const& parameters, int numIterations);

    // Member access.  The texels are (velocity.xyz, density).
    inline std::shared_ptr<gte::Texture3> const& GetPoisson() const;

    // Update the state for the fluid simulation.
    void Execute(std::shared_ptr<GraphicsEngine> const& engine,
        std::shared_ptr<Texture3> const& divergence);

private:
    int mNumXGroups, mNumYGroups, mNumZGroups;
    std::shared_ptr<ComputeProgram> mZeroPoisson;
    std::shared_ptr<ComputeProgram> mSolvePoisson;
    std::shared_ptr<ComputeProgram> mWriteXFace;
    std::shared_ptr<ComputeProgram> mWriteYFace;
    std::shared_ptr<ComputeProgram> mWriteZFace;
    std::shared_ptr<Texture3> mPoisson0;
    std::shared_ptr<Texture3> mPoisson1;
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

inline std::shared_ptr<Texture3> const& Fluid3SolvePoisson::GetPoisson() const
{
    return mPoisson0;
}

}

