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

class GTE_IMPEXP Fluid2ComputeDivergence
{
public:
    // Construction.  Compute the divergence of the velocity vector field:
    // Divergence(V(x,y)) = dV/dx + dV/dy.  The derivatives are estimated
    // using centered finite differences,
    //   dV/dx = (V(x+dx,y) - V(x-dx,y))/(2*dx)
    //   dV/dy = (V(x,y+dy) - V(x,y-dy))/(2*dy)
    Fluid2ComputeDivergence(std::shared_ptr<ProgramFactory> const& factory,
        int xSize, int ySize, int numXThreads, int numYThreads,
        std::shared_ptr<ConstantBuffer> const& parameters);

    // Member access.
    inline std::shared_ptr<Texture2> const& GetDivergence() const;

    // Compute the divergence of the velocity vector field for the input
    // state.
    void Execute(std::shared_ptr<GraphicsEngine> const& engine,
        std::shared_ptr<Texture2> const& state);

private:
    int mNumXGroups, mNumYGroups;
    std::shared_ptr<ComputeProgram> mComputeDivergence;
    std::shared_ptr<Texture2> mDivergence;

    // Shader source code as strings.
    static std::string const msGLSLSource;
    static std::string const msHLSLSource;
    static std::string const* msSource[ProgramFactory::PF_NUM_API];
};

inline std::shared_ptr<Texture2> const& Fluid2ComputeDivergence::GetDivergence() const
{
    return mDivergence;
}

}
