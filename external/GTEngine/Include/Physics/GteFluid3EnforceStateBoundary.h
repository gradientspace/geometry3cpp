// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
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

class GTE_IMPEXP Fluid3EnforceStateBoundary
{
public:
    // Construction.  Clamp the velocity vectors at the image boundary to
    // ensure they do not point outside the domain.  For example, let
    // image(x,y,z) have velocity (vx,vy,vz).  At the boundary pixel
    // (x,y,0), the velocity is clamped to (vx,vy,0).  The density is
    // set to zero on the image boundary.
    Fluid3EnforceStateBoundary(std::shared_ptr<ProgramFactory> const& factory,
        int xSize, int ySize, int zSize, int numXThreads, int numYThreads, int numZThreads);

    // Set the density and velocity values at the image boundary as described
    // in the comments for the constructor.  The state texture has texels
    // (velocity.x, velocity.y, velocity.z, density).
    void Execute(std::shared_ptr<GraphicsEngine> const& engine,
        std::shared_ptr<Texture3> const& state);

private:
    int mNumXGroups, mNumYGroups, mNumZGroups;
    std::shared_ptr<ComputeProgram> mCopyXFace;
    std::shared_ptr<ComputeProgram> mWriteXFace;
    std::shared_ptr<ComputeProgram> mCopyYFace;
    std::shared_ptr<ComputeProgram> mWriteYFace;
    std::shared_ptr<ComputeProgram> mCopyZFace;
    std::shared_ptr<ComputeProgram> mWriteZFace;
    std::shared_ptr<Texture2> mXMin;
    std::shared_ptr<Texture2> mXMax;
    std::shared_ptr<Texture2> mYMin;
    std::shared_ptr<Texture2> mYMax;
    std::shared_ptr<Texture2> mZMin;
    std::shared_ptr<Texture2> mZMax;

    // Shader source code as strings.
    static std::string const msGLSLSource;
    static std::string const msHLSLSource;
    static std::string const* msSource[ProgramFactory::PF_NUM_API];
};

}
