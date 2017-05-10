// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <GTEngineDEF.h>

// Expose this define if you want GPGPU support in computing any algorithms
// that have a GPU implemetnation.  Alternatively, your application can
// define this in the project settings so that you do not have to edit this
// source file.
//
//#define GTE_COMPUTE_MODEL_ALLOW_GPGPU

// The ComputeModel class allows you to select the type of hardware to use
// in your computational algorithms.
//
// If your computational algorithm requires the GPU, set 'inEngine' to the
// object that will execute the compute shaders.  If your algorithm
// requires CPU multithreading, set the 'inNumThreads' to the desired
// number of threads, presumably 2 or larger.  You can query for the number
// of concurrent hardware threads using std::thread::hardware_concurrency().
// If you want single-threaded computations (on the main thread), set
// inNumThreads to 1.  An example of using this class is
//
//  ComputeModel cmodel(...);
//  if (cmodel.engine)
//  {
//      ComputeUsingGPU(...);
//  }
//  else if (cmodel.numThreads > 1)
//  {
//      ComputeUsingCPUMultipleThreads();
//  }
//  else
//  {
//      ComputeUsingCPUSingleThread();
//  }
// See GenerateMeshUV<Real>::SolveSystem(...) for a concrete example.
//
// Of course, your algorithm can interpret cmodel anyway it likes.  For
// example, you might ignore cmodel.engine if all you care about is
// multithreading on the CPU.

#if defined(GTE_COMPUTE_MODEL_ALLOW_GPGPU)
#include <memory>
#endif

namespace gte
{

#if defined(GTE_COMPUTE_MODEL_ALLOW_GPGPU)
class GraphicsEngine;
class ProgramFactory;
#endif

class GTE_IMPEXP ComputeModel
{
public:
    // Construction and destruction.  You may derive from this class to make
    // additional behavior available to your algorithms.  For example, you
    // might add a callback function to report progress of the algorithm.
    virtual ~ComputeModel()
    {
    }

    ComputeModel()
        :
        numThreads(1)
    {
    }

    ComputeModel(unsigned int inNumThreads)
        :
        numThreads(inNumThreads > 0 ? inNumThreads : 1)
    {
    }

#if defined(GTE_COMPUTE_MODEL_ALLOW_GPGPU)
    ComputeModel(unsigned int inNumThreads,
        std::shared_ptr<GraphicsEngine> const& inEngine,
        std::shared_ptr<ProgramFactory> const& inFactory)
        :
        numThreads(inNumThreads > 0 ? inNumThreads : 1),
        engine(inEngine),
        factory(inFactory)
    {
    }
#endif

    unsigned int numThreads;
#if defined(GTE_COMPUTE_MODEL_ALLOW_GPGPU)
    std::shared_ptr<GraphicsEngine> engine;
    std::shared_ptr<ProgramFactory> factory;
#endif
};

}
