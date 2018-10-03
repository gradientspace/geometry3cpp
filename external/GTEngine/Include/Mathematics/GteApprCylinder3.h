// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.3.3 (2017/08/10)

#pragma once

#include <Mathematics/GteCylinder3.h>
#include <Mathematics/GteMatrix3x3.h>
#include <Mathematics/GteSymmetricEigensolver3x3.h>
#include <Mathematics/GteConstants.h>
#include <algorithm>
#include <vector>
#include <thread>

// The algorithm for least-squares fitting of a point set by a cylinder is
// described in
//   http://www.geometrictools.com/Documentation/CylinderFitting.pdf
// This document shows how to compute the cylinder radius r and the cylinder
// axis as a line C+t*W with origin C, unit-length direction W, and any
// real-valued t.  The implementation here adds one addition step.  It
// projects the point set onto the cylinder axis, computes the bounding
// t-interval [tmin,tmax] for the projections, and sets the cylinder center
// to C + 0.5*(tmin+tmax)*W and the cylinder height to tmax-tmin.

namespace gte
{

template <typename Real>
class ApprCylinder3
{
public:
    // Search the hemisphere for a minimum, choose numThetaSamples and
    // numPhiSamples to be positive (and preferably large).  These are used
    // to generate a hemispherical grid of samples to be evaluated to find
    // the cylinder axis-direction W.  If the grid samples is quite large
    // and the number of points to be fitted is large, you most likely will
    // want to run multithreaded.  Set numThreads to 0 to run single-threaded
    // in the main process.  Set numThreads > 0 to run multithreaded.  If
    // either of numThetaSamples or numPhiSamples is zero, the operator() sets
    // the cylinder origin and axis to the zero vectors, the radius and height
    // to zero, and returns std::numeric_limits<Real>::max().
    ApprCylinder3(unsigned int numThreads, unsigned int numThetaSamples,
        unsigned int numPhiSamples);

    // TODO: Add a constructor that will use a numerical minimizer to search
    // the hemisphere efficiently.

    // Choose one of the eigenvectors for the covariance matrix as the
    // cylinder axis direction.  If eigenIndex is 0, the eigenvector
    // associated with the smallest eigenvalue is chosen.  If eigenIndex is 2,
    // the eigenvector associated with the largest eigenvalue is chosen.  If
    // eigenIndex is 1, the eigenvector associated with the median eigenvalue
    // is chosen; keep in mind that this could be the minimum or maximum
    // eigenvalue if the eigenspace has dimension 2 or 3.  If eigenIndex is
    // 3 or larger, the operator() sets the cylinder origin and axis to the
    // zero vectors, the radius and height to zero, and returns
    // std::numeric_limits<Real>::max().
    ApprCylinder3(unsigned int eigenIndex);

    // Choose the cylinder axis.  If cylinderAxis is not the zero vector,
    // the constructor will normalize it.  If cylinderAxis is the zero vector,
    // the operator() sets the cylinder origin and axis to the zero vectors,
    // the radius and height to zero, and returns
    // std::numeric_limits<Real>::max().
    ApprCylinder3(Vector3<Real> const& cylinderAxis);

    // The algorithm must estimate 6 parameters, so the number of points must
    // be at least 6 but preferably larger.  The returned value is the
    // root-mean-square of the least-squares error.  If numPoints is less than
    // 6 or if points is a null pointer, the operator() sets the cylinder origin
    // and axis to the zero vectors, the radius and height to zero, and returns
    // std::numeric_limits<Real>::max().
    Real operator()(unsigned int numPoints, Vector3<Real> const* points,
        Cylinder3<Real>& cylinder);

private:
    enum ConstructorType
    {
        FIT_BY_HEMISPHERE_SEARCH,
        FIT_USING_COVARIANCE_EIGENVECTOR,
        FIT_USING_SPECIFIED_AXIS
    };

    Real ComputeUsingDirection(Vector3<Real>& minPC, Vector3<Real>& minW, Real& minRSqr);
    Real ComputeUsingCovariance(Vector3<Real>& minPC, Vector3<Real>& minW, Real& minRSqr);
    Real ComputeSingleThreaded(Vector3<Real>& minPC, Vector3<Real>& minW, Real& minRSqr);
    Real ComputeMultiThreaded(Vector3<Real>& minPC, Vector3<Real>& minW, Real& minRSqr);
    Real Error(Vector3<Real> const& W, Vector3<Real>& PC, Real& rsqr) const;

    ConstructorType mConstructorType;

    // Parameters for the hemisphere-search constructor.
    unsigned int mNumThreads;
    unsigned int mNumThetaSamples;
    unsigned int mNumPhiSamples;

    // Parameters for the eigenvector-index constructor.
    unsigned int mEigenIndex;

    // Parameters for the specified-axis constructor.
    Vector3<Real> mCylinderAxis;

    // A copy of the input points but translated by their average for
    // numerical robustness.
    std::vector<Vector3<Real>> mX;
    Real mInvNumPoints;
};

template <typename Real>
ApprCylinder3<Real>::ApprCylinder3(unsigned int numThreads, unsigned int numThetaSamples,
    unsigned int numPhiSamples)
    :
    mConstructorType(FIT_BY_HEMISPHERE_SEARCH),
    mNumThreads(numThreads),
    mNumThetaSamples(numThetaSamples),
    mNumPhiSamples(numPhiSamples),
    mEigenIndex(0),
    mInvNumPoints((Real)0)
{
    mCylinderAxis = { 0.0f, 0.0f, 0.0f };
}

template <typename Real>
ApprCylinder3<Real>::ApprCylinder3(unsigned int eigenIndex)
    :
    mConstructorType(FIT_USING_COVARIANCE_EIGENVECTOR),
    mNumThreads(0),
    mNumThetaSamples(0),
    mNumPhiSamples(0),
    mEigenIndex(eigenIndex),
    mInvNumPoints((Real)0)
{
    mCylinderAxis = { 0.0f, 0.0f, 0.0f };
}

template <typename Real>
ApprCylinder3<Real>::ApprCylinder3(Vector3<Real> const& cylinderAxis)
    :
    mConstructorType(FIT_USING_SPECIFIED_AXIS),
    mNumThreads(0),
    mNumThetaSamples(0),
    mNumPhiSamples(0),
    mEigenIndex(0),
    mCylinderAxis(cylinderAxis),
    mInvNumPoints((Real)0)
{
    Normalize(mCylinderAxis, true);
}

template <typename Real>
Real ApprCylinder3<Real>::operator()(unsigned int numPoints, Vector3<Real> const* points,
    Cylinder3<Real>& cylinder)
{
    mX.clear();
    mInvNumPoints = (Real)0;
    cylinder.axis.origin = Vector3<Real>::Zero();
    cylinder.axis.direction = Vector3<Real>::Zero();
    cylinder.radius = (Real)0;
    cylinder.height = (Real)0;

    // Validate the input parameters.
    if (numPoints < 6 || !points)
    {
        return std::numeric_limits<Real>::max();
    }

    mX.resize(numPoints);
    mInvNumPoints = (Real)1 / (Real)numPoints;

    // Copy the points and translate by the average for numerical robustness.
    Vector3<Real> average{ (Real)0, (Real)0, (Real)0 };
    for (unsigned int i = 0; i < numPoints; ++i)
    {
        average += points[i];
    }
    average *= mInvNumPoints;
    for (unsigned int i = 0; i < numPoints; ++i)
    {
        mX[i] = points[i] - average;
    }

    // Fit the points based on which constructor the caller used.  The
    // direction is either estimated or selected directly or indirectly
    // by the caller.  The center and squared radius are estimated.
    Vector3<Real> minPC, minW;
    Real minRSqr, minError;

    if (mConstructorType == FIT_BY_HEMISPHERE_SEARCH)
    {
        // Validate the relevant internal parameters.
        if (mNumThetaSamples == 0 || mNumPhiSamples == 0)
        {
            return std::numeric_limits<Real>::max();
        }

        // Search the hemisphere for the vector that leads to minimum error
        // and use it for the cylinder axis.
        if (mNumThreads == 0)
        {
            // Execute the algorithm in the main process.
            minError = ComputeSingleThreaded(minPC, minW, minRSqr);
        }
        else
        {
            // Execute the algorithm in multiple threads.
            minError = ComputeMultiThreaded(minPC, minW, minRSqr);
        }
    }
    else if (mConstructorType == FIT_USING_COVARIANCE_EIGENVECTOR)
    {
        // Validate the relevant internal parameters.
        if (mEigenIndex >= 3)
        {
            return std::numeric_limits<Real>::max();
        }

        // Use the eigenvector corresponding to the mEigenIndex of the
        // eigenvectors of the covariance matrix as the cylinder axis
        // direction.  The eigenvectors are sorted from smallest
        // eigenvalue (mEigenIndex = 0) to largest eigenvalue
        // (mEigenIndex = 2).
        minError = ComputeUsingCovariance(minPC, minW, minRSqr);
    }
    else  // mConstructorType == FIT_USING_SPECIFIED_AXIS
    {
        // Validate the relevant internal parameters.
        if (mCylinderAxis == Vector3<Real>::Zero())
        {
            return std::numeric_limits<Real>::max();
        }

        minError = ComputeUsingDirection(minPC, minW, minRSqr);
    }

    // Translate back to the original space by the average of the points.
    cylinder.axis.origin = minPC + average;
    cylinder.axis.direction = minW;

    // Compute the cylinder radius.
    cylinder.radius = sqrt(minRSqr);

    // Project the points onto the cylinder axis and choose the cylinder
    // center and cylinder height as described in the comments at the top of
    // this header file.
    Real tmin = (Real)0, tmax = (Real)0;
    for (unsigned int i = 0; i < numPoints; ++i)
    {
        Real t = Dot(cylinder.axis.direction, points[i] - cylinder.axis.origin);
        tmin = std::min(t, tmin);
        tmax = std::max(t, tmax);
    }

    cylinder.axis.origin += ((tmin + tmax) * (Real)0.5) * cylinder.axis.direction;
    cylinder.height = tmax - tmin;
    return minError;
}

template <typename Real>
Real ApprCylinder3<Real>::ComputeUsingDirection(Vector3<Real>& minPC, Vector3<Real>& minW, Real& minRSqr)
{
    minW = mCylinderAxis;
    return Error(minW, minPC, minRSqr);
}

template <typename Real>
Real ApprCylinder3<Real>::ComputeUsingCovariance(Vector3<Real>& minPC, Vector3<Real>& minW, Real& minRSqr)
{
    Matrix3x3<Real> covar = Matrix3x3<Real>::Zero();
    for (auto const& X : mX)
    {
        covar += OuterProduct(X, X);
    }
    covar *= mInvNumPoints;
    std::array<Real, 3> eval;
    std::array<std::array<Real, 3>, 3> evec;
    SymmetricEigensolver3x3<Real>()(
        covar(0, 0), covar(0, 1), covar(0, 2), covar(1, 1), covar(1, 2), covar(2, 2),
        true, +1, eval, evec);
    minW = evec[mEigenIndex];
    return Error(minW, minPC, minRSqr);
}

template <typename Real>
Real ApprCylinder3<Real>::ComputeSingleThreaded(Vector3<Real>& minPC, Vector3<Real>& minW, Real& minRSqr)
{
    Real const iMultiplier = (Real)GTE_C_TWO_PI / (Real)mNumThetaSamples;
    Real const jMultiplier = (Real)GTE_C_HALF_PI / (Real)mNumPhiSamples;

    // Handle the north pole (0,0,1) separately.
    minW = { (Real)0, (Real)0, (Real)1 };
    Real minError = Error(minW, minPC, minRSqr);

    for (unsigned int j = 1; j <= mNumPhiSamples; ++j)
    {
        Real phi = jMultiplier * static_cast<Real>(j);  // in [0,pi/2]
        Real csphi = cos(phi);
        Real snphi = sin(phi);
        for (unsigned int i = 0; i < mNumThetaSamples; ++i)
        {
            Real theta = iMultiplier * static_cast<Real>(i);  // in [0,2*pi)
            Real cstheta = cos(theta);
            Real sntheta = sin(theta);
            Vector3<Real> W{ cstheta * snphi, sntheta * snphi, csphi };
            Vector3<Real> PC;
            Real rsqr;
            Real error = Error(W, PC, rsqr);
            if (error < minError)
            {
                minError = error;
                minRSqr = rsqr;
                minW = W;
                minPC = PC;
            }
        }
    }

    return minError;
}

template <typename Real>
Real ApprCylinder3<Real>::ComputeMultiThreaded(Vector3<Real>& minPC, Vector3<Real>& minW, Real& minRSqr)
{
    Real const iMultiplier = (Real)GTE_C_TWO_PI / (Real)mNumThetaSamples;
    Real const jMultiplier = (Real)GTE_C_HALF_PI / (Real)mNumPhiSamples;

    // Handle the north pole (0,0,1) separately.
    minW = { (Real)0, (Real)0, (Real)1 };
    Real minError = Error(minW, minPC, minRSqr);

    struct Local
    {
        Real error;
        Real rsqr;
        Vector3<Real> W;
        Vector3<Real> PC;
        unsigned int jmin;
        unsigned int jmax;
    };

    std::vector<Local> local(mNumThreads);
    unsigned int numPhiSamplesPerThread = mNumPhiSamples / mNumThreads;
    for (unsigned int t = 0; t < mNumThreads; ++t)
    {
        local[t].error = std::numeric_limits<Real>::max();
        local[t].rsqr = (Real)0;
        local[t].W = Vector3<Real>::Zero();
        local[t].PC = Vector3<Real>::Zero();
        local[t].jmin = numPhiSamplesPerThread * t;
        local[t].jmax = numPhiSamplesPerThread * (t + 1);
    }
    local[mNumThreads - 1].jmax = mNumPhiSamples + 1;

    std::vector<std::thread> process(mNumThreads);
    for (unsigned int t = 0; t < mNumThreads; ++t)
    {
        process[t] = std::thread
        (
            [this, t, iMultiplier, jMultiplier, &local]()
        {
            for (unsigned int j = local[t].jmin; j < local[t].jmax; ++j)
            {
                Real phi = jMultiplier * static_cast<Real>(j);  // in [0,pi/2]
                Real csphi = cos(phi);
                Real snphi = sin(phi);
                for (unsigned int i = 0; i < mNumThetaSamples; ++i)
                {
                    Real theta = iMultiplier * static_cast<Real>(i);  // in [0,2*pi)
                    Real cstheta = cos(theta);
                    Real sntheta = sin(theta);
                    Vector3<Real> W{ cstheta * snphi, sntheta * snphi, csphi };
                    Vector3<Real> PC;
                    Real rsqr;
                    Real error = Error(W, PC, rsqr);
                    if (error < local[t].error)
                    {
                        local[t].error = error;
                        local[t].rsqr = rsqr;
                        local[t].W = W;
                        local[t].PC = PC;
                    }
                }
            }
        }
        );
    }

    for (unsigned int t = 0; t < mNumThreads; ++t)
    {
        process[t].join();

        if (local[t].error < minError)
        {
            minError = local[t].error;
            minRSqr = local[t].rsqr;
            minW = local[t].W;
            minPC = local[t].PC;
        }
    }

    return minError;
}

template <typename Real>
Real ApprCylinder3<Real>::Error(Vector3<Real> const& W, Vector3<Real>& PC, Real& rsqr) const
{
    Matrix3x3<Real> P = Matrix3x3<Real>::Identity() - OuterProduct(W, W);
    Matrix3x3<Real> S
    {
        (Real)0, -W[2], W[1],
        W[2], (Real)0, -W[0],
        -W[1], W[0], (Real)0
    };

    std::vector<Vector3<Real>> Y(mX.size());
    std::vector<Real> sqrLength(mX.size());

    Matrix3x3<Real> A = Matrix3x3<Real>::Zero();
    Vector3<Real>B = Vector3<Real>::Zero();
    Real qform = (Real)0;
    for (size_t i = 0; i < mX.size(); ++i)
    {
        Vector3<Real> projection = P * mX[i];
        Y[i] = projection;
        sqrLength[i] = Dot(projection, projection);
        A += OuterProduct(projection, projection);
        B += sqrLength[i] * projection;
        qform += sqrLength[i];
    }
    A *= mInvNumPoints;
    B *= mInvNumPoints;
    qform *= mInvNumPoints;

    Matrix3x3<Real> Ahat = -S * A * S;
    PC = (Ahat * B) / Trace<Real>(Ahat * A);

    Real error = (Real)0;
    rsqr = (Real)0;
    for (size_t i = 0; i < mX.size(); ++i)
    {
        Real term = sqrLength[i] - Dot(Y[i], PC) * (Real)2 - qform;
        error += term * term;
        Vector3<Real> diff = PC - Y[i];
        rsqr += Dot(diff, diff);
    }
    error *= mInvNumPoints;
    rsqr *= mInvNumPoints;
    return error;
}

}
