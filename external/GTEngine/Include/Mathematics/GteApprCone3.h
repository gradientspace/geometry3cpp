// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2016
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.14.0 (2018/07/16)

#pragma once

#include <Mathematics/GteFunctions.h>
#include <Mathematics/GteGaussNewtonMinimizer.h>
#include <Mathematics/GteLevenbergMarquardtMinimizer.h>

// The cone vertex is V, the unit-length axis direction is U and the
// cone angle is A in (0,pi/2).  The cone is defined algebraically by
// those points X for which
//     Dot(U,X-V)/Length(X-V) = cos(A)
// This can be written as a quadratic equation
//     (V-X)^T * (cos(A)^2 - U * U^T) * (V-X) = 0
// with the implicit constraint that Dot(U, X-V) > 0 (X is on the
// "positive" cone).  Define W = U/cos(A), so Length(W) > 1 and
//     F(X;V,W) = (V-X)^T * (I - W * W^T) * (V-X) = 0
// The nonlinear least squares fitting of points {X[i]}_{i=0}^{n-1}
// computes V and W to minimize the error function
//     E(V,W) = sum_{i=0}^{n-1} F(X[i];V,W)^2
// I recommend using the Gauss-Newton minimizer when your cone points
// are truly nearly a cone; otherwise, try the Levenberg-Marquardt
// minimizer.

namespace gte
{
    template <typename Real>
    class ApprCone3
    {
    public:
        ApprCone3()
            :
            mNumPoints(0),
            mPoints(nullptr)
        {
            // F[i](V,W) = D^T * (I - W * W^T) * D, D = V - X[i], P = (V,W)
            mFFunction = [this](GVector<Real> const& P, GVector<Real>& F)
            {
                Vector<3, Real> V = { P[0], P[1], P[2] };
                Vector<3, Real> W = { P[3], P[4], P[5] };
                for (int i = 0; i < mNumPoints; ++i)
                {
                    Vector<3, Real> delta = V - mPoints[i];
                    Real deltaDotW = Dot(delta, W);
                    F[i] = Dot(delta, delta) - deltaDotW * deltaDotW;
                }
            };

            // dF[i]/dV = 2 * (D - Dot(W, D) * W)
            // dF[i]/dW = -2 * Dot(W, D) * D
            mJFunction = [this](GVector<Real> const& P, GMatrix<Real>& J)
            {
                Vector<3, Real> V = { P[0], P[1], P[2] };
                Vector<3, Real> W = { P[3], P[4], P[5] };
                for (int row = 0; row < mNumPoints; ++row)
                {
                    Vector<3, Real> delta = V - mPoints[row];
                    Real deltaDotW = Dot(delta, W);
                    Vector<3, Real> temp0 = delta - deltaDotW * W;
                    Vector<3, Real> temp1 = deltaDotW * delta;
                    for (int col = 0; col < 3; ++col)
                    {
                        J(row, col) = (Real)2 * temp0[col];
                        J(row, col + 3) = (Real)-2 * temp1[col];
                    }
                }
            };
        }

        // The parameters coneVertex, coneAxis and coneAngle are in/out
        // variables.  The caller must provide initial guesses for these.
        // The function estimates the cone parameters and returns them.  See
        // GteGaussNewtonMinimizer.h for a description of the least-squares
        // algorithm and the parameters that it requires.
        typename GaussNewtonMinimizer<Real>::Result
        operator()(int numPoints, Vector<3, Real> const* points,
            size_t maxIterations, Real updateLengthTolerance, Real errorDifferenceTolerance,
            Vector<3, Real>& coneVertex, Vector<3, Real>& coneAxis, Real& coneAngle)
        {
            mNumPoints = numPoints;
            mPoints = points;
            GaussNewtonMinimizer<Real> minimizer(6, mNumPoints, mFFunction, mJFunction);

            // The initial guess for the cone vertex.
            GVector<Real> initial(6);
            initial[0] = coneVertex[0];
            initial[1] = coneVertex[1];
            initial[2] = coneVertex[2];

            // The initial guess for the weighted cone axis.
            Normalize(coneAxis);
            coneAxis /= Function<Real>::Cos(coneAngle);
            initial[3] = coneAxis[0];
            initial[4] = coneAxis[1];
            initial[5] = coneAxis[2];

            auto result = minimizer(initial, maxIterations, updateLengthTolerance,
                errorDifferenceTolerance);

            // No test is made for result.converged so that we return some
            // estimates of the cone.  The caller can decide how to respond
            // when result.converged is false.
            for (int i = 0; i < 3; ++i)
            {
                coneVertex[i] = result.minLocation[i];
                coneAxis[i] = result.minLocation[i + 3];
            }
            Real cosConeAngle = std::min((Real)1 / Normalize(coneAxis), (Real)1);
            coneAngle = Function<Real>::ACos(cosConeAngle);

            mNumPoints = 0;
            mPoints = nullptr;
            return result;
        }

        // The parameters coneVertex, coneAxis and coneAngle are in/out
        // variables.  The caller must provide initial guesses for these.
        // The function estimates the cone parameters and returns them.  See
        // GteGaussNewtonMinimizer.h for a description of the least-squares
        // algorithm and the parameters that it requires.
        typename LevenbergMarquardtMinimizer<Real>::Result
        operator()(int numPoints, Vector<3, Real> const* points,
            size_t maxIterations, Real updateLengthTolerance, Real errorDifferenceTolerance,
            Real lambdaFactor, Real lambdaAdjust, size_t maxAdjustments,
            Vector<3, Real>& coneVertex, Vector<3, Real>& coneAxis, Real& coneAngle)
        {
            mNumPoints = numPoints;
            mPoints = points;
            LevenbergMarquardtMinimizer<Real> minimizer(6, mNumPoints, mFFunction, mJFunction);

            // The initial guess for the cone vertex.
            GVector<Real> initial(6);
            initial[0] = coneVertex[0];
            initial[1] = coneVertex[1];
            initial[2] = coneVertex[2];

            // The initial guess for the weighted cone axis.
            Normalize(coneAxis);
            coneAxis /= Function<Real>::Cos(coneAngle);
            initial[3] = coneAxis[0];
            initial[4] = coneAxis[1];
            initial[5] = coneAxis[2];

            auto result = minimizer(initial, maxIterations, updateLengthTolerance,
                errorDifferenceTolerance, lambdaFactor, lambdaAdjust, maxAdjustments);

            // No test is made for result.converged so that we return some
            // estimates of the cone.  The caller can decide how to respond
            // when result.converged is false.
            for (int i = 0; i < 3; ++i)
            {
                coneVertex[i] = result.minLocation[i];
                coneAxis[i] = result.minLocation[i + 3];
            }
            Real cosConeAngle = std::min((Real)1 / Normalize(coneAxis), (Real)1);
            coneAngle = Function<Real>::ACos(cosConeAngle);

            mNumPoints = 0;
            mPoints = nullptr;
            return result;
        }

    private:
        int mNumPoints;
        Vector<3, Real> const* mPoints;
        std::function<void(GVector<Real> const&, GVector<Real>&)> mFFunction;
        std::function<void(GVector<Real> const&, GMatrix<Real>&)> mJFunction;
    };
}
