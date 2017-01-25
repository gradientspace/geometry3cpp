// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <Mathematics/GteBasisFunction.h>
#include <Mathematics/GteParametricCurve.h>

namespace gte
{

template <int N, typename Real>
class NURBSCurve : public ParametricCurve<N, Real>
{
public:
    // Construction and destruction.  This object makes copies of the input
    // arrays.  The domain is t in [t[d],t[n]], where t[d] and t[n] are knots
    // with d the degree and n the number of control points.  To validate
    // construction, create an object as shown:
    //     NURBSCurve<N, Real> curve(parameters);
    //     if (!curve) { <constructor failed, handle accordingly>; }
    virtual ~NURBSCurve();
    NURBSCurve(BasisFunctionInput<Real> const& input,
        Vector<N, Real> const* controls, Real const* weights);

    // Member access.
    inline int GetNumControls() const;
    inline Vector<N, Real> const* GetControls() const;
    inline Real const* GetWeights() const;
    inline BasisFunction<Real> const& GetBasisFunction() const;
    void SetControl(int i, Vector<N, Real> const& control);
    Vector<N, Real> const& GetControl(int i) const;
    void SetWeight(int i, Real weight);
    Real const& GetWeight(int i) const;

    // Evaluation of the curve.  The function supports derivative calculation
    // through order 3; that is, maxOrder <= 3 is required.  If you want
    // only the position, pass in maxOrder of 0.  If you want the position and
    // first derivative, pass in maxOrder of 1, and so on.  The output
    // 'values' are ordered as: position, first derivative, second derivative,
    // third derivative.
    virtual void Evaluate(Real t, unsigned int maxOrder,
        Vector<N, Real> values[4]) const;

private:
    // Support for Evaluate(...).
    void Compute(unsigned int order, int imin, int imax, Vector<N, Real>& X,
        Real& w) const;

    BasisFunction<Real> mBasisFunction;
    std::vector<Vector<N, Real>> mControls;
    std::vector<Real> mWeights;
};


template <int N, typename Real>
NURBSCurve<N, Real>::~NURBSCurve()
{
}

template <int N, typename Real>
NURBSCurve<N, Real>::NURBSCurve(BasisFunctionInput<Real> const& input,
    Vector<N, Real> const* controls, Real const* weights)
    :
    ParametricCurve<N, Real>((Real)0, (Real)1),
    mBasisFunction(input)
{
    if (!controls || !weights)
    {
        LogError("Invalid controls or weights pointer.");
        return;
    }

    if (!mBasisFunction)
    {
        // Errors were already generated during construction of the
        // basis function.
        return;
    }

    // The mBasisFunction stores the domain but so does ParametricCurve.
    this->mTime.front() = mBasisFunction.GetMinDomain();
    this->mTime.back() = mBasisFunction.GetMaxDomain();

    // The replication of control points for periodic splines is avoided
    // by wrapping the i-loop index in Evaluate.
    mControls.resize(input.numControls);
    mWeights.resize(input.numControls);
    std::copy(controls, controls + input.numControls, mControls.begin());
    std::copy(weights, weights + input.numControls, mWeights.begin());
    this->mConstructed = true;
}

template <int N, typename Real>
int NURBSCurve<N, Real>::GetNumControls() const
{
    return static_cast<int>(mControls.size());
}

template <int N, typename Real>
Vector<N, Real> const* NURBSCurve<N, Real>::GetControls() const
{
    return &mControls[0];
}

template <int N, typename Real>
Real const* NURBSCurve<N, Real>::GetWeights() const
{
    return &mWeights[0];
}

template <int N, typename Real>
BasisFunction<Real> const& NURBSCurve<N, Real>::GetBasisFunction() const
{
    return mBasisFunction;
}

template <int N, typename Real>
void NURBSCurve<N, Real>::SetControl(int i, Vector<N, Real> const& control)
{
    if (0 <= i && i < GetNumControls())
    {
        mControls[i] = control;
    }
}

template <int N, typename Real>
Vector<N, Real> const& NURBSCurve<N, Real>::GetControl(int i) const
{
    if (0 <= i && i < GetNumControls())
    {
        return mControls[i];
    }
    else
    {
        // Invalid index, return something.
        return mControls[0];
    }
}

template <int N, typename Real>
void NURBSCurve<N, Real>::SetWeight(int i, Real weight)
{
    if (0 <= i && i < GetNumControls())
    {
        mWeights[i] = weight;
    }
}

template <int N, typename Real>
Real const& NURBSCurve<N, Real>::GetWeight(int i) const
{
    if (0 <= i && i < GetNumControls())
    {
        return mWeights[i];
    }
    else
    {
        // Invalid index, return something.
        return mWeights[0];
    }
}

template <int N, typename Real>
void NURBSCurve<N, Real>::Evaluate(Real t, unsigned int maxOrder,
    Vector<N, Real> values[4]) const
{
    if (!this->mConstructed)
    {
        // Errors were already generated during construction.
        for (unsigned int order = 0; order < 4; ++order)
        {
            values[order].MakeZero();
        }
        return;
    }

    int imin, imax;
    mBasisFunction.Evaluate(t, maxOrder, imin, imax);

    // Compute position.
    Vector<N, Real> X;
    Real w;
    Compute(0, imin, imax, X, w);
    Real invW = ((Real)1) / w;
    values[0] = invW * X;

    if (maxOrder >= 1)
    {
        // Compute first derivative.
        Vector<N, Real> XDer1;
        Real wDer1;
        Compute(1, imin, imax, XDer1, wDer1);
        values[1] = invW * (XDer1 - wDer1 * values[0]);

        if (maxOrder >= 2)
        {
            // Compute second derivative.
            Vector<N, Real> XDer2;
            Real wDer2;
            Compute(2, imin, imax, XDer2, wDer2);
            values[2] = invW * (XDer2 - ((Real)2) * wDer1 * values[1] -
                wDer2 * values[0]);

            if (maxOrder == 3)
            {
                // Compute third derivative.
                Vector<N, Real> XDer3;
                Real wDer3;
                Compute(3, imin, imax, XDer3, wDer3);
                values[3] = invW * (XDer3 - ((Real)3) * wDer1 * values[2] -
                    ((Real)3) * wDer2 * values[1] - wDer3 * values[0]);
            }
            else
            {
                values[3].MakeZero();
            }
        }
    }
}

template <int N, typename Real>
void NURBSCurve<N, Real>::Compute(unsigned int order, int imin, int imax,
    Vector<N, Real>& X, Real& w) const
{
    // The j-index introduces a tiny amount of overhead in order to handle
    // both aperiodic and periodic splines.  For aperiodic splines, j = i
    // always.

    int numControls = GetNumControls();
    X.MakeZero();
    w = (Real)0;
    for (int i = imin; i <= imax; ++i)
    {
        int j = (i >= numControls ? i - numControls : i);
        Real tmp = mBasisFunction.GetValue(order, i) * mWeights[j];
        X += tmp * mControls[j];
        w += tmp;
    }
}


}
