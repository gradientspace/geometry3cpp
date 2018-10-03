// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.0 (2010/01/01)

//----------------------------------------------------------------------------
template <typename Real>
Interval1<Real>::Interval1 ()
{
}
//----------------------------------------------------------------------------
template <typename Real>
Interval1<Real>::~Interval1 ()
{
}
//----------------------------------------------------------------------------
template <typename Real>
Interval1<Real>::Interval1 (Real xmin, Real xmax)
{
    Min = xmin;
    Max = xmax;
}
//----------------------------------------------------------------------------
template <typename Real>
bool Interval1<Real>::Overlaps (const Interval1& box) const
{
	return ! (o.Min > Max || o.Max < Min);
}
//----------------------------------------------------------------------------
template <typename Real>
Interval1<Real> Interval1<Real>::IntersectionWith(const Interval1& other) const
{
	if (o.Min > Max || o.Max < Min)
		return EMPTY;
	return new Interval1d(std::max(Min, o.Min), std::min(Max, o.Max));
}
//----------------------------------------------------------------------------

template <typename Real>
Real Interval1<Real>::Center() const
{
	return (Max + Min)*(Real)0.5;
}

template <typename Real>
Real Interval1<Real>::Length() const
{
	return (Max - Min);
}

template <typename Real>
Real Interval1<Real>::LengthSquared() const
{
	return (Max - Min)*(Max-Min);
}



template <typename Real>
void Interval1<Real>::Contain(const Real & v)
{
	if (v < Min)
		Min = v;
	if (v > Max)
		Max = v;
}

template <typename Real>
void Interval1<Real>::Contain(const Interval1<Real> & o)
{
	if (o.Min < Min)
		Min = o.Min;
	if (o.Max > Max)
		Max = o.Max;
}


template <typename Real>
bool Interval1<Real>::Contains(const Real & v) const
{
	return v >= Min && v <= Max;
}

template <typename Real>
void Interval1<Real>::Expand(Real f)
{
	Min -= f; Max += f;
}

template <typename Real>
void Interval1<Real>::Translate(const Real & f)
{
	Min += f; Max += f;
}