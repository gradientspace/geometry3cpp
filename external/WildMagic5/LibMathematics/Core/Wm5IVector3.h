// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.1 (2010/10/01)

#ifndef WM5IVECTOR3_H
#define WM5IVECTOR3_H

#include "Wm5MathematicsLIB.h"
#include "Wm5IVector.h"

namespace Wm5
{

class WM5_MATHEMATICS_ITEM IVector3 : public IVector<3>
{
public:
    // Construction.
    IVector3 ();
    IVector3 (const IVector3& vec);
    IVector3 (const IVector<3>& vec);
    IVector3 (const int& x, const int& y, const int& z);

    // Member access.
    inline int X () const;
    inline int& X ();
    inline int Y () const;
    inline int& Y ();
    inline int Z () const;
    inline int& Z ();

    // Assignment.
    IVector3& operator= (const IVector3& vec);
    IVector3& operator= (const IVector<3>& vec);

    // returns Dot(this,V)
    int Dot (const IVector3& vec) const;

    // returns Cross(this,V)
    IVector3 Cross (const IVector3& vec) const;

    // returns Dot(this,Cross(U,V))
    int TripleScalar (const IVector3& U, const IVector3& V) const;

protected:
    using IVector<3>::mTuple;
};

#include "Wm5IVector3.inl"

}

#endif
