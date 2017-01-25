// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.1 (2010/10/01)

#ifndef WM5IVECTOR4_H
#define WM5IVECTOR4_H

#include "Wm5MathematicsLIB.h"
#include "Wm5IVector.h"

namespace Wm5
{

class WM5_MATHEMATICS_ITEM IVector4 : public IVector<4>
{
public:
    // Construction.
    IVector4 ();
    IVector4 (const IVector4& vec);
    IVector4 (const IVector<4>& vec);
    IVector4 (const int& x, const int& y, const int& z, const int &w);
	IVector4 (const int * pValues);

    // Assignment.
    IVector4& operator= (const IVector4& vec);
    IVector4& operator= (const IVector<4>& vec);

    // returns Dot(this,V)
    int Dot (const IVector4& vec) const;

protected:
    using IVector<4>::mTuple;
};

}

#endif
