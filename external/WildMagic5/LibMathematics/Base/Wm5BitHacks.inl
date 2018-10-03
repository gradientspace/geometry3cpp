// Geometric Tools, LLC
// Copyright (c) 1998-2014
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.1 (2017/10/12)

//----------------------------------------------------------------------------
inline int ScaledFloatToInt (float value, int power)
{
    int result;

    union { float fValue; int pattern; } u;
    u.fValue = value;

    int shift = 150 - power - ((u.pattern >> 23) & 0xFF);
    if (shift < 24)
    {
        result = ((u.pattern & 0x007FFFFF) | 0x00800000) >> shift;
        if (result == (1 << power))
        {
            --result;
        }
    }
    else
    {
        result = 0;
    }

    return result;
}
//----------------------------------------------------------------------------
