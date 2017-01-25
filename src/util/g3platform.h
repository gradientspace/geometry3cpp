#pragma once

#ifdef __APPLE__
    #define memcpy_s(dest, dest_size, src, count) memcpy((dest), (src), (count))
#endif

// if this is defined, then we do not have TBB dependency
//#define G3_DISABLE_TBB