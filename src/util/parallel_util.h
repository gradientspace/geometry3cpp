#pragma once

#include <atomic>
#include <g3platform.h>

#ifdef G3_ENABLE_TBB	
#include <tbb/parallel_for.h>
#endif

namespace g3 {


// Non-TBB versions of these functions must use portable (eg C++11)
//   multi-threading, or do serial computations
#ifndef G3_ENABLE_TBB


// evaluate f[k] = f(k)
template<typename vector_type, typename ValueFunc>
void parallel_fill( vector_type & v, const ValueFunc & f )
{
	size_t nCount = v.size();
	for ( unsigned int k = 0; k < nCount; ++k )
		v[k] = f(k);
}


// evaluate f[k] = f(k), if valid(k)
template<typename vector_type, typename ValueFunc, typename ValidFunc>
void parallel_fill( vector_type & v, const ValueFunc & f, const ValidFunc & valid )
{
	size_t nCount = v.size();
	for (unsigned int k = 0; k < nCount; ++k) {
		if ( valid(k) == true )
			v[k] = f( k );
	}
}



#else
// TBB versions of these functions


// evaluate f[k] = f(k)
template<typename vector_type, typename ValueFunc>
void parallel_fill( vector_type & v, const ValueFunc & f )
{
	unsigned int nCount = (unsigned int)v.size();
	tbb::parallel_for( tbb::blocked_range<unsigned int>(0, nCount),  
					   [&] (const tbb::blocked_range<unsigned int>& r) {
		for (unsigned int k = r.begin(); k != r.end(); ++k)
			v[k] = f(k);
	});
}


// evaluate f[k] = f(k), if valid(k)
template<typename vector_type, typename ValueFunc, typename ValidFunc>
void parallel_fill( vector_type & v, const ValueFunc & f, const ValidFunc & valid )
{
	unsigned int nCount = (unsigned int)v.size();
	tbb::parallel_for( tbb::blocked_range<unsigned int>(0, nCount),  
					   [&] (const tbb::blocked_range<unsigned int>& r) {
		for (unsigned int k = r.begin(); k != r.end(); ++k) {
			if ( valid(k) == true )
				v[k] = f( k );
		}
	});
}


#endif


} // end namespace g3