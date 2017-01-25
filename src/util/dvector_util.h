#pragma once

#include <atomic>
#include <g3platform.h>
#include <dvector.h>

#ifndef G3_DISABLE_TBB	
#include <tbb/parallel_for.h>
#endif


namespace g3 {


// Non-TBB versions of these functions must use portable (eg C++11)
//   multi-threading, or do serial computations
#ifdef G3_DISABLE_TBB		


// apply f() to each element of v, parallelized by segment blocks
template<typename Type, typename Func>
void parallel_apply( dvector<Type> & v, const Func & f )
{
	v.apply(f);
}



#else		
// TBB versions of these functions


// apply f() to each element of v, parallelized by segment blocks
template<typename Type, typename Func>
void parallel_apply( dvector<Type> & v, const Func & f )
{
	unsigned int nSegments = (unsigned int)v.m_vSegments.size();
	tbb::parallel_for( tbb::blocked_range<unsigned int>(0, nSegments),  
				  [&] (const tbb::blocked_range<unsigned int>& r) {
		for (unsigned int si = r.begin(); si != r.end(); ++si) {
			dvector_segment<Type> & seg = v.m_vSegments[si];
			for (unsigned int i = 0; i < seg.nCur; ++i) {
				f( seg.pData[i] );
			}
		}
	});
}


#endif


} // end namespace g3