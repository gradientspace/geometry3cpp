#pragma once

#include <vector>
#include <g3Debug.h>

namespace g3
{

template<typename T>
void assert_same( const std::vector<T>& v0, const std::vector<T>& v1 )
{
	gDevAssert( v0.size() == v1.size() );
	for ( unsigned int k = 0; k < v0.size(); ++k )
		gDevAssert( v0[k] == v1[k] );
}



}


