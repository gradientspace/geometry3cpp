#pragma once

#include "cinterface.h"

#ifdef __cplusplus
extern "C"
{
#endif

g3ExternalC CMeshHandle triangulate_polygon( int nVertices, double * vertexBuffer );


#ifdef __cplusplus
}
#endif
