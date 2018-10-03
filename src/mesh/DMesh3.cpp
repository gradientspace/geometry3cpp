#include <geometry3PCH.h>
#include "DMesh3.h"

using namespace g3;

const Vector3d DMesh3::InvalidVertex = Vector3d(max_double, 0, 0);
const Index3i DMesh3::InvalidTriangle = Index3i(InvalidID, InvalidID, InvalidID);
const Index2i DMesh3::InvalidEdge = Index2i(InvalidID, InvalidID);
