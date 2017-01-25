#include "Wm5MathematicsPCH.h"
#include "Wm5Box3.h"

namespace Wm5
{
	template<> const Box3<float> Box3<float>::EMPTY (
		Vector3<float>::ZERO, 
		Vector3<float>::UNIT_X, Vector3<float>::UNIT_Y, Vector3<float>::UNIT_Z,
		(float)0, (float)0, (float)0 );

	template<> const Box3<double> Box3<double>::EMPTY (
		Vector3<double>::ZERO, 
		Vector3<double>::UNIT_X, Vector3<double>::UNIT_Y, Vector3<double>::UNIT_Z,
		(double)0, (double)0, (double)0 );

}
