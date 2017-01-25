#include "Wm5MathematicsPCH.h"
#include "Wm5AxisAlignedBox3.h"

namespace Wm5
{
	template<> const AxisAlignedBox3<float> AxisAlignedBox3<float>::EMPTY (
		std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(),
		std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(),
		std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());

	template<> const AxisAlignedBox3<double> AxisAlignedBox3<double>::EMPTY (
		std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(),
		std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(),
		std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());

}
