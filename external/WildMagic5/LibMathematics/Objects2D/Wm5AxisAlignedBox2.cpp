#include "Wm5MathematicsPCH.h"
#include "Wm5AxisAlignedBox2.h"

namespace Wm5
{
	template<> const AxisAlignedBox2<float> AxisAlignedBox2<float>::EMPTY(
		std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(),
		std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());

	template<> const AxisAlignedBox2<double> AxisAlignedBox2<double>::EMPTY(
		std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(),
		std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());
}
