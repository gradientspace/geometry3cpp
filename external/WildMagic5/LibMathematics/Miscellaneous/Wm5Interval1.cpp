#include "Wm5MathematicsPCH.h"
#include "Wm5Interval1.h"

namespace Wm5
{
	template<> const Interval1<float> Interval1<float>::EMPTY(
		std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());

	template<> const Interval1<double> Interval1<double>::EMPTY(
		std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());
}
