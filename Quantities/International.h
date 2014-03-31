#pragma once

#include "Quantities.h"
#include "NamedQuantities.h"
#include "Numbers.h"
#include "SI.h"

namespace Principia {
// This namespace contains length units using the foot as defined
// from the metric system by the international yard and pound agreement.
// Only length units used both in the US and UK are listed here.
namespace International {
Quantities::Length Yard = 0.9144 * SI::Metre;
Quantities::Length Foot = Yard / 3;
Quantities::Length Inch = Foot / 12;
}
}
