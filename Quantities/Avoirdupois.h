#pragma once

#include "Quantities.h"
#include "NamedQuantities.h"
#include "Numbers.h"
#include "SI.h"

namespace Principia {
// This namespace contains the avoirdupois mass units using the pound as defined
// from the metric system by the international yard and pound agreement.
// One pounds weighs 7000 grains.
namespace Avoirdupois {
Quantities::Mass const Pound  = 0.45359237 * SI::Kilogram;
Quantities::Mass const Ounce  = Pound / 16;
Quantities::Mass const Drachm = Pound / 256;
Quantities::Mass const Grain  = Pound / 7000;

Quantities::Mass const Stone = 14 * Pound;

// Imperial quarter, hundredweight and pound, yielding the 'long ton'.
namespace Imperial {
Quantities::Mass const Quarter       = 2 * Stone;
Quantities::Mass const Hundredweight = 4 * Quarter;
Quantities::Mass const Ton           = 20 * Hundredweight;
}
// US quarter, hundredweight and pound, yielding the 'short ton'.
namespace US {
Quantities::Mass const Quarter       = 25 * Pound;
Quantities::Mass const Hundredweight = 4 * Quarter;
Quantities::Mass const Ton           = 20 * Hundredweight;
}
}
}