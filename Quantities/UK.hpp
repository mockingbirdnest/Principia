#pragma once

#include "Constants.hpp"
#include "NamedQuantities.hpp"
#include "Numbers.hpp"
#include "Quantities.hpp"
#include "SI.hpp"

namespace Principia {
// This namespace contains the imperial units as defined by the international
// yard and pound agreement, as well as the units of the English Engineering
// system.
namespace UK {
Quantities::Mass const Pound  = 0.45359237 * SI::Kilogram;
Quantities::Mass const Ounce  = Pound / 16;
Quantities::Mass const Drachm = Pound / 256;
Quantities::Mass const Grain  = Pound / 7000;

Quantities::Mass const Stone = 14 * Pound;
// Imperial quarter, hundredweight and pound, yielding the 'long ton'.
Quantities::Mass const Quarter       = 2 * Stone;
Quantities::Mass const Hundredweight = 4 * Quarter;
Quantities::Mass const Ton           = 20 * Hundredweight;

Quantities::Length const Yard = 0.9144 * SI::Metre;
Quantities::Length const Foot = Yard / 3;
Quantities::Length const Inch = Foot / 12;
Quantities::Length const Thou = Foot / 1000;

Quantities::Length const Chain   = 22 * Yard;
Quantities::Length const Furlong = 10 * Chain;
Quantities::Length const Mile    = 8 * Furlong;
Quantities::Length const League  = 3 * Mile;

Quantities::Length const Link = Chain / 100;
Quantities::Length const Rod  = Chain / 4;

namespace Admiralty {
Quantities::Length const NauticalMile = 6080 * Foot;
Quantities::Length const Cable        = NauticalMile / 10;
Quantities::Length const Fathom       = Cable / 100;
}  // namespace Admiralty

Quantities::Area const Perch = Rod.Pow<2>();
Quantities::Area const Rood  = Furlong * Rod;
Quantities::Area const Acre  = Furlong * Chain;

Quantities::Volume const FluidOunce = 28.4130625 * SI::Milli(SI::Litre);
Quantities::Volume const Gill       = 5 * FluidOunce;
Quantities::Volume const Pint       = 4 * Gill;
Quantities::Volume const Quart      = 2 * Pint;
Quantities::Volume const Gallon     = 4 * Quart;

Quantities::Force    const PoundForce = Pound * Constants::StandardGravity;
Quantities::Pressure const PoundPerSquareInch = PoundForce / Inch.Pow<2>();
Quantities::Power    const HorsePower = 550 * PoundForce * Foot / SI::Second;
}  // namespace UK
}  // namespace Principia
