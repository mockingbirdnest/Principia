#pragma once

#include "quantities/constants.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace quantities {

// This namespace contains the imperial units as defined by the international
// yard and pound agreement, as well as the units of the English Engineering
// system.
namespace uk {

Mass const Pound  = 0.45359237 * si::Kilogram;
Mass const Ounce  = Pound / 16;
Mass const Drachm = Pound / 256;
Mass const Grain  = Pound / 7000;

Mass const Stone = 14 * Pound;
// Imperial quarter, hundredweight and pound, yielding the 'long ton'.
Mass const Quarter       = 2 * Stone;
Mass const Hundredweight = 4 * Quarter;
Mass const Ton           = 20 * Hundredweight;

Length const Yard = 0.9144 * si::Metre;
Length const Foot = Yard / 3;
Length const Inch = Foot / 12;
Length const Thou = Foot / 1000;

Length const Chain   = 22 * Yard;
Length const Furlong = 10 * Chain;
Length const Mile    = 8 * Furlong;
Length const League  = 3 * Mile;

Length const Link = Chain / 100;
Length const Rod  = Chain / 4;

namespace admiralty {
Length const NauticalMile = 6080 * Foot;
Length const Cable        = NauticalMile / 10;
Length const Fathom       = Cable / 100;
}  // namespace admiralty

Area const Perch = Pow<2>(Rod);
Area const Rood  = Furlong * Rod;
Area const Acre  = Furlong * Chain;

Volume const FluidOunce = 28.4130625 * si::Milli(si::Litre);
Volume const Gill       = 5 * FluidOunce;
Volume const Pint       = 4 * Gill;
Volume const Quart      = 2 * Pint;
Volume const Gallon     = 4 * Quart;

Force    const PoundForce         = Pound * constants::StandardGravity;
Power    const HorsePower         = 550 * PoundForce * Foot / si::Second;
Pressure const PoundPerSquareInch = PoundForce / Pow<2>(Inch);

}  // namespace uk
}  // namespace quantities
}  // namespace principia
