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

Mass constexpr Pound  = 0.45359237 * si::Kilogram;
Mass constexpr Ounce  = Pound / 16;
Mass constexpr Drachm = Pound / 256;
Mass constexpr Grain  = Pound / 7000;

Mass constexpr Stone = 14 * Pound;
// Imperial quarter, hundredweight and pound, yielding the 'long ton'.
Mass constexpr Quarter       = 2 * Stone;
Mass constexpr Hundredweight = 4 * Quarter;
Mass constexpr Ton           = 20 * Hundredweight;

Length constexpr Yard = 0.9144 * si::Metre;
Length constexpr Foot = Yard / 3;
Length constexpr Inch = Foot / 12;
Length constexpr Thou = Foot / 1000;

Length constexpr Chain   = 22 * Yard;
Length constexpr Furlong = 10 * Chain;
Length constexpr Mile    = 8 * Furlong;
Length constexpr League  = 3 * Mile;

Length constexpr Link = Chain / 100;
Length constexpr Rod  = Chain / 4;

namespace admiralty {
Length constexpr NauticalMile = 6080 * Foot;
Length constexpr Cable        = NauticalMile / 10;
Length constexpr Fathom       = Cable / 100;
}  // namespace admiralty

Area constexpr Perch = Pow<2>(Rod);
Area constexpr Rood  = Furlong * Rod;
Area constexpr Acre  = Furlong * Chain;

Volume constexpr FluidOunce = 28.4130625 * si::Milli(si::Litre);
Volume constexpr Gill       = 5 * FluidOunce;
Volume constexpr Pint       = 4 * Gill;
Volume constexpr Quart      = 2 * Pint;
Volume constexpr Gallon     = 4 * Quart;

Force    constexpr PoundForce         = Pound * constants::StandardGravity;
Power    constexpr HorsePower         = 550 * PoundForce * Foot / si::Second;
Pressure constexpr PoundPerSquareInch = PoundForce / Pow<2>(Inch);

}  // namespace uk
}  // namespace quantities
}  // namespace principia
