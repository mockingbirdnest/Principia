#pragma once

#include "quantities/constants.hpp"
#include "quantities/elementary_functions.hpp"
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

constexpr Mass Pound  = 0.45359237 * si::Kilogram;
constexpr Mass Ounce  = Pound / 16;
constexpr Mass Drachm = Pound / 256;
constexpr Mass Grain  = Pound / 7000;

constexpr Mass Stone = 14 * Pound;
// Imperial quarter, hundredweight and pound, yielding the 'long ton'.
constexpr Mass Quarter       = 2 * Stone;
constexpr Mass Hundredweight = 4 * Quarter;
constexpr Mass Ton           = 20 * Hundredweight;

constexpr Length Yard = 0.9144 * si::Metre;
constexpr Length Foot = Yard / 3;
constexpr Length Inch = Foot / 12;
constexpr Length Thou = Inch / 1000;

constexpr Length Chain   = 22 * Yard;
constexpr Length Furlong = 10 * Chain;
constexpr Length Mile    = 8 * Furlong;
constexpr Length League  = 3 * Mile;

constexpr Length Link = Chain / 100;
constexpr Length Rod  = Chain / 4;

namespace admiralty {
constexpr Length NauticalMile = 6080 * Foot;
constexpr Length Cable        = NauticalMile / 10;
constexpr Length Fathom       = Cable / 100;
}  // namespace admiralty

constexpr Area Perch = Pow<2>(Rod);
constexpr Area Rood  = Furlong * Rod;
constexpr Area Acre  = Furlong * Chain;

constexpr Volume FluidOunce = 28.4130625 * si::Milli(si::Litre);
constexpr Volume Gill       = 5 * FluidOunce;
constexpr Volume Pint       = 4 * Gill;
constexpr Volume Quart      = 2 * Pint;
constexpr Volume Gallon     = 4 * Quart;

constexpr Force    PoundForce         = Pound * constants::StandardGravity;
constexpr Power    HorsePower         = 550 * PoundForce * Foot / si::Second;
constexpr Pressure PoundPerSquareInch = PoundForce / Pow<2>(Inch);

}  // namespace uk
}  // namespace quantities
}  // namespace principia
