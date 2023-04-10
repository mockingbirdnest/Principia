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
namespace _uk {
namespace internal {

using namespace principia::quantities::_constants;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

constexpr Mass Pound  = 0.45359237 * Kilogram;
constexpr Mass Ounce  = Pound / 16;
constexpr Mass Drachm = Pound / 256;
constexpr Mass Grain  = Pound / 7000;

constexpr Mass Stone = 14 * Pound;
// Imperial quarter, hundredweight and pound, yielding the 'long ton'.
constexpr Mass Quarter       = 2 * Stone;
constexpr Mass Hundredweight = 4 * Quarter;
constexpr Mass Ton           = 20 * Hundredweight;

constexpr Length Yard = 0.9144 * Metre;
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

constexpr Volume FluidOunce = 28.4130625 * Milli(Litre);
constexpr Volume Gill       = 5 * FluidOunce;
constexpr Volume Pint       = 4 * Gill;
constexpr Volume Quart      = 2 * Pint;
constexpr Volume Gallon     = 4 * Quart;

constexpr Force    PoundForce         = Pound * StandardGravity;
constexpr Power    HorsePower         = 550 * PoundForce * Foot / Second;
constexpr Pressure PoundPerSquareInch = PoundForce / Pow<2>(Inch);

}  // namespace internal

using internal::Acre;
using internal::Chain;
using internal::Drachm;
using internal::FluidOunce;
using internal::Foot;
using internal::Furlong;
using internal::Gallon;
using internal::Gill;
using internal::Grain;
using internal::HorsePower;
using internal::Hundredweight;
using internal::Inch;
using internal::League;
using internal::Link;
using internal::Mile;
using internal::Ounce;
using internal::Perch;
using internal::Pint;
using internal::Pound;
using internal::PoundForce;
using internal::PoundPerSquareInch;
using internal::Quart;
using internal::Quarter;
using internal::Rod;
using internal::Rood;
using internal::Stone;
using internal::Thou;
using internal::Ton;
using internal::Yard;
namespace admiralty = internal::admiralty;

}  // namespace _uk
}  // namespace quantities
}  // namespace principia
