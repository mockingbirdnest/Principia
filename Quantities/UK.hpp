#pragma once

#include "Quantities/Constants.hpp"
#include "Quantities/NamedQuantities.hpp"
#include "Quantities/Numbers.hpp"
#include "Quantities/Quantities.hpp"
#include "Quantities/SI.hpp"

namespace principia {
// This namespace contains the imperial units as defined by the international
// yard and pound agreement, as well as the units of the English Engineering
// system.
namespace uk {
quantities::Mass const Pound  = 0.45359237 * si::Kilogram;
quantities::Mass const Ounce  = Pound / 16;
quantities::Mass const Drachm = Pound / 256;
quantities::Mass const Grain  = Pound / 7000;

quantities::Mass const Stone = 14 * Pound;
// Imperial quarter, hundredweight and pound, yielding the 'long ton'.
quantities::Mass const Quarter       = 2 * Stone;
quantities::Mass const Hundredweight = 4 * Quarter;
quantities::Mass const Ton           = 20 * Hundredweight;

quantities::Length const Yard = 0.9144 * si::Metre;
quantities::Length const Foot = Yard / 3;
quantities::Length const Inch = Foot / 12;
quantities::Length const Thou = Foot / 1000;

quantities::Length const Chain   = 22 * Yard;
quantities::Length const Furlong = 10 * Chain;
quantities::Length const Mile    = 8 * Furlong;
quantities::Length const League  = 3 * Mile;

quantities::Length const Link = Chain / 100;
quantities::Length const Rod  = Chain / 4;

namespace admiralty {
quantities::Length const NauticalMile = 6080 * Foot;
quantities::Length const Cable        = NauticalMile / 10;
quantities::Length const Fathom       = Cable / 100;
}  // namespace admiralty

quantities::Area const Perch = Rod.Pow<2>();
quantities::Area const Rood  = Furlong * Rod;
quantities::Area const Acre  = Furlong * Chain;

quantities::Volume const FluidOunce = 28.4130625 * si::Milli(si::Litre);
quantities::Volume const Gill       = 5 * FluidOunce;
quantities::Volume const Pint       = 4 * Gill;
quantities::Volume const Quart      = 2 * Pint;
quantities::Volume const Gallon     = 4 * Quart;

quantities::Force    const PoundForce = Pound * constants::StandardGravity;
quantities::Pressure const PoundPerSquareInch = PoundForce / Inch.Pow<2>();
quantities::Power    const HorsePower = 550 * PoundForce * Foot / si::Second;
}  // namespace uk
}  // namespace principia
