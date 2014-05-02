#pragma once

#include "Quantities.hpp"
#include "NamedQuantities.hpp"
#include "Numbers.hpp"
#include "SI.hpp"
#include "Constants.hpp"

namespace principia {
// This namespace contains units commonly used in astronomy that are not
// accepted for use with the SI.
namespace Astronomy {
Quantities::Mass   const SolarMass     = 1.98855e30 * SI::Kilogram;
Quantities::Mass   const JupiterMass   = 1.8986e27 * SI::Kilogram;
Quantities::Mass   const EarthMass     = 5.9742e24 * SI::Kilogram;
Quantities::Time   const JulianYear    = 365.25 * SI::Day;
Quantities::Length const Parsec        = 648000 / π * SI::AstronomicalUnit;
Quantities::Length const LightYear     = Constants::SpeedOfLight * JulianYear;
Quantities::Length const LunarDistance = 384400000 * SI::Metre;
}  // namespace Astronomy
}  // namespace principia
