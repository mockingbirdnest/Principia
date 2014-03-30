#pragma once

#include "Quantities.h"
#include "NamedQuantities.h"
#include "Numbers.h"
#include "SI.h"
#include "Constants.h"

namespace Principia {
// This namespace contains units commonly used in astronomy that are not
// accepted for use with the SI.
namespace Astronomy {
Quantities::Mass   const SolarMass   = 1.98855e30 * SI::Kilogram;
Quantities::Mass   const JupiterMass = 1.8986e27 * SI::Kilogram;
Quantities::Mass   const EarthMass   = 5.9742e24 * SI::Kilogram;
Quantities::Time   const JulianYear  = 365.25 * SI::Day;
Quantities::Length const Parsec      = 648000 / π * SI::AstronomicalUnit;
Quantities::Length const LightYear   = Constants::SpeedOfLight * JulianYear;
}
}
