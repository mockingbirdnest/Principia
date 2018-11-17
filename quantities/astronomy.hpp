
#pragma once

#include "quantities/constants.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace quantities {

// This namespace contains units commonly used in astronomy that are not
// accepted for use with the SI.  See Résolution B2 "Re-définition de l’unité
// astronomique de longueur " adopted at the  XXVIII General Assembly of the IAU
// in 2012 and Résolution B3 "Sur les valeurs recommandées de constantes de
// conversion pour une sélection de propriétés solaires et planétaires" adopted
// at the XXIX General Assembly of the IAU in 2015.
namespace astronomy {

constexpr Mass   SolarMass             = 1.98855e30 * si::Kilogram;
constexpr Mass   JupiterMass           = 1.8986e27 * si::Kilogram;
constexpr Mass   EarthMass             = 5.9742e24 * si::Kilogram;
constexpr Time   JulianYear            = 365.25 * si::Day;
constexpr Length AstronomicalUnit      = 149597870700 * si::Metre;
constexpr Length EarthEquatorialRadius = 6.3781e6 * si::Metre;
constexpr Length SolarEquatorialRadius = 6.957e8 * si::Metre;
constexpr Length Parsec                = 648000 / π * AstronomicalUnit;
constexpr Length LightYear             = constants::SpeedOfLight * JulianYear;
constexpr Length LunarDistance         = 384400000 * si::Metre;

}  // namespace astronomy
}  // namespace quantities
}  // namespace principia
