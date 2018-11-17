
#pragma once

#include "quantities/constants.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace quantities {

// This namespace contains units commonly used in astronomy.
namespace astronomy {

// Résolution B2 "Re-définition de l’unité astronomique de longueur" adopted
// at the XXVIIIth General Assembly of the IAU in 2012.
constexpr Length AstronomicalUnit      = 149597870700 * si::Metre;

// See note 4 of Résolution B2 "Sur la recommandation du "point zéro" des
// échelles de magnitude bolométrique absolue et apparente" adopted at the
// XXIXth General Assembly of the IAU in 2015.
constexpr Length Parsec                = 648000 / π * AstronomicalUnit;

// System of nominal solar and planetary conversion constants, Résolution B3
// "Sur les valeurs  recommandées de constantes de conversion pour une sélection
// de propriétés solaires et planétaires" adopted at the XXIXth General Assembly
// of the IAU in 2015.
constexpr Mass   SolarMass             = 1.98855e30 * si::Kilogram;
constexpr Mass   JupiterMass           = 1.8986e27 * si::Kilogram;
constexpr Mass   EarthMass             = 5.9742e24 * si::Kilogram;
constexpr Length EarthEquatorialRadius = 6.3781e6 * si::Metre;
constexpr Length SolarRadius           = 6.957e8 * si::Metre;

constexpr Time   JulianYear            = 365.25 * si::Day;
constexpr Length LightYear             = constants::SpeedOfLight * JulianYear;

}  // namespace astronomy
}  // namespace quantities
}  // namespace principia
