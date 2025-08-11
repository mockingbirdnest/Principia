#pragma once

#include "quantities/constants.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace quantities {

// This namespace contains units commonly used in astronomy.
namespace _astronomy {
namespace internal {

using namespace principia::quantities::_constants;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

// Résolution B2 "Re-définition de l’unité astronomique de longueur" adopted
// at the XXVIIIth General Assembly of the IAU in 2012.
constexpr Length AstronomicalUnit = 149597870700 * Metre;

// See note 4 of Résolution B2 "Sur la recommandation du "point zéro" des
// échelles de magnitude bolométrique absolue et apparente" adopted at the
// XXIXth General Assembly of the IAU in 2015.
constexpr Length Parsec = 648000 / π * AstronomicalUnit;

// System of nominal solar and planetary conversion constants, Résolution B3
// "Sur les valeurs  recommandées de constantes de conversion pour une sélection
// de propriétés solaires et planétaires" adopted at the XXIXth General Assembly
// of the IAU in 2015.

// Solar conversion constants.
constexpr Length      SolarRadius               = 6.957e8 * Metre;
constexpr Irradiance  TotalSolarIrradiance      = 1361 * (Watt /
                                                          Pow<2>(Metre));
constexpr Power       SolarLuminosity           = 3.828e26 * Watt;
constexpr Temperature SolarEffectiveTemperature = 5772 * Kelvin;
constexpr GravitationalParameter SolarGravitationalParameter =
    1.327'124'4e20 * (Pow<3>(Metre) / Pow<2>(Second));

// Planetary conversion constants.
// “If equatorial vs. polar radius is not explicitly specified, it should be
// understood that nominal terrestrial [or jovian] radius refers specifically to
// [the nominal equatorial radius], following common usage.”
constexpr Length TerrestrialEquatorialRadius = 6.3781e6 * Metre;
constexpr Length TerrestrialPolarRadius      = 6.3568e6 * Metre;
constexpr Length JovianEquatorialRadius      = 7.1492e7 * Metre;
constexpr Length JovianPolarRadius           = 6.6854e7 * Metre;
constexpr GravitationalParameter TerrestrialGravitationalParameter =
    3.986'004e14 * (Pow<3>(Metre) / Pow<2>(Second));
constexpr GravitationalParameter JovianGravitationalParameter      =
    1.266'865'3e17 * (Pow<3>(Metre) / Pow<2>(Second));

constexpr Time   JulianYear = 365.25 * Day;
constexpr Length LightYear  = SpeedOfLight * JulianYear;

}  // namespace internal

using internal::AstronomicalUnit;
using internal::JovianEquatorialRadius;
using internal::JovianGravitationalParameter;
using internal::JovianPolarRadius;
using internal::JulianYear;
using internal::LightYear;
using internal::Parsec;
using internal::SolarEffectiveTemperature;
using internal::SolarGravitationalParameter;
using internal::SolarLuminosity;
using internal::SolarRadius;
using internal::TerrestrialEquatorialRadius;
using internal::TerrestrialGravitationalParameter;
using internal::TerrestrialPolarRadius;
using internal::TotalSolarIrradiance;

}  // namespace _astronomy
}  // namespace quantities
}  // namespace principia
