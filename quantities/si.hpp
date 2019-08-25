
#pragma once

#include <string>

#include "quantities/named_quantities.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace quantities {

// This namespace contains the units and prefixes of the SI (except the
// Becquerel, Gray and Sievert), as well as the Non-SI units accepted for use
// with the SI.
namespace si {

// Prefixes
template<typename D> constexpr Quantity<D> Yotta(Quantity<D>);
template<typename D> constexpr Quantity<D> Zetta(Quantity<D>);
template<typename D> constexpr Quantity<D> Exa(Quantity<D>);
template<typename D> constexpr Quantity<D> Peta(Quantity<D>);
template<typename D> constexpr Quantity<D> Tera(Quantity<D>);
template<typename D> constexpr Quantity<D> Giga(Quantity<D>);
template<typename D> constexpr Quantity<D> Mega(Quantity<D>);
template<typename D> constexpr Quantity<D> Kilo(Quantity<D>);

template<typename D> constexpr Quantity<D> Hecto(Quantity<D>);
template<typename D> constexpr Quantity<D> Deca(Quantity<D>);

template<typename D> constexpr Quantity<D> Deci(Quantity<D>);
template<typename D> constexpr Quantity<D> Centi(Quantity<D>);

template<typename D> constexpr Quantity<D> Milli(Quantity<D>);
template<typename D> constexpr Quantity<D> Micro(Quantity<D>);
template<typename D> constexpr Quantity<D> Nano(Quantity<D>);
template<typename D> constexpr Quantity<D> Pico(Quantity<D>);
template<typename D> constexpr Quantity<D> Femto(Quantity<D>);
template<typename D> constexpr Quantity<D> Atto(Quantity<D>);
template<typename D> constexpr Quantity<D> Zepto(Quantity<D>);
template<typename D> constexpr Quantity<D> Yocto(Quantity<D>);

// SI base units
// From the BIPM's SI brochure 8, section 2.1.2, table 1,
// http://www.bipm.org/en/si/si_brochure/chapter2/2-1/.
constexpr Length            Metre    = SIUnit<Length>();
constexpr Mass              Kilogram = SIUnit<Mass>();
constexpr Time              Second   = SIUnit<Time>();
constexpr Current           Ampere   = SIUnit<Current>();
constexpr Temperature       Kelvin   = SIUnit<Temperature>();
constexpr Amount            Mole     = SIUnit<Amount>();
constexpr LuminousIntensity Candela  = SIUnit<LuminousIntensity>();
// Not a base unit in the SI.
constexpr Angle Radian = SIUnit<Angle>();

// Gram, for use with prefixes.
constexpr Mass Gram = 1e-3 * Kilogram;

// Coherent derived units in the SI with special names and symbols
// From the BIPM's SI brochure 8, section 2.2.2, table 3,
// http://www.bipm.org/en/si/si_brochure/chapter2/2-2/table3.html.
// We exclude the Becquerel, Gray and Sievert as they are weakly typed.
// The Celsius only really makes sense as an affine temperature and is not taken
// care of here.

constexpr SolidAngle          Steradian = Radian * Radian;
constexpr Frequency           Hertz     = 1 / Second;
constexpr Force               Newton    = Metre * Kilogram / (Second * Second);
constexpr Pressure            Pascal    = Newton / (Metre * Metre);
constexpr Energy              Joule     = Newton * Metre;
constexpr Power               Watt      = Joule / Second;
constexpr Charge              Coulomb   = Ampere * Second;
constexpr Voltage             Volt      = Watt / Ampere;
constexpr Capacitance         Farad     = Coulomb / Volt;
constexpr Resistance          Ohm       = Volt / Ampere;
constexpr Conductance         Siemens   = Ampere / Volt;
constexpr MagneticFlux        Weber     = Volt * Second;
constexpr MagneticFluxDensity Tesla     = Weber / (Metre * Metre);
constexpr Inductance          Henry     = Weber / Ampere;
constexpr LuminousFlux        Lumen     = Candela * Steradian;
constexpr CatalyticActivity   Katal     = Mole / Second;

// Non-SI units accepted for use with the SI
// From the BIPM's SI brochure 8, section 4.1, table 6,
// http://www.bipm.org/en/si/si_brochure/chapter4/table6.html
constexpr Time Minute = 60 * Second;
constexpr Time Hour   = 60 * Minute;
constexpr Time Day    = 24 * Hour;

constexpr Angle  Degree    = π / 180 * Radian;
constexpr Angle  ArcMinute = π / 10800 * Radian;
constexpr Angle  ArcSecond = π / 648000 * Radian;
constexpr Area   Hectare   = 1e4 * Metre * Metre;
constexpr Volume Litre     = 1e-3 * Metre * Metre * Metre;
constexpr Mass   Tonne     = 1e3 * Kilogram;

}  // namespace si
}  // namespace quantities
}  // namespace principia

#include "quantities/si_body.hpp"
