#pragma once

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
Length            constexpr Metre    = SIUnit<Length>();
Mass              constexpr Kilogram = SIUnit<Mass>();
Time              constexpr Second   = SIUnit<Time>();
Current           constexpr Ampere   = SIUnit<Current>();
Temperature       constexpr Kelvin   = SIUnit<Temperature>();
Amount            constexpr Mole     = SIUnit<Amount>();
LuminousIntensity constexpr Candela  = SIUnit<LuminousIntensity>();
// Nonstandard.
Winding constexpr Cycle = SIUnit<Winding>();
// Not base units in the SI. We make these quantities rather than units as they
// are natural.
Angle constexpr Radian = SIUnit<Angle>();
SolidAngle constexpr Steradian = SIUnit<SolidAngle>();

// Gram, for use with prefixes.
Mass constexpr Gram = 1e-3 * Kilogram;

// Coherent derived units in the SI with special names and symbols
// From the BIPM's SI brochure 8, section 2.2.2, table 3,
// http://www.bipm.org/en/si/si_brochure/chapter2/2-2/table3.html.
// We exclude the Becquerel, Gray and Sievert as they are weakly typed.
// The Celsius only really makes sense as an affine temperature and is not taken
// care of here.
// Note the nonstandard definition of the Hertz, with a dimensionful cycle.

// The uno was proposed but never accepted.
double              constexpr Uno     = 1;
Frequency           constexpr Hertz   = Cycle / Second;
Force               constexpr Newton  = Metre * Kilogram / (Second * Second);
Pressure            constexpr Pascal  = Newton / (Metre * Metre);
Energy              constexpr Joule   = Newton * Metre;
Power               constexpr Watt    = Joule / Second;
Charge              constexpr Coulomb = Ampere * Second;
Voltage             constexpr Volt    = Watt / Ampere;
Capacitance         constexpr Farad   = Coulomb / Volt;
Resistance          constexpr Ohm     = Volt / Ampere;
Conductance         constexpr Siemens = Ampere / Volt;
MagneticFlux        constexpr Weber   = Volt * Second;
MagneticFluxDensity constexpr Tesla   = Weber / (Metre * Metre);
Inductance          constexpr Henry   = Weber / Ampere;
LuminousFlux        constexpr Lumen   = Candela * Steradian;
CatalyticActivity   constexpr Katal   = Mole / Second;

// Non-SI units accepted for use with the SI
// From the BIPM's SI brochure 8, section 4.1, table 6,
// http://www.bipm.org/en/si/si_brochure/chapter4/table6.html
Time constexpr Minute = 60 * Second;
Time constexpr Hour   = 60 * Minute;
Time constexpr Day    = 24 * Hour;

Angle  constexpr Degree    = π / 180 * Radian;
Angle  constexpr ArcMinute = π / 10800 * Radian;
Angle  constexpr ArcSecond = π / 648000 * Radian;
Area   constexpr Hectare   = 1e4 * Metre * Metre;
Volume constexpr Litre     = 1e-3 * Metre * Metre * Metre;
Mass   constexpr Tonne     = 1e3 * Kilogram;

// Non-SI units whose values must be obtained experimentally
// From the BIPM's SI brochure 8, section 4.1, table 7,
// Units accepted for use with the SI.
Energy constexpr ElectronVolt     = 1.602176565e-19 * Joule;
Mass   constexpr Dalton           = 1.660538921e-27 * Kilogram;
Length constexpr AstronomicalUnit = 149597870700 * si::Metre;

}  // namespace si
}  // namespace quantities
}  // namespace principia

#include "quantities/si_body.hpp"
