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
template<typename D> Quantity<D> Yotta(Quantity<D>);
template<typename D> Quantity<D> Zetta(Quantity<D>);
template<typename D> Quantity<D> Exa(Quantity<D>);
template<typename D> Quantity<D> Peta(Quantity<D>);
template<typename D> Quantity<D> Tera(Quantity<D>);
template<typename D> Quantity<D> Giga(Quantity<D>);
template<typename D> Quantity<D> Mega(Quantity<D>);
template<typename D> Quantity<D> Kilo(Quantity<D>);

template<typename D> Quantity<D> Hecto(Quantity<D>);
template<typename D> Quantity<D> Deca(Quantity<D>);

template<typename D> Quantity<D> Deci(Quantity<D>);
template<typename D> Quantity<D> Centi(Quantity<D>);

template<typename D> Quantity<D> Milli(Quantity<D>);
template<typename D> Quantity<D> Micro(Quantity<D>);
template<typename D> Quantity<D> Nano(Quantity<D>);
template<typename D> Quantity<D> Pico(Quantity<D>);
template<typename D> Quantity<D> Femto(Quantity<D>);
template<typename D> Quantity<D> Atto(Quantity<D>);
template<typename D> Quantity<D> Zepto(Quantity<D>);
template<typename D> Quantity<D> Yocto(Quantity<D>);

// SI base units
// From the BIPM's SI brochure 8, section 2.1.2, table 1,
// http://www.bipm.org/en/si/si_brochure/chapter2/2-1/.
Length            const Metre    = SIUnit<Length>();
Mass              const Kilogram = SIUnit<Mass>();
Time              const Second   = SIUnit<Time>();
Current           const Ampere   = SIUnit<Current>();
Temperature       const Kelvin   = SIUnit<Temperature>();
Amount            const Mole     = SIUnit<Amount>();
LuminousIntensity const Candela  = SIUnit<LuminousIntensity>();
// Nonstandard.
Winding const Cycle = SIUnit<Winding>();
// Not base units in the SI. We make these quantities rather than units as they
// are natural.
Angle const Radian = SIUnit<Angle>();
SolidAngle const Steradian = SIUnit<SolidAngle>();

// Gram, for use with prefixes.
Mass const Gram = 1e-3 * Kilogram;

// Coherent derived units in the SI with special names and symbols
// From the BIPM's SI brochure 8, section 2.2.2, table 3,
// http://www.bipm.org/en/si/si_brochure/chapter2/2-2/table3.html.
// We exclude the Becquerel, Gray and Sievert as they are weakly typed.
// The Celsius only really makes sense as an affine temperature and is not taken
// care of here.
// Note the nonstandard definition of the Hertz, with a dimensionful cycle.

// The uno was proposed but never accepted.
double              const Uno     = 1;
Frequency           const Hertz   = Cycle / Second;
Force               const Newton  = Metre * Kilogram / (Second * Second);
Pressure            const Pascal  = Newton / (Metre * Metre);
Energy              const Joule   = Newton * Metre;
Power               const Watt    = Joule / Second;
Charge              const Coulomb = Ampere * Second;
Voltage             const Volt    = Watt / Ampere;
Capacitance         const Farad   = Coulomb / Volt;
Resistance          const Ohm     = Volt / Ampere;
Conductance         const Siemens = Ampere / Volt;
MagneticFlux        const Weber   = Volt * Second;
MagneticFluxDensity const Tesla   = Weber / (Metre * Metre);
Inductance          const Henry   = Weber / Ampere;
LuminousFlux        const Lumen   = Candela * Steradian;
CatalyticActivity   const Katal   = Mole / Second;

// Non-SI units accepted for use with the SI
// From the BIPM's SI brochure 8, section 4.1, table 6,
// http://www.bipm.org/en/si/si_brochure/chapter4/table6.html
Time const Minute = 60 * Second;
Time const Hour   = 60 * Minute;
Time const Day    = 24 * Hour;

Angle  const Degree    = π / 180 * Radian;
Angle  const ArcMinute = π / 10800 * Radian;
Angle  const ArcSecond = π / 648000 * Radian;
Area   const Hectare   = 1e4 * Metre * Metre;
Volume const Litre     = Pow<3>(Deci(Metre));
Mass   const Tonne     = 1e3 * Kilogram;

// Non-SI units whose values must be obtained experimentally
// From the BIPM's SI brochure 8, section 4.1, table 7,
// Units accepted for use with the SI.
Energy const ElectronVolt     = 1.602176565e-19 * Joule;
Mass   const Dalton           = 1.660538921e-27 * Kilogram;
Length const AstronomicalUnit = 149597870700 * si::Metre;

}  // namespace si
}  // namespace quantities
}  // namespace principia

#include "quantities/si_body.hpp"
