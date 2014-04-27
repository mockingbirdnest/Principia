#pragma once

#include "Quantities.hpp"
#include "NamedQuantities.hpp"
#include "Numbers.hpp"

namespace Principia {
// This namespace contains the units and prefixes of the SI (except the
// Becquerel, Gray and Sievert), as well as the Non-SI units accepted for use
// with the SI.
namespace SI {
#pragma region Prefixes
template<typename D> Quantities::Quantity<D> Yotta(Quantities::Quantity<D>);
template<typename D> Quantities::Quantity<D> Zetta(Quantities::Quantity<D>);
template<typename D> Quantities::Quantity<D> Exa(Quantities::Quantity<D>);
template<typename D> Quantities::Quantity<D> Peta(Quantities::Quantity<D>);
template<typename D> Quantities::Quantity<D> Tera(Quantities::Quantity<D>);
template<typename D> Quantities::Quantity<D> Giga(Quantities::Quantity<D>);
template<typename D> Quantities::Quantity<D> Mega(Quantities::Quantity<D>);
template<typename D> Quantities::Quantity<D> Kilo(Quantities::Quantity<D>);

template<typename D> Quantities::Quantity<D> Hecto(Quantities::Quantity<D>);
template<typename D> Quantities::Quantity<D> Deca(Quantities::Quantity<D>);

template<typename D> Quantities::Quantity<D> Deci(Quantities::Quantity<D>);
template<typename D> Quantities::Quantity<D> Centi(Quantities::Quantity<D>);

template<typename D> Quantities::Quantity<D> Milli(Quantities::Quantity<D>);
template<typename D> Quantities::Quantity<D> Micro(Quantities::Quantity<D>);
template<typename D> Quantities::Quantity<D> Nano(Quantities::Quantity<D>);
template<typename D> Quantities::Quantity<D> Pico(Quantities::Quantity<D>);
template<typename D> Quantities::Quantity<D> Femto(Quantities::Quantity<D>);
template<typename D> Quantities::Quantity<D> Atto(Quantities::Quantity<D>);
template<typename D> Quantities::Quantity<D> Zepto(Quantities::Quantity<D>);
template<typename D> Quantities::Quantity<D> Yocto(Quantities::Quantity<D>);
#pragma endregion

#pragma region SI base units
// From the BIPM's SI brochure 8, section 2.1.2, table 1,
// http://www.bipm.org/en/si/si_brochure/chapter2/2-1/.
Quantities::Length      const Metre    = Quantities::Factories::Metres(1);
Quantities::Mass        const Kilogram = Quantities::Factories::Kilograms(1);
Quantities::Time        const Second   = Quantities::Factories::Seconds(1);
Quantities::Current     const Ampere   = Quantities::Factories::Amperes(1);
Quantities::Temperature const Kelvin   = Quantities::Factories::Kelvins(1);
Quantities::Amount      const Mole     = Quantities::Factories::Moles(1);
Quantities::LuminousIntensity
                        const Candela  = Quantities::Factories::Candelas(1);
// Nonstandard.
Quantities::Winding const Cycle = Quantities::Factories::Cycles(1);
// Not base units in the SI. We make these quantities rather than units as they
// are natural.
Quantities::Angle      const Radian    = Quantities::Factories::Radians(1);
Quantities::SolidAngle const Steradian = Quantities::Factories::Steradians(1);
#pragma endregion

// Gram, for use with prefixes.
Quantities::Mass const Gram = 1e-3 * Kilogram;

#pragma region Coherent derived units in the SI with special names and symbols
// From the BIPM's SI brochure 8, section 2.2.2, table 3,
// http://www.bipm.org/en/si/si_brochure/chapter2/2-2/table3.html.
// We exclude the Becquerel, Gray and Sievert as they are weakly typed.
// The Celsius only really makes sense as an affine temperature and is not taken
// care of here.
// Note the nonstandard definition of the Hertz, with a dimensionful cycle.

// The uno was proposed but never accepted.
Quantities::Dimensionless       const Uno = 1;
Quantities::Frequency           const Hertz   = Cycle / Second;
Quantities::Force               const Newton  = Metre * Kilogram / 
                                                (Second * Second);
Quantities::Pressure            const Pascal  = Newton / (Metre * Metre);
Quantities::Energy              const Joule   = Newton * Metre;
Quantities::Power               const Watt    = Joule / Second;
Quantities::Charge              const Coulomb = Ampere * Second;
Quantities::Voltage             const Volt    = Watt / Ampere;
Quantities::Capacitance         const Farad   = Coulomb / Volt;
Quantities::Resistance          const Ohm     = Volt / Ampere;
Quantities::Conductance         const Siemens = Ampere / Volt;
Quantities::MagneticFlux        const Weber   = Volt * Second;
Quantities::MagneticFluxDensity const Tesla   = Weber / (Metre * Metre);
Quantities::Inductance          const Henry   = Weber / Ampere;
Quantities::LuminousFlux        const Lumen   = Candela * Steradian;
Quantities::CatalyticActivity   const Katal   = Mole / Second;
#pragma endregion

#pragma region Non-SI units accepted for use with the SI
// From the BIPM's SI brochure 8, section 4.1, table 6,
// http://www.bipm.org/en/si/si_brochure/chapter4/table6.html
Quantities::Time const Minute = 60 * Second;
Quantities::Time const Hour = 60 * Minute;
Quantities::Time const Day = 24 * Hour;

Quantities::Angle  const Degree    = π / 180 * Radian;
Quantities::Angle  const ArcMinute = π / 10800 * Radian;
Quantities::Angle  const ArcSecond = π / 648000 * Radian;
Quantities::Area   const Hectare   = 1e4 * Metre * Metre;
Quantities::Volume const Litre     = Deci(Metre).Pow<3>();
Quantities::Mass   const Tonne     = 1e3 * Kilogram;
#pragma endregion
#pragma region Non-SI units whose values must be obtained experimentally
// From the BIPM's SI brochure 8, section 4.1, table 7,
// Units accepted for use with the SI.
Quantities::Energy const ElectronVolt     = 1.602176565e-19 * Joule;
Quantities::Mass   const Dalton           = 1.660538921e-27 * Kilogram;
Quantities::Length const AstronomicalUnit = 149597870700 * SI::Metre;
}  // namespace SI
}  // namespace Principia

#include "SI-body.hpp"
