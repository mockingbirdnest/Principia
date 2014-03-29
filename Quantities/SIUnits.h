// Units.h

#pragma once

#include "Quantities.h"
#include "NamedQuantities.h"
#include "MathematicalConstants.h"

namespace Quantities {
#pragma region Prefixes
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
#pragma endregion

#pragma region SI base units
// From the BIPM's SI brochure 8, section 2.1.2, table 1,
// http://www.bipm.org/en/si/si_brochure/chapter2/2-1/.
Length            const Metre     = Factories::Metres(1);
Mass              const Kilogram  = Factories::Kilograms(1);
Time              const Second    = Factories::Seconds(1);
Current           const Ampere    = Factories::Amperes(1);
Temperature       const Kelvin    = Factories::Kelvins(1);
Amount            const Mole      = Factories::Moles(1);
LuminousIntensity const Candela   = Factories::Candelas(1);
// Nonstandard.
Winding           const Cycle     = Factories::Cycles(1);
// Not base units in the SI. We make these quantities rather than units as they
// are natural.
Angle             const Radian    = Factories::Radians(1);
SolidAngle        const Steradian = Factories::Steradians(1);
#pragma endregion

#pragma region Coherent derived units in the SI with special names and symbols
// From the BIPM's SI brochure 8, section 2.2.2, table 3,
// http://www.bipm.org/en/si/si_brochure/chapter2/2-2/table3.html.
// We exclude the Becquerel, Gray and Sievert as they are weakly typed.
// The Celsius only really makes sense as an affine temperature and is not taken
// care of here.
// Note the nonstandard definition of the Hertz, with a dimensionful cycle.

// The uno was proposed but never accepted.
Dimensionless       const Uno     = 1;
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
#pragma endregion

#pragma region Non-SI units accepted for use with the SI
// From the BIPM's SI brochure 8, section 4.1, table 6,
// http://www.bipm.org/en/si/si_brochure/chapter4/table6.html
Time const Minute = 60 * Second;
Time const Hour = 60 * Minute;
Time const Day = 24 * Hour;

Angle  const Degree    = π / 180 * Radian;
Angle  const ArcMinute = π / 10800 * Radian;
Angle  const ArcSecond = π / 648000 * Radian;
Area   const Hectare   = 1e4 * Metre * Metre;
Volume const Litre     = 1e-3 * Metre.Pow<3>();
Mass   const Tonne     = 1e3 * Kilogram;
#pragma endregion
}

#include "SIUnits.inl"
