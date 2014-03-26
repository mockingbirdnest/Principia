// Units.h

#pragma once

#include "PhysicalQuantities.h"
#include "NamedQuantities.h"

namespace PhysicalQuantities {
#pragma region SI base units
// From the BIPM's SI brochure 8, section 2.1.2, table 1,
// http://www.bipm.org/en/si/si_brochure/chapter2/2-1/.
Unit<Length>            const Metre     = Unit<Length>(Metres(1));
Unit<Mass>              const Kilogram  = Unit<Mass>(Kilograms(1));
Unit<Time>              const Second    = Unit<Time>(Seconds(1));
Unit<Current>           const Ampere    = Unit<Current>(Amperes(1));
Unit<Temperature>       const Kelvin    = Unit<Temperature>(Kelvins(1));
Unit<Amount>            const Mole      = Unit<Amount>(Moles(1));
Unit<LuminousIntensity> const Candela   = Unit<LuminousIntensity>(Candelas(1));
// Nonstandard.
Unit<Winding>           const Cycle     = Unit<Winding>(Cycles(1));
// Not base units in the SI. We make these quantities rather than units as they
// are natural.
Angle                   const Radian    = Radians(1);
SolidAngle              const Steradian = Steradians(1);
#pragma endregion
#pragma region Coherent derived units in the SI with special names and symbols
// From the BIPM's SI brochure 8, section 2.2.2, table 3,
// http://www.bipm.org/en/si/si_brochure/chapter2/2-2/table3.html.
// We exclude the Becquerel, Gray and Sievert as they are weakly typed.
// The Celsius only really makes sense as an affine temperature and is not taken
// care of here.
// Note the nonstandard definition of the Hertz, with a dimensionful cycle.

// The uno was proposed but never accepted. We use it because uno/mol is better 
// than writing rad/mol and allowing 1/mol is complicated.
Unit<Dimensionless> const Uno =
  Unit<Dimensionless>(Dimensionless(1));
// Dimensionful units.
Unit<Frequency>           const Hertz   = Cycle / Second;
Unit<Force>               const Newton  = Metre * Kilogram / (Second * Second);
Unit<Pressure>            const Pascal  = Newton / (Metre * Metre);
Unit<Energy>              const Joule   = Newton * Metre;
Unit<Power>               const Watt    = Joule / Second;
Unit<Charge>              const Coulomb = Ampere * Second;
Unit<Voltage>             const Volt    = Watt / Ampere;
Unit<Capacitance>         const Farad   = Coulomb / Volt;
Unit<Resistance>          const Ohm     = Volt / Ampere;
Unit<Conductance>         const Siemens = Ampere / Volt;
Unit<MagneticFlux>        const Weber   = Volt * Second;
Unit<MagneticFluxDensity> const Tesla   = Weber / (Metre * Metre);
Unit<Inductance>          const Henry   = Weber / Ampere;
Unit<LuminousFlux>        const Lumen   = Candela * Steradian;
Unit<CatalyticActivity>   const Katal   = Mole / Second;
#pragma endregion
#pragma region Non-SI units accepted for use with the SI
// From the BIPM's SI brochure 8, section 4.1, table 6,
// http://www.bipm.org/en/si/si_brochure/chapter4/table6.html
Unit<Time> const Minute = Unit<Time>(60 * Second);
Unit<Time> const Hour = Unit<Time>(60 * Minute);
Unit<Time> const Day = Unit<Time>(24 * Hour);
// TODO(robin): Move this declaration somewhere else.
double π = 3.14159265358979323846264338327950288419716939937511;
Unit<Angle>  const Degree    = Unit<Angle>(π / 180 * Radian);
Unit<Angle>  const ArcMinute = Unit<Angle>(π / 10800 * Radian);
Unit<Angle>  const ArcSecond = Unit<Angle>(π / 648000 * Radian);
Unit<Area>   const Hectare   = Unit<Area>(1e4 * (Metre * Metre));
Unit<Volume> const Litre     = Unit<Volume>(1e-3 * (Metre * Metre * Metre));
Unit<Mass>   const Tonne     = Unit<Mass>(1e3 * Kilogram);
#pragma endregion
}