// Units.h

#pragma once

#include "PhysicalQuantities.h"
#include "NamedQuantities.h"

namespace PhysicalQuantities {

#pragma region Base SI units
Unit<Length>      const Metre    = Unit<Length>(Metres(1.0));
Unit<Time>        const Second   = Unit<Time>(Seconds(1.0));
Unit<Mass>        const Kilogram = Unit<Mass>(Kilograms(1.0));
Unit<Temperature> const Kelvin   = Unit<Temperature>(Kelvins(1.0));
#pragma endregion
#pragma region General mechanics
Unit<Force>  const Newton = Metre * Kilogram / (Second * Second);
Unit<Energy> const Joule  = Newton * Metre;
#pragma endregion
#pragma region Thermodynamics
Unit<Pressure> const Pascal = Newton / (Metre * Metre);
Unit<Volume>   const Litre  = Unit<Volume>(1e-3 * (Metre * Metre * Metre));
#pragma endregion

}