// Units.h

#pragma once

#include "PhysicalQuantities.h"
#include "NamedQuantities.h"

namespace PhysicalQuantities {

#pragma region Base SI units
const Unit<Length> Metre = Unit<Length>(Metres(1.0));
const Unit<Time> Second = Unit<Time>(Seconds(1.0));
const Unit<Mass> Kilogram = Unit<Mass>(Kilograms(1.0));
const Unit<Temperature> Kelvin = Unit<Temperature>(Kelvins(1.0));
#pragma endregion
#pragma region General mechanics
const Unit<Force> Newton = Metre * Kilogram / (Second * Second);
const Unit<Energy> Joule = Newton * Metre;
#pragma endregion
#pragma region Thermodynamics
const Unit<Pressure> Pascal = Newton / (Metre * Metre);
const Unit<Volume> Litre = Unit<Volume>(1e-3 * (Metre * Metre * Metre));
#pragma endregion

}