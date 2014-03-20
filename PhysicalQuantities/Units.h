// Units.h

#pragma once

#include "PhysicalQuantities.h"

namespace PhysicalQuantities {

#pragma region Base SI units
inline Length Metres(double number) { return Length(number); }
inline Time Seconds(double number) { return Time(number); }
inline Mass Kilograms(double number) { return Mass(number); }
inline Temperature Kelvins(double number) { return Temperature(number); }
const Unit<Length> Metre = Unit<Length>(Metres(1.0));
const Unit<Time> Second = Unit<Time>(Seconds(1.0));
const Unit<Mass> Kilogram = Unit<Mass>(Kilograms(1.0));
const Unit<Temperature> Kelvin = Unit<Temperature>(Kelvins(1.0));
#pragma endregion
#pragma region Further units for base quantities
inline Temperature Celsius(double number) {
  return Kelvins(number) + Kelvins(273.15);
}
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