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
namespace _si {
namespace internal {

// Returns the base or derived SI unit of |Q|.
// For instance, |si::Unit<Action>() == Joule * Second|.
template<typename Q>
constexpr Q Unit = _quantities::internal::SIUnit<Q>();
template<>
inline constexpr double Unit<double> = 1;

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
constexpr Length            Metre    = Unit<Length>;
constexpr Mass              Kilogram = Unit<Mass>;
constexpr Time              Second   = Unit<Time>;
constexpr Current           Ampere   = Unit<Current>;
constexpr Temperature       Kelvin   = Unit<Temperature>;
constexpr Amount            Mole     = Unit<Amount>;
constexpr LuminousIntensity Candela  = Unit<LuminousIntensity>;
// Not a base unit in the SI.
constexpr Angle Radian = Unit<Angle>;

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

}  // namespace internal

using internal::Ampere;
using internal::ArcMinute;
using internal::ArcSecond;
using internal::Atto;
using internal::Candela;
using internal::Centi;
using internal::Coulomb;
using internal::Day;
using internal::Deca;
using internal::Deci;
using internal::Degree;
using internal::Exa;
using internal::Farad;
using internal::Femto;
using internal::Giga;
using internal::Gram;
using internal::Hectare;
using internal::Hecto;
using internal::Henry;
using internal::Hertz;
using internal::Hour;
using internal::Joule;
using internal::Katal;
using internal::Kelvin;
using internal::Kilo;
using internal::Kilogram;
using internal::Litre;
using internal::Lumen;
using internal::Mega;
using internal::Metre;
using internal::Micro;
using internal::Milli;
using internal::Minute;
using internal::Mole;
using internal::Nano;
using internal::Newton;
using internal::Ohm;
using internal::Pascal;
using internal::Peta;
using internal::Pico;
using internal::Radian;
using internal::Second;
using internal::Siemens;
using internal::Steradian;
using internal::Tera;
using internal::Tesla;
using internal::Tonne;
using internal::Unit;
using internal::Volt;
using internal::Watt;
using internal::Weber;
using internal::Yocto;
using internal::Yotta;
using internal::Zepto;
using internal::Zetta;
namespace si = _si;

}  // namespace _si
}  // namespace quantities
}  // namespace principia

namespace principia::quantities {
using namespace principia::quantities::_si;
}  // namespace principia::quantities

#include "quantities/si_body.hpp"
