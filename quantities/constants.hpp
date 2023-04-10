#pragma once

#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace quantities {
namespace _constants {
namespace internal {

using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

// Defining constants for the SI, by resolution A of the 26th meeting of the
// CGPM.
constexpr Speed SpeedOfLight = 299'792'458 * (Metre / Second);
constexpr Action PlanckConstant = 6.626'070'15e-34 * Joule * Second;
constexpr Charge ElementaryCharge = 1.602'176'634e-19 * Coulomb;
constexpr Entropy BoltzmannConstant =  1.380'649e-23 * (Joule / Kelvin);
constexpr Inverse<Amount> AvogadroConstant = 6.022'140'76e23 * (1 / Mole);

// We use the 2014 CODATA recommended values.  We do not support uncertainties.
constexpr double FineStructureConstant = 7.297'352'5664e-3;
constexpr Quotient<GravitationalParameter, Mass> GravitationalConstant =
    6.674'08e-11 * Newton * Pow<2>(Metre) / Pow<2>(Kilogram);

constexpr Mass ElectronMass = 9.109'382'91e-31 * Kilogram;
constexpr Mass ProtonMass   = 1.672'621'777e-27 * Kilogram;

constexpr Acceleration StandardGravity =
    9.80665 * Metre / Pow<2>(Second);

constexpr Permeability VacuumPermeability =
    2 * FineStructureConstant * PlanckConstant * Steradian /
    (Pow<2>(ElementaryCharge) * SpeedOfLight);
constexpr Permittivity VacuumPermittivity =
    Pow<2>(ElementaryCharge) /
    (2 * FineStructureConstant * PlanckConstant * SpeedOfLight * Steradian);

// Units derived from the defining constants.
constexpr Energy ElectronVolt = ElementaryCharge * Volt;
constexpr Mass Dalton         = (Gram / Mole) / AvogadroConstant;

}  // namespace internal

using internal::AvogadroConstant;
using internal::BoltzmannConstant;
using internal::Dalton;
using internal::ElectronMass;
using internal::ElectronVolt;
using internal::ElementaryCharge;
using internal::FineStructureConstant;
using internal::GravitationalConstant;
using internal::PlanckConstant;
using internal::ProtonMass;
using internal::SpeedOfLight;
using internal::StandardGravity;
using internal::VacuumPermeability;
using internal::VacuumPermittivity;

}  // namespace _constants
}  // namespace quantities
}  // namespace principia
