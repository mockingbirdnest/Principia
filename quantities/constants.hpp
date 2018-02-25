
#pragma once

#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace quantities {
namespace constants {

constexpr Speed SpeedOfLight = 299792458 * (si::Metre / si::Second);
constexpr Permeability VacuumPermeability =
    4e-7 * π * si::Steradian * si::Henry / si::Metre;
constexpr Permittivity VacuumPermittivity =
    1 / (VacuumPermeability * Pow<2>(SpeedOfLight));
// We use the 2010 CODATA recommended values.  We do not support uncertainties.
constexpr AngularMomentum ReducedPlanckConstant =
    1.054571726e-34 * si::Joule * si::Second * si::Radian;
constexpr Quotient<GravitationalParameter, Mass> GravitationalConstant =
    6.67384e-11 * si::Newton * Pow<2>(si::Metre) / Pow<2>(si::Kilogram);
constexpr Entropy BoltzmannConstant = 1.3806488e-23 * (si::Joule / si::Kelvin);
constexpr Amount::Inverse AvogadroConstant = 6.02214129 * (1 / si::Mole);

constexpr Mass   ElectronMass     = 9.10938291e-31 * si::Kilogram;
constexpr Mass   ProtonMass       = 1.672621777e-27 * si::Kilogram;
constexpr Charge ElementaryCharge = si::ElectronVolt / si::Volt;

constexpr double FineStructureConstant = 7.2973525698e-3;

constexpr Acceleration StandardGravity =
    9.80665 * si::Metre / Pow<2>(si::Second);

}  // namespace constants
}  // namespace quantities
}  // namespace principia
