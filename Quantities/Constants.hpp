#pragma once

#include "Quantities/NamedQuantities.hpp"
#include "Quantities/Quantities.hpp"
#include "Quantities/SI.hpp"

namespace principia {
namespace constants {
quantities::Speed const SpeedOfLight =
  299792458 * (si::Metre / si::Second);
quantities::Permeability const VacuumPermeability =
  4e-7*π * si::Steradian * si::Henry / si::Metre;
quantities::Permittivity const VacuumPermittivity =
  1 / (VacuumPermeability*SpeedOfLight.Pow<2>());
// We use the 2010 CODATA recommended values. We do not support uncertainties.
quantities::AngularMomentum const ReducedPlanckConstant =
  1.054571726e-34 * si::Joule * si::Second / si::Radian;
auto const GravitationalConstant =
  6.67384e-11 * si::Newton * si::Metre.Pow<2>() / si::Kilogram.Pow<2>();

quantities::Entropy const BoltzmannConstant =
  1.3806488e-23 * (si::Joule / si::Kelvin);
quantities::Inverse<quantities::Amount> const AvogadroConstant =
  6.02214129 * (1 / si::Mole);

quantities::Mass   const ElectronMass     = 9.10938291e-31 * si::Kilogram;
quantities::Mass   const ProtonMass       = 1.672621777e-27 * si::Kilogram;
quantities::Charge const ElementaryCharge = si::ElectronVolt / si::Volt;

quantities::Dimensionless const FineStructureConstant = 7.2973525698e-3;

quantities::Acceleration const StandardGravity = 9.80665 * si::Metre /
                                                 si::Second.Pow<2>();
}  // namespace constants
}  // namespace principia
