#pragma once

#include "Quantities.hpp"
#include "NamedQuantities.hpp"
#include "SI.hpp"

namespace principia {
namespace Constants {
Quantities::Speed const SpeedOfLight =
  299792458 * (SI::Metre / SI::Second);
Quantities::Permeability const VacuumPermeability =
  4e-7*π * SI::Steradian * SI::Henry / SI::Metre;
Quantities::Permittivity const VacuumPermittivity =
  1 / (VacuumPermeability*SpeedOfLight.Pow<2>());
// We use the 2010 CODATA recommended values. We do not support uncertainties.
Quantities::AngularMomentum const ReducedPlanckConstant =
  1.054571726e-34 * SI::Joule * SI::Second / SI::Radian;
auto const GravitationalConstant =
  6.67384e-11 * SI::Newton * SI::Metre.Pow<2>() / SI::Kilogram.Pow<2>();

Quantities::Entropy const BoltzmannConstant =
  1.3806488e-23 * (SI::Joule / SI::Kelvin);
Quantities::Inverse<Quantities::Amount> const AvogadroConstant =
  6.02214129 * (1 / SI::Mole);

Quantities::Mass   const ElectronMass     = 9.10938291e-31 * SI::Kilogram;
Quantities::Mass   const ProtonMass       = 1.672621777e-27 * SI::Kilogram;
Quantities::Charge const ElementaryCharge = SI::ElectronVolt / SI::Volt;

Quantities::Dimensionless const FineStructureConstant = 7.2973525698e-3;

Quantities::Acceleration const StandardGravity = 9.80665 * SI::Metre /
                                                 SI::Second.Pow<2>();
}  // namespace Constants
}  // namespace principia
