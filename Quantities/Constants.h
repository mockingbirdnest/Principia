// Constants.h

#pragma once

#include "Quantities.h"
#include "NamedQuantities.h"
#include "SI.h"

namespace Principia {
namespace Constants {
auto const SpeedOfLight       = 299792458 * (SI::Metre / SI::Second);
auto const VacuumPermeability = 4e-7*π * SI::Steradian * SI::Henry / SI::Metre;
auto const VacuumPermittivity = 1 / (VacuumPermeability*SpeedOfLight.Pow<2>());
// We use the 2010 CODATA recommended values. We do not support uncertainties.
auto const ReducedPlanckConstant = 1.054571726e-34 * (SI::Joule * SI::Second);
auto const GravitationalConstant =
  6.67384e-11 * SI::Newton * SI::Metre.Pow<2>() / SI::Kilogram.Pow<2>();

auto const BoltzmannConstant = 1.3806488e-23 * (SI::Joule / SI::Kelvin);
auto const AvogadroConstant  = 6.02214129 * (1 / SI::Mole);

auto const ElectronMass     = 9.10938291e-31 * SI::Kilogram;
auto const ProtonMass       = 1.672621777e-27 * SI::Kilogram;
auto const ElementaryCharge = SI::ElectronVolt / SI::Volt;

Quantities::Dimensionless const FineStructureConstant = 7.2973525698e-3;

Quantities::Acceleration const StandardGravity = 9.80665 * SI::Metre /
                                                 SI::Second.Pow<2>();
}
}
