#pragma once

#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace quantities {
namespace constants {

Speed const SpeedOfLight = 299792458 * (si::Metre / si::Second);
Permeability const VacuumPermeability =
    4e-7 * π * si::Steradian * si::Henry / si::Metre;
Permittivity const VacuumPermittivity =
    1 / (VacuumPermeability * Pow<2>(SpeedOfLight));
// We use the 2010 CODATA recommended values.  We do not support uncertainties.
AngularMomentum const ReducedPlanckConstant =
    1.054571726e-34 * si::Joule * si::Second / si::Radian;
Quotient<GravitationalParameter, Mass> const GravitationalConstant =
    6.67384e-11 * si::Newton * Pow<2>(si::Metre) / Pow<2>(si::Kilogram);
Entropy const BoltzmannConstant = 1.3806488e-23 * (si::Joule / si::Kelvin);
Amount::Inverse const AvogadroConstant = 6.02214129 * (1 / si::Mole);

Mass   const ElectronMass     = 9.10938291e-31 * si::Kilogram;
Mass   const ProtonMass       = 1.672621777e-27 * si::Kilogram;
Charge const ElementaryCharge = si::ElectronVolt / si::Volt;

double const FineStructureConstant = 7.2973525698e-3;

Acceleration const StandardGravity = 9.80665 * si::Metre / Pow<2>(si::Second);

}  // namespace constants
}  // namespace quantities
}  // namespace principia
