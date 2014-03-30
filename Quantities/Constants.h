// Constants.h

#pragma once

#include "Quantities.h"
#include "NamedQuantities.h"
#include "SI.h"

namespace Quantities {
Speed const SpeedOfLight = 299792458 * (Metre / Second);
Permeability const MagneticConstant = 4 * π * Steradian * 1e-7 * Henry / Metre;
Permittivity const ElectricConstant = 1 /
  (MagneticConstant * SpeedOfLight * SpeedOfLight);
// We use the 2010 CODATA recommended values. We do not support uncertainties.
Action const ReducedPlanckConstant = 1.054571726e-34 * (Joule * Second);
auto const GravitationalConstant = 
  6.67384e-11 * Newton * Metre * Metre / (Kilogram * Kilogram);

Entropy         const BoltzmannConstant = 1.3806488e-23 * (Joule / Kelvin);
Inverse<Amount> const AvogadroConstant  = 6.02214129 * (1 / Mole);

Mass   const ElectronMass     = 9.10938291e-31 * Kilogram;
Mass   const ProtonMass       = 1.672621777e-27 * Kilogram;
Charge const ElementaryCharge = 1.602176565e-19 * Coulomb;

Dimensionless const FineStructureConstant = 7.2973525698e-3;
}
