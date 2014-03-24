// Constants.h

#pragma once

#include "PhysicalQuantities.h"
#include "NamedQuantities.h"
#include "SIUnits.h"

namespace PhysicalQuantities {
Speed const SpeedOfLight = 299792458 * (Metre / Second);
Permeability const MagneticConstant = 4e-7 * π * (Henry / Metre);
Permittivity const ElectricConstant = 
  Dimensionless(1) / (MagneticConstant * SpeedOfLight * SpeedOfLight);
// We use the 2010 CODATA recommended values. We do not support uncertainties.
Action const ReducedPlanckConstant = 1.054571726e-34 * (Joule * Second);
auto const GravitationalConstant = 
  6.67384e-11 * (Newton * Metre * Metre / (Kilogram * Kilogram));

Entropy         const BoltzmannConstant = 1.3806488e-23 * (Joule / Kelvin);
Inverse<Amount> const AvogadroConstant  = 6.02214129 * (Uno / Mole);

Mass   const ElectronMass     = 9.10938291e-31 * Kilogram;
Mass   const ProtonMass       = 1.672621777e-27 * Kilogram;
Charge const ElementaryCharge = 1.602176565e-19 * Coulomb;

DimensionlessScalar const FineStructureConstant = 
  Dimensionless(7.2973525698e-3);
}