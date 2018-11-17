
#pragma once

#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace quantities {
namespace constants {

// Resolution A of the 26th meeting of the CGPM.
constexpr Speed SpeedOfLight = 299'792'458 * (si::Metre / si::Second);
constexpr Action PlanckConstant = 6.626'070'15e-34 * si::Joule * si::Second;
constexpr Charge ElementaryCharge = 1.602'176'634e-19 * si::Coulomb;
constexpr Entropy BoltzmannConstant =  1.380'649e-23 * (si::Joule / si::Kelvin);
constexpr Inverse<Amount> AvogadroConstant = 6.022'140'76e23 * (1 / si::Mole);

// We use the 2014 CODATA recommended values.  We do not support uncertainties.
constexpr double FineStructureConstant = 7.297'352'5664e-3;
constexpr Quotient<GravitationalParameter, Mass> GravitationalConstant =
    6.674'08e-11 * si::Newton * Pow<2>(si::Metre) / Pow<2>(si::Kilogram);

constexpr Mass ElectronMass = 9.109'382'91e-31 * si::Kilogram;
constexpr Mass ProtonMass   = 1.672'621'777e-27 * si::Kilogram;

constexpr Permeability VacuumPermeability =
    2 * FineStructureConstant * PlanckConstant * si::Steradian /
    (Pow<2>(ElementaryCharge) * SpeedOfLight);
constexpr Permittivity VacuumPermittivity =
    Pow<2>(ElementaryCharge) /
    (2 * FineStructureConstant * PlanckConstant * SpeedOfLight * si::Steradian);

constexpr Acceleration StandardGravity =
    9.80665 * si::Metre / Pow<2>(si::Second);

}  // namespace constants
}  // namespace quantities
}  // namespace principia
