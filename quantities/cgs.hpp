#pragma once

#include "quantities/named_quantities.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace quantities {

// This namespace contains the non-SI units associated with the CGS and the CGS-
// Gaussian system of units listed in the BIPM's SI brochure 8, section 4.1,
// table 9, http://www.bipm.org/en/si/si_brochure/chapter4/table9.html.
namespace cgs {

Length constexpr Centimetre = si::Centi(si::Metre);
using si::Gram;
using si::Second;

Energy       constexpr Erg  = 1e-7 * si::Joule;
Force        constexpr Dyne = 1e-5 * si::Newton;
Acceleration constexpr Gal  = Centimetre / Pow<2>(Second);

Pressure constexpr Barye = 1 * Dyne / Pow<2>(Centimetre);

DynamicViscosity   constexpr Poise  = Barye * Second;
KinematicViscosity constexpr Stokes = Pow<2>(Centimetre) / Second;

Luminance   constexpr Stilb = si::Candela / Pow<2>(Centimetre);
Illuminance constexpr Phot  = Stilb * si::Steradian;

MagneticFluxDensity constexpr Gauss   = 1e-4 * si::Tesla;
MagneticFlux        constexpr Maxwell = Gauss * Pow<2>(Centimetre);
MagneticField       constexpr Œrsted  =
    1e3 / (4 * π * si::Steradian) * si::Ampere / si::Metre;

SpectroscopicWavenumber constexpr Kayser = si::Cycle / Centimetre;

}  // namespace cgs
}  // namespace quantities
}  // namespace principia
