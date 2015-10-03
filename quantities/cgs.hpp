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

Length const Centimetre = si::Centi(si::Metre);
using si::Gram;
using si::Second;

Energy       const Erg  = 1e-7 * si::Joule;
Force        const Dyne = 1e-5 * si::Newton;
Acceleration const Gal  = Centimetre / Pow<2>(Second);

Pressure const Barye = 1 * Dyne / Pow<2>(Centimetre);

DynamicViscosity   const Poise  = Barye * Second;
KinematicViscosity const Stokes = Pow<2>(Centimetre) / Second;

Luminance   const Stilb = si::Candela / Pow<2>(Centimetre);
Illuminance const Phot  = Stilb * si::Steradian;

MagneticFluxDensity const Gauss   = 1e-4 * si::Tesla;
MagneticFlux        const Maxwell = Gauss * Pow<2>(Centimetre);
MagneticField       const Œrsted  =
    1e3 / (4 * π * si::Steradian) * si::Ampere / si::Metre;

SpectroscopicWavenumber const Kayser = si::Cycle / Centimetre;

}  // namespace cgs
}  // namespace quantities
}  // namespace principia
