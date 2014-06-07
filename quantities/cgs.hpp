#pragma once

#include "quantities/named_quantities.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
// This namespace contains the non-SI units associated with the CGS and the CGS-
// Gaussian system of units listed in the BIPM's SI brochure 8, section 4.1,
// table 9, http://www.bipm.org/en/si/si_brochure/chapter4/table9.html.
namespace cgs {
quantities::Length const Centimetre = si::Centi(si::Metre);
using si::Gram;
using si::Second;

quantities::Energy const Erg  = 1e-7 * si::Joule;
quantities::Force  const Dyne = 1e-5 * si::Newton;

quantities::Pressure const Barye = 1 * Dyne / quantities::Pow<2>(Centimetre);

quantities::DynamicViscosity const Poise = Barye * Second;
quantities::KinematicViscosity const Stokes =
    quantities::Pow<2>(Centimetre) / Second;

quantities::Luminance const Stilb =
    si::Candela / quantities::Pow<2>(Centimetre);
quantities::Illuminance const Phot = Stilb * si::Steradian;

quantities::Acceleration const Gal = Centimetre / quantities::Pow<2>(Second);

quantities::MagneticFluxDensity const Gauss = 1e-4 * si::Tesla;
quantities::MagneticFlux const Maxwell = Gauss * quantities::Pow<2>(Centimetre);
quantities::MagneticField const Œrsted =
    1e3 / (4 * π * si::Steradian) * si::Ampere / si::Metre;

quantities::SpectroscopicWavenumber const Kayser = si::Cycle / Centimetre;
}  // namespace cgs
}  // namespace principia
