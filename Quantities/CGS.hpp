#pragma once

#include "Quantities.hpp"
#include "NamedQuantities.hpp"
#include "Numbers.hpp"
#include "SI.hpp"

namespace principia {
// This namespace contains the non-SI units associated with the CGS and the CGS-
// Gaussian system of units listed in the BIPM's SI brochure 8, section 4.1,
// table 9, http://www.bipm.org/en/si/si_brochure/chapter4/table9.html.
namespace CGS {
quantities::Length const Centimetre = SI::Centi(SI::Metre);
using SI::Gram;
using SI::Second;

quantities::Energy const Erg  = 1e-7 * SI::Joule;
quantities::Force  const Dyne = 1e-5 * SI::Newton;

quantities::Pressure const Barye = 1 * Dyne / Centimetre.Pow<2>();

quantities::DynamicViscosity const Poise = Barye * Second;
quantities::KinematicViscosity const Stokes = Centimetre.Pow<2>() / Second;

quantities::Luminance   const Stilb = SI::Candela * Centimetre.Pow<-2>();
quantities::Illuminance const Phot  = Stilb * SI::Steradian ;

quantities::Acceleration const Gal = Centimetre / Second.Pow<2>();

quantities::MagneticFluxDensity const Gauss   = 1e-4 * SI::Tesla;
quantities::MagneticFlux        const Maxwell = Gauss * Centimetre.Pow<2>();
quantities::MagneticField       const Œrsted  = 1e3 / (4 * π * SI::Steradian) *
                                                SI::Ampere / SI::Metre;

quantities::SpectroscopicWavenumber const Kayser = SI::Cycle / Centimetre;
}  // namespace CGS
}  // namespace principia
