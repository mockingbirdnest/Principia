#pragma once

#include "Quantities.hpp"
#include "NamedQuantities.hpp"
#include "Numbers.hpp"
#include "SI.hpp"

namespace Principia {
// This namespace contains the non-SI units associated with the CGS and the CGS-
// Gaussian system of units listed in the BIPM's SI brochure 8, section 4.1,
// table 9, http://www.bipm.org/en/si/si_brochure/chapter4/table9.html.
namespace CGS {
Quantities::Length const Centimetre = SI::Centi(SI::Metre);
using SI::Gram;
using SI::Second;

Quantities::Energy const Erg  = 1e-7 * SI::Joule;
Quantities::Force  const Dyne = 1e-5 * SI::Newton;

Quantities::Pressure const Barye = 1 * Dyne / Centimetre.Pow<2>();

Quantities::DynamicViscosity const Poise = Barye * Second;
Quantities::KinematicViscosity const Stokes = Centimetre.Pow<2>() / Second;

Quantities::Luminance   const Stilb = SI::Candela * Centimetre.Pow<-2>();
Quantities::Illuminance const Phot  = Stilb * SI::Steradian ;

Quantities::Acceleration const Gal = Centimetre / Second.Pow<2>();

Quantities::MagneticFluxDensity const Gauss   = 1e-4 * SI::Tesla;
Quantities::MagneticFlux        const Maxwell = Gauss * Centimetre.Pow<2>();
Quantities::MagneticField       const Œrsted  = 1e3 / (4 * π * SI::Steradian) *
                                                SI::Ampere / SI::Metre;

Quantities::SpectroscopicWavenumber const Kayser = SI::Cycle / Centimetre;
}  // namespace CGS
}  // namespace Principia
