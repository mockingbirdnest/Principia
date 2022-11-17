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

constexpr Length Centimetre = si::Centi(si::Metre);
using si::Gram;
using si::Second;

constexpr Energy       Erg  = 1e-7 * si::Joule;
constexpr Force        Dyne = 1e-5 * si::Newton;
constexpr Acceleration Gal  = Centimetre / Pow<2>(Second);

constexpr Pressure Barye = 1 * Dyne / Pow<2>(Centimetre);

constexpr DynamicViscosity   Poise  = Barye * Second;
constexpr KinematicViscosity Stokes = Pow<2>(Centimetre) / Second;

constexpr Luminance   Stilb = si::Candela / Pow<2>(Centimetre);
constexpr Illuminance Phot  = Stilb * si::Steradian;

constexpr MagneticFluxDensity Gauss   = 1e-4 * si::Tesla;
constexpr MagneticFlux        Maxwell = Gauss * Pow<2>(Centimetre);
constexpr MagneticField       Œrsted  =
    1e3 / (4 * π * si::Steradian) * si::Ampere / si::Metre;

constexpr SpectroscopicWavenumber Kayser = 1 / Centimetre;

}  // namespace cgs
}  // namespace quantities
}  // namespace principia
