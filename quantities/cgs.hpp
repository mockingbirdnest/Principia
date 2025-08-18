#pragma once

#include "numerics/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace quantities {

// This namespace contains the non-SI units associated with the CGS and the CGS-
// Gaussian system of units listed in the BIPM's SI brochure 8, section 4.1,
// table 9, http://www.bipm.org/en/si/si_brochure/chapter4/table9.html.
namespace _cgs {
namespace internal {

using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

constexpr Length Centimetre = Centi(Metre);

constexpr Energy       Erg  = 1e-7 * Joule;
constexpr Force        Dyne = 1e-5 * Newton;
constexpr Acceleration Gal  = Centimetre / Pow<2>(Second);

constexpr Pressure Barye = 1 * Dyne / Pow<2>(Centimetre);

constexpr DynamicViscosity   Poise  = Barye * Second;
constexpr KinematicViscosity Stokes = Pow<2>(Centimetre) / Second;

constexpr Luminance   Stilb = Candela / Pow<2>(Centimetre);
constexpr Illuminance Phot  = Stilb * Steradian;

constexpr MagneticFluxDensity Gauss   = 1e-4 * Tesla;
constexpr MagneticFlux        Maxwell = Gauss * Pow<2>(Centimetre);
constexpr MagneticField       Œrsted  =
    1e3 / (4 * π * Steradian) * Ampere / Metre;

constexpr SpectroscopicWavenumber Kayser = 1 / Centimetre;

}  // namespace internal

using internal::Barye;
using internal::Centimetre;
using internal::Dyne;
using internal::Erg;
using internal::Gal;
using internal::Gauss;
using internal::Gram;
using internal::Kayser;
using internal::Maxwell;
using internal::Œrsted;
using internal::Phot;
using internal::Poise;
using internal::Second;
using internal::Stilb;
using internal::Stokes;

}  // namespace _cgs
}  // namespace quantities
}  // namespace principia
