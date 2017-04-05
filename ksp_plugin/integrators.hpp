#pragma once

#include "ksp_plugin/frames.hpp"
#include "physics/ephemeris.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_integrators {

using physics::Ephemeris;
using quantities::Length;
using quantities::si::Metre;
using quantities::si::Milli;

constexpr Length default_ephemeris_fitting_tolerance = 1 * Milli(Metre);

// Factories for parameters used to control integration.
Ephemeris<Barycentric>::FixedStepParameters DefaultEphemerisParameters();
Ephemeris<Barycentric>::FixedStepParameters DefaultHistoryParameters();
Ephemeris<Barycentric>::AdaptiveStepParameters DefaultPredictionParameters();
Ephemeris<Barycentric>::AdaptiveStepParameters DefaultProlongationParameters();

}  // namespace internal_integrators

using internal_integrators::default_ephemeris_fitting_tolerance;
using internal_integrators::DefaultEphemerisParameters;
using internal_integrators::DefaultHistoryParameters;
using internal_integrators::DefaultPredictionParameters;
using internal_integrators::DefaultProlongationParameters;

}  // namespace ksp_plugin
}  // namespace principia
