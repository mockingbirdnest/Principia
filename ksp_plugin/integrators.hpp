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

// Parameters for downsamplng after fixed-step integration.
constexpr std::int64_t MaxDenseIntervals = 10'000;
constexpr Length DownsamplingTolerance = 10 * Metre;

// Factories for parameters used to control integration.
Ephemeris<Barycentric>::AccuracyParameters
DefaultEphemerisAccuracyParameters();
Ephemeris<Barycentric>::FixedStepParameters
DefaultEphemerisFixedStepParameters();
Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters
DefaultBurnParameters();
Ephemeris<Barycentric>::FixedStepParameters DefaultHistoryParameters();
Ephemeris<Barycentric>::AdaptiveStepParameters DefaultPredictionParameters();
Ephemeris<Barycentric>::AdaptiveStepParameters DefaultPsychohistoryParameters();

}  // namespace internal_integrators

using internal_integrators::DefaultBurnParameters;
using internal_integrators::DefaultEphemerisAccuracyParameters;
using internal_integrators::DefaultEphemerisFixedStepParameters;
using internal_integrators::DefaultHistoryParameters;
using internal_integrators::DefaultPredictionParameters;
using internal_integrators::DefaultPsychohistoryParameters;
using internal_integrators::DownsamplingTolerance;
using internal_integrators::MaxDenseIntervals;

}  // namespace ksp_plugin
}  // namespace principia
