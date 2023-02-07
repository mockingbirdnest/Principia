#pragma once

#include "ksp_plugin/frames.hpp"
#include "physics/discrete_trajectory_segment.hpp"
#include "physics/ephemeris.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_integrators {

using physics::DiscreteTrajectorySegment;
using physics::Ephemeris;
using quantities::Length;
using quantities::si::Metre;
using quantities::si::Milli;

// Parameters for downsampling after fixed-step integration.
DiscreteTrajectorySegment<Barycentric>::DownsamplingParameters
DefaultDownsamplingParameters();

// Parameters for the orbit analyser.  Finer-grained than the default to obtain
// reasonable element.  The 10 times smaller tolerance results in trajectories
// that are about twice as big.
DiscreteTrajectorySegment<Barycentric>::DownsamplingParameters
OrbitAnalyserDownsamplingParameters();

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
using internal_integrators::DefaultDownsamplingParameters;
using internal_integrators::DefaultEphemerisAccuracyParameters;
using internal_integrators::DefaultEphemerisFixedStepParameters;
using internal_integrators::DefaultHistoryParameters;
using internal_integrators::DefaultPredictionParameters;
using internal_integrators::DefaultPsychohistoryParameters;
using internal_integrators::OrbitAnalyserDownsamplingParameters;

}  // namespace ksp_plugin
}  // namespace principia
