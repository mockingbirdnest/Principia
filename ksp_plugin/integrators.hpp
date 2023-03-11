#pragma once

#include "ksp_plugin/frames.hpp"
#include "physics/discrete_trajectory_segment.hpp"
#include "physics/ephemeris.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace ksp_plugin {
namespace _integrators {
namespace internal {

using namespace principia::physics::_discrete_trajectory_segment;
using namespace principia::physics::_ephemeris;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

// Parameters for downsampling after fixed-step integration.
DiscreteTrajectorySegment<Barycentric>::DownsamplingParameters
DefaultDownsamplingParameters();

// Parameters for the orbit analyser.  Finer-grained than the default to obtain
// reasonable element.  The 10'000 times smaller tolerance results in
// trajectories that are about 10 times as big.
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

}  // namespace internal

using internal::DefaultBurnParameters;
using internal::DefaultDownsamplingParameters;
using internal::DefaultEphemerisAccuracyParameters;
using internal::DefaultEphemerisFixedStepParameters;
using internal::DefaultHistoryParameters;
using internal::DefaultPredictionParameters;
using internal::DefaultPsychohistoryParameters;
using internal::OrbitAnalyserDownsamplingParameters;

}  // namespace _integrators
}  // namespace ksp_plugin
}  // namespace principia

namespace principia::ksp_plugin {
using namespace principia::ksp_plugin::_integrators;
}  // namespace principia::ksp_plugin
