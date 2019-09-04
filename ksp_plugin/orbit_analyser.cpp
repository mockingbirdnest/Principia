#include "ksp_plugin/orbit_analyser.hpp"

#include "physics/discrete_trajectory.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_orbit_analyser {

using physics::DiscreteTrajectory;

OrbitAnalyser::OrbitAnalyser(not_null<Ephemeris<Barycentric>*> const ephemeris,
                             Ephemeris<Barycentric>::FixedStepParameters const
                                 analysed_trajectory_parameters)
    : ephemeris_(ephemeris),
      analysed_trajectory_parameters_(analysed_trajectory_parameters) {}

void OrbitAnalyser::RepeatedlyAnalyseOrbit() {
  for (;;) {
    // No point in going faster than 50 Hz.
    std::chrono::steady_clock::time_point const wakeup_time =
        std::chrono::steady_clock::now() + std::chrono::milliseconds(20);

    std::optional<Parameters> parameters;
    {
      absl::ReaderMutexLock l(&lock_);
      if (!parameters_.has_value()) {
        // No parameters, let's wait for them to appear.
        continue;
      }
      std::swap(parameters, parameters_);
    }

    if (parameters_->shutdown) {
      break;
    }

    DiscreteTrajectory<Barycentric> trajectory;
    std::vector<not_null<DiscreteTrajectory<Barycentric>*>> trajectories = {
        trajectory};
    auto instance = Ephemeris<Barycentric>::NewInstance(
        trajectories,
        Ephemeris<Barycentric>::NoIntrinsicAccelerations,
        analysed_trajectory_parameters_);
    for (trajectory) {
      ephemeris_->FlowWithFixedStep(t, instance);
    }
    {
      absl::MutexLock l(&prognosticator_lock_);
      SwapPrognostication(prognostication, status);
    }

    std::this_thread::sleep_until(wakeup_time);
  }
}

}  // namespace internal_orbit_analyser
}  // namespace ksp_plugin
}  // namespace principia
