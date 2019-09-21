
#include "ksp_plugin/orbit_analyser.hpp"

#include "physics/body_centred_non_rotating_dynamic_frame.hpp"
#include "physics/discrete_trajectory.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_orbit_analyser {

using geometry::Frame;
using physics::DiscreteTrajectory;
using physics::MasslessBody;
using physics::BodyCentredNonRotatingDynamicFrame;
using quantities::IsFinite;

OrbitAnalyser::OrbitAnalyser(not_null<Ephemeris<Barycentric>*> const ephemeris,
                             Ephemeris<Barycentric>::FixedStepParameters const
                                 analysed_trajectory_parameters)
    : ephemeris_(ephemeris),
      analysed_trajectory_parameters_(analysed_trajectory_parameters),
      analyser_([this] { RepeatedlyAnalyseOrbit(); }) {}

OrbitAnalyser::~OrbitAnalyser() {
  if (analyser_.joinable()) {
    keep_analysing_ = false;
    analyser_.join();
  }
}

void OrbitAnalyser::RequestAnalysis(
    Instant const& first_time,
    DegreesOfFreedom<Barycentric> const& first_degrees_of_freedom,
    Time const& mission_duration,
    not_null<RotatingBody<Barycentric> const*> primary) {
  Ephemeris<Barycentric>::Guard guard(ephemeris_);
  if (ephemeris_->t_min() > first_time) {
    // Too much has been forgotten; we cannot perform this analysis.
    return;
  }
  absl::MutexLock l(&lock_);
  parameters_ = {std::move(guard),
                 first_time,
                 first_degrees_of_freedom,
                 mission_duration,
                 primary};
}

void OrbitAnalyser::RefreshAnalysis() {
  absl::MutexLock l(&lock_);
  if (next_analysis_.has_value()) {
    analysis_ = next_analysis_;
    next_analysis_.reset();
  }
}

std::optional<OrbitAnalyser::Analysis>& OrbitAnalyser::analysis() {
  return analysis_;
}

int8_t OrbitAnalyser::next_analysis_percentage() const {
  return next_analysis_percentage_;
}

void OrbitAnalyser::RepeatedlyAnalyseOrbit() {
  while (keep_analysing_) {
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

    Analysis analysis{parameters->first_time,
                      parameters->mission_duration,
                      parameters->primary};
    DiscreteTrajectory<Barycentric> trajectory;
    trajectory.Append(parameters->first_time, parameters->first_degrees_of_freedom);
    std::vector<not_null<DiscreteTrajectory<Barycentric>*>> trajectories = {
        &trajectory};
    auto instance = ephemeris_->NewInstance(
        trajectories,
        Ephemeris<Barycentric>::NoIntrinsicAccelerations,
        analysed_trajectory_parameters_);
    for (Instant t =
             parameters->first_time + parameters->mission_duration / 100;
         trajectory.back().time <
         parameters->first_time + parameters->mission_duration;
         t += parameters->mission_duration / 100) {
      if (!ephemeris_->FlowWithFixedStep(t, *instance).ok()) {
        break;
      }
      next_analysis_percentage_ =
          100 * (trajectory.back().time - parameters->first_time) /
          parameters->mission_duration;
      if (!keep_analysing_) {
        return;
      }
    }

    enum class PrimaryCentredTag { tag };
    using PrimaryCentred = Frame<PrimaryCentredTag,
                                 PrimaryCentredTag::tag,
                                 /*frame_is_inertial=*/false>;
    BodyCentredNonRotatingDynamicFrame<Barycentric, PrimaryCentred>
        primary_centred(ephemeris_, parameters->primary);
    DiscreteTrajectory<PrimaryCentred> primary_centred_trajectory;
    for (auto const& [time, degrees_of_freedom] : trajectory) {
      primary_centred_trajectory.Append(
          time, primary_centred.ToThisFrameAtTime(time)(degrees_of_freedom));
    }

    auto const elements = OrbitalElements::ForTrajectory(
        primary_centred_trajectory, *parameters->primary, MasslessBody{});
    if (elements.ok()) {
      analysis.elements = elements.ValueOrDie();
      if (IsFinite(analysis.elements->nodal_period()) &&
          IsFinite(analysis.elements->nodal_precession())) {
        analysis.auto_detected_recurrence = OrbitRecurrence::ClosestRecurrence(
            analysis.elements->nodal_period(),
            analysis.elements->nodal_precession(),
            *parameters->primary,
            /*max_abs_Cᴛₒ=*/100);
        analysis.ground_track =
            OrbitGroundTrack::ForTrajectory(primary_centred_trajectory,
                                            *parameters->primary,
                                            /*mean_sun=*/std::nullopt);
      }
    }

    {
      absl::MutexLock l(&lock_);
      next_analysis_ = std::move(analysis);
    }

    std::this_thread::sleep_until(wakeup_time);
  }
}

}  // namespace internal_orbit_analyser
}  // namespace ksp_plugin
}  // namespace principia
