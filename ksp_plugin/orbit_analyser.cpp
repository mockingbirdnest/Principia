
#include "ksp_plugin/orbit_analyser.hpp"

#include <algorithm>
#include <vector>

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
                             Ephemeris<Barycentric>::FixedStepParameters const&
                                 analysed_trajectory_parameters)
    : ephemeris_(ephemeris),
      analysed_trajectory_parameters_(analysed_trajectory_parameters) {}

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
  if (!analyser_.joinable()) {
    analyser_ = std::thread([this] { RepeatedlyAnalyseOrbit(); });
  }
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

OrbitAnalyser::Analysis* OrbitAnalyser::analysis() {
  return analysis_.has_value() ? &*analysis_ : nullptr;
}

double OrbitAnalyser::progress_of_next_analysis() const {
  return progress_of_next_analysis_;
}

void OrbitAnalyser::RepeatedlyAnalyseOrbit() {
  for (;;) {
    // No point in going faster than 50 Hz.
    std::chrono::steady_clock::time_point const wakeup_time =
        std::chrono::steady_clock::now() + std::chrono::milliseconds(20);

    if (!keep_analysing_) {
      return;
    }

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
                      parameters->primary};
    DiscreteTrajectory<Barycentric> trajectory;
    trajectory.Append(parameters->first_time,
                      parameters->first_degrees_of_freedom);
    std::vector<not_null<DiscreteTrajectory<Barycentric>*>> trajectories = {
        &trajectory};
    auto instance = ephemeris_->NewInstance(
        trajectories,
        Ephemeris<Barycentric>::NoIntrinsicAccelerations,
        analysed_trajectory_parameters_);
    for (Instant t =
             parameters->first_time + parameters->mission_duration / 0x1p10;
         trajectory.back().time <
         parameters->first_time + parameters->mission_duration;
         t += parameters->mission_duration / 0x1p10) {
      if (!ephemeris_->FlowWithFixedStep(t, *instance).ok()) {
        break;
      }
      progress_of_next_analysis_ =
          (trajectory.back().time - parameters->first_time) /
          parameters->mission_duration;
      if (!keep_analysing_) {
        return;
      }
    }
    analysis.mission_duration_ =
        trajectory.back().time - parameters->first_time;

    // TODO(egg): |next_analysis_percentage_| only reflects the progress of the
    // integration, but the analysis itself can take a while; this results in
    // the progress bar being stuck at 100% while the elements and nodes are
    // being computed.

    using PrimaryCentred = Frame<enum class PrimaryCentredTag>;
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
      analysis.elements_ = elements.ValueOrDie();
      // TODO(egg): max_abs_Cᴛₒ should probably depend on the number of
      // revolutions.
      analysis.closest_recurrence_ = OrbitRecurrence::ClosestRecurrence(
          analysis.elements_->nodal_period(),
          analysis.elements_->nodal_precession(),
          *parameters->primary,
          /*max_abs_Cᴛₒ=*/100);
      analysis.ground_track_ =
          OrbitGroundTrack::ForTrajectory(primary_centred_trajectory,
                                          *parameters->primary,
                                          /*mean_sun=*/std::nullopt);
      analysis.ResetRecurrence();
    }

    {
      absl::MutexLock l(&lock_);
      next_analysis_ = std::move(analysis);
    }

    std::this_thread::sleep_until(wakeup_time);
  }
}

Instant const& OrbitAnalyser::Analysis::first_time() const {
  return first_time_;
}

Time const& OrbitAnalyser::Analysis::mission_duration() const {
  return mission_duration_;
}

RotatingBody<Barycentric> const& OrbitAnalyser::Analysis::primary() const {
  return *primary_;
}

std::optional<OrbitalElements> const& OrbitAnalyser::Analysis::elements()
    const {
  return elements_;
}

std::optional<OrbitRecurrence> const& OrbitAnalyser::Analysis::recurrence()
    const {
  return recurrence_;
}

std::optional<OrbitGroundTrack> const& OrbitAnalyser::Analysis::ground_track()
    const {
  return ground_track_;
}

std::optional<OrbitGroundTrack::EquatorCrossingLongitudes> const&
OrbitAnalyser::Analysis::equatorial_crossings() const {
  return equatorial_crossings_;
}

void OrbitAnalyser::Analysis::SetRecurrence(
    OrbitRecurrence const& recurrence) {
  if (recurrence_ != recurrence) {
    recurrence_ = recurrence;
    if (ground_track_.has_value()) {
      equatorial_crossings_ = ground_track_->equator_crossing_longitudes(
          recurrence, /*first_ascending_pass_index=*/1);
    }
  }
}

void OrbitAnalyser::Analysis::ResetRecurrence() {
  if (closest_recurrence_.has_value()) {
    SetRecurrence(*closest_recurrence_);
  } else {
    recurrence_.reset();
    equatorial_crossings_.reset();
  }
}

OrbitAnalyser::Analysis::Analysis(
    Instant const& first_time,
    not_null<RotatingBody<Barycentric> const*> const primary)
    : first_time_(first_time), primary_(primary) {}

}  // namespace internal_orbit_analyser
}  // namespace ksp_plugin
}  // namespace principia
