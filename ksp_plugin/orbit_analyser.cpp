
#include "ksp_plugin/orbit_analyser.hpp"

#include <algorithm>
#include <utility>
#include <vector>

#include "ksp_plugin/integrators.hpp"
#include "physics/body_centred_non_rotating_dynamic_frame.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/kepler_orbit.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_orbit_analyser {

using base::dynamic_cast_not_null;
using base::MakeStoppableThread;
using geometry::Frame;
using geometry::NonRotating;
using physics::BodyCentredNonRotatingDynamicFrame;
using physics::DiscreteTrajectory;
using physics::KeplerOrbit;
using physics::MasslessBody;
using quantities::IsFinite;
using quantities::Infinity;

OrbitAnalyser::OrbitAnalyser(
    not_null<Ephemeris<Barycentric>*> const ephemeris,
    Ephemeris<Barycentric>::FixedStepParameters analysed_trajectory_parameters)
    : ephemeris_(ephemeris),
      analysed_trajectory_parameters_(
          std::move(analysed_trajectory_parameters)) {}

OrbitAnalyser::~OrbitAnalyser() {
  // Ensure that we do not have a thread still running with references to the
  // members of this class when those are destroyed.
  Interrupt();
}

void OrbitAnalyser::Interrupt() {
  analyser_ = jthread();
  // We are single-threaded here, no need to lock.
  analyser_idle_ = true;
}

void OrbitAnalyser::RequestAnalysis(Parameters const& parameters) {
  if (ephemeris_->t_min() > parameters.first_time) {
    // Too much has been forgotten; we cannot perform this analysis.
    return;
  }
  last_parameters_ = parameters;
  absl::MutexLock l(&lock_);
  // Only process this request if there is no analysis in progress.
  if (analyser_idle_) {
    analyser_idle_ = false;
    analyser_ = MakeStoppableThread(
        [this](Parameters parameters) { AnalyseOrbit(parameters); },
        parameters);
  }
}

std::optional<OrbitAnalyser::Parameters> const& OrbitAnalyser::last_parameters()
    const {
  return last_parameters_;
}

void OrbitAnalyser::RefreshAnalysis() {
  absl::MutexLock l(&lock_);
  if (next_analysis_.has_value()) {
    analysis_ = std::move(next_analysis_);
    next_analysis_.reset();
  }
}

OrbitAnalyser::Analysis* OrbitAnalyser::analysis() {
  return analysis_.has_value() ? &*analysis_ : nullptr;
}

double OrbitAnalyser::progress_of_next_analysis() const {
  return progress_of_next_analysis_;
}

absl::Status OrbitAnalyser::AnalyseOrbit(Parameters const parameters) {
  Analysis analysis{parameters.first_time};
  DiscreteTrajectory<Barycentric> trajectory;
  trajectory.SetDownsampling(DefaultDownsamplingParameters());
  trajectory.Append(parameters.first_time, parameters.first_degrees_of_freedom);

  RotatingBody<Barycentric> const* primary = nullptr;
  auto smallest_osculating_period = Infinity<Time>;
  for (auto const body : ephemeris_->bodies()) {
    RETURN_IF_STOPPED;
    auto const initial_osculating_elements =
        KeplerOrbit<Barycentric>{
            *body,
            MasslessBody{},
            parameters.first_degrees_of_freedom -
                ephemeris_->trajectory(body)->EvaluateDegreesOfFreedom(
                    parameters.first_time),
            parameters.first_time}
            .elements_at_epoch();
    if (initial_osculating_elements.period.has_value() &&
        initial_osculating_elements.period < smallest_osculating_period) {
      smallest_osculating_period = *initial_osculating_elements.period;
      primary = dynamic_cast_not_null<RotatingBody<Barycentric> const*>(body);
    }
  }
  if (primary != nullptr) {
    std::vector<not_null<DiscreteTrajectory<Barycentric>*>> trajectories = {
        &trajectory};
    auto instance = ephemeris_->StoppableNewInstance(
        trajectories,
        Ephemeris<Barycentric>::NoIntrinsicAccelerations,
        analysed_trajectory_parameters_);
    RETURN_IF_STOPPED;

    Time const analysis_duration = std::min(
        parameters.extended_mission_duration.value_or(
            parameters.mission_duration),
        std::max(2 * smallest_osculating_period, parameters.mission_duration));
    constexpr double progress_bar_steps = 0x1p10;
    for (double n = 0; n <= progress_bar_steps; ++n) {
      Instant const t =
          parameters.first_time + n / progress_bar_steps * analysis_duration;
      if (!ephemeris_->FlowWithFixedStep(t, *instance.value()).ok()) {
        // TODO(egg): Report that the integration failed.
        break;
      }
      progress_of_next_analysis_ =
          (trajectory.back().time - parameters.first_time) / analysis_duration;
      RETURN_IF_STOPPED;
    }
    analysis.mission_duration_ = trajectory.back().time - parameters.first_time;

    // TODO(egg): |next_analysis_percentage_| only reflects the progress of
    // the integration, but the analysis itself can take a while; this results
    // in the progress bar being stuck at 100% while the elements and nodes
    // are being computed.

    using PrimaryCentred = Frame<enum class PrimaryCentredTag, NonRotating>;
    DiscreteTrajectory<PrimaryCentred> primary_centred_trajectory;
    BodyCentredNonRotatingDynamicFrame<Barycentric, PrimaryCentred>
        body_centred(ephemeris_, primary);
    for (auto const& [time, degrees_of_freedom] : trajectory) {
      RETURN_IF_STOPPED;
      primary_centred_trajectory.Append(
          time, body_centred.ToThisFrameAtTime(time)(degrees_of_freedom));
    }
    analysis.primary_ = primary;
    auto elements = OrbitalElements::ForTrajectory(
        primary_centred_trajectory, *primary, MasslessBody{});
    // We do not RETURN_IF_ERROR as ForTrajectory can return non-CANCELLED
    // statuses.
    RETURN_IF_STOPPED;
    if (elements.ok()) {
      analysis.elements_ = std::move(elements).value();
      // TODO(egg): max_abs_Cᴛₒ should probably depend on the number of
      // revolutions.
      analysis.closest_recurrence_ = OrbitRecurrence::ClosestRecurrence(
          analysis.elements_->nodal_period(),
          analysis.elements_->nodal_precession(),
          *primary,
          /*max_abs_Cᴛₒ=*/100);
      auto ground_track =
          OrbitGroundTrack::ForTrajectory(primary_centred_trajectory,
                                          *primary,
                                          /*mean_sun=*/std::nullopt);
      RETURN_IF_ERROR(ground_track);
      analysis.ground_track_ = std::move(ground_track).value();
      analysis.ResetRecurrence();
    }
  }

  absl::MutexLock l(&lock_);
  next_analysis_ = std::move(analysis);
  analyser_idle_ = true;
  return absl::OkStatus();
}

Instant const& OrbitAnalyser::Analysis::first_time() const {
  return first_time_;
}

Time const& OrbitAnalyser::Analysis::mission_duration() const {
  return mission_duration_;
}

RotatingBody<Barycentric> const* OrbitAnalyser::Analysis::primary() const {
  return primary_;
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

OrbitAnalyser::Analysis::Analysis(Instant const& first_time)
    : first_time_(first_time) {}

}  // namespace internal_orbit_analyser
}  // namespace ksp_plugin
}  // namespace principia
