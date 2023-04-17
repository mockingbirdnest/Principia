#include "ksp_plugin/orbit_analyser.hpp"

#include <algorithm>
#include <utility>
#include <vector>

#include "ksp_plugin/integrators.hpp"
#include "physics/body_centred_non_rotating_reference_frame.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/kepler_orbit.hpp"
#include "quantities/astronomy.hpp"

namespace principia {
namespace ksp_plugin {
namespace _orbit_analyser {
namespace internal {

using namespace principia::base::_jthread;
using namespace principia::base::_not_null;
using namespace principia::geometry::_frame;
using namespace principia::ksp_plugin::_integrators;
using namespace principia::physics::_body_centred_non_rotating_reference_frame;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_kepler_orbit;
using namespace principia::physics::_massive_body;
using namespace principia::physics::_massless_body;
using namespace principia::quantities::_astronomy;
using namespace principia::quantities::_quantities;

// TODO(egg): This could be implemented using ComputeApsides.
template<typename PrimaryCentred>
Interval<Length> RadialDistanceInterval(
    DiscreteTrajectory<PrimaryCentred> const& trajectory) {
  std::vector<Length> radial_distances;
  radial_distances.reserve(trajectory.size());
  DegreesOfFreedom<PrimaryCentred> const primary_dof{PrimaryCentred::origin,
                                                     PrimaryCentred::unmoving};
  for (auto const& [time, degrees_of_freedom] : trajectory) {
    radial_distances.push_back(
        (degrees_of_freedom.position() - primary_dof.position()).Norm());
  }

  Interval<Length> radial_distance_interval;
  for (auto const& r : radial_distances) {
    radial_distance_interval.Include(r);
  }
  return radial_distance_interval;
}

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
        [this, parameters]() {
          AnalyseOrbit(parameters).IgnoreError();
        });
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

absl::Status OrbitAnalyser::AnalyseOrbit(Parameters const& parameters) {
  Analysis analysis{parameters.first_time};

  RotatingBody<Barycentric> const* primary = nullptr;
  auto smallest_osculating_period = Infinity<Time>;
  auto const primary_status =
      FindBodyWithSmallestOsculatingPeriod(parameters,
                                           primary,
                                           smallest_osculating_period);
  RETURN_IF_ERROR(primary_status);

  if (primary != nullptr) {
    DiscreteTrajectory<Barycentric> trajectory;
    trajectory.segments().front().SetDownsampling(
        OrbitAnalyserDownsamplingParameters());

    Time const analysis_duration = std::min(
        parameters.extended_mission_duration.value_or(
            parameters.mission_duration),
        std::max(2 * smallest_osculating_period, parameters.mission_duration));
    RETURN_IF_ERROR(
        FlowWithProgressBar(parameters, analysis_duration, trajectory));
    analysis.mission_duration_ = trajectory.back().time - parameters.first_time;

    // TODO(egg): |next_analysis_percentage_| only reflects the progress of
    // the integration, but the analysis itself can take a while; this results
    // in the progress bar being stuck at 100% while the elements and nodes
    // are being computed.

    BodyCentredNonRotatingReferenceFrame<Barycentric, PrimaryCentred> const
        primary_centred(ephemeris_, primary);
    auto const status_or_primary_centred_trajectory =
        ToPrimaryCentred(primary_centred, trajectory);
    RETURN_IF_ERROR(status_or_primary_centred_trajectory);
    auto const& primary_centred_trajectory =
        status_or_primary_centred_trajectory.value();

    analysis.primary_ = primary;
    analysis.radial_distance_interval_ =
        RadialDistanceInterval(primary_centred_trajectory);
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
      if (analysis.closest_recurrence_->number_of_revolutions() == 0) {
        analysis.closest_recurrence_.reset();
      }

      std::optional<OrbitGroundTrack::MeanSun> mean_sun;
      RETURN_IF_ERROR(ComputeMeanSunIfPossible(parameters,
                                               primary_centred,
                                               mean_sun));

      auto ground_track =
          OrbitGroundTrack::ForTrajectory(primary_centred_trajectory,
                                          *primary,
                                          mean_sun);
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

absl::Status OrbitAnalyser::FindBodyWithSmallestOsculatingPeriod(
    Parameters const& parameters,
    RotatingBody<Barycentric> const*& primary,
    Time& smallest_osculating_period) {
  primary = nullptr;
  smallest_osculating_period = Infinity<Time>;
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

  return absl::OkStatus();
}

absl::Status OrbitAnalyser::FlowWithProgressBar(
    Parameters const& parameters,
    Time const& analysis_duration,
    DiscreteTrajectory<Barycentric>& trajectory) {
  trajectory.Append(parameters.first_time,
                    parameters.first_degrees_of_freedom).IgnoreError();

  std::vector<not_null<DiscreteTrajectory<Barycentric>*>> trajectories = {
      &trajectory};
  auto instance = ephemeris_->StoppableNewInstance(
      trajectories,
      Ephemeris<Barycentric>::NoIntrinsicAccelerations,
      analysed_trajectory_parameters_);
  RETURN_IF_STOPPED;

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
  return absl::OkStatus();
}

absl::Status
OrbitAnalyser::ComputeMeanSunIfPossible(
    Parameters const& parameters,
    BodyCentredNonRotatingReferenceFrame<Barycentric, PrimaryCentred> const&
        primary_centred,
    std::optional<OrbitGroundTrack::MeanSun>& mean_sun) {
  mean_sun = std::nullopt;
  auto const primary = primary_centred.centre();
  MassiveBody const* sun = nullptr;

  // Find the sun.  If there is none, we are done.
  for (auto const body : ephemeris_->bodies()) {
    if (body->name() == "Sun") {
      sun = body;
      break;
    }
  }

  if (primary != sun && sun != nullptr) {
    auto const sun_osculating_elements =
        KeplerOrbit<Barycentric>{
            *primary,
            *sun,
            ephemeris_->trajectory(sun)->EvaluateDegreesOfFreedom(
                parameters.first_time) -
                ephemeris_->trajectory(primary)->EvaluateDegreesOfFreedom(
                    parameters.first_time),
            parameters.first_time}
            .elements_at_epoch();
    Time const ephemeris_span = ephemeris_->t_max() - ephemeris_->t_min();
    if (ephemeris_span < 1.5 * *sun_osculating_elements.period &&
        ephemeris_span < 20 * JulianYear) {
      RETURN_IF_ERROR(
          ephemeris_->Prolong(ephemeris_->t_max() + 0.5 * JulianYear));
    }
    auto const sun_elements = OrbitalElements::ForTrajectory(
        *ephemeris_->trajectory(sun), primary_centred, *primary, *sun);
    if (sun_elements.ok()) {
      auto const& sun_mean_elements = sun_elements->mean_elements().front();
      mean_sun = OrbitGroundTrack::MeanSun{
          .epoch = sun_mean_elements.time,
          .mean_longitude_at_epoch =
              sun_mean_elements.longitude_of_ascending_node +
              sun_mean_elements.argument_of_periapsis +
              sun_mean_elements.mean_anomaly,
          .year = sun_elements->sidereal_period()};
    }
  }

  return absl::OkStatus();
}

absl::StatusOr<DiscreteTrajectory<OrbitAnalyser::PrimaryCentred>>
OrbitAnalyser::ToPrimaryCentred(
    BodyCentredNonRotatingReferenceFrame<Barycentric, PrimaryCentred> const&
        primary_centred,
    DiscreteTrajectory<Barycentric> const& trajectory) {
  DiscreteTrajectory<PrimaryCentred> primary_centred_trajectory;
  for (auto const& [time, degrees_of_freedom] : trajectory) {
    RETURN_IF_STOPPED;
    primary_centred_trajectory
        .Append(time,
                primary_centred.ToThisFrameAtTime(time)(degrees_of_freedom))
        .IgnoreError();
  }
  return primary_centred_trajectory;
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

std::optional<Interval<Length>>
OrbitAnalyser::Analysis::radial_distance_interval() const {
  return radial_distance_interval_;
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

}  // namespace internal
}  // namespace _orbit_analyser
}  // namespace ksp_plugin
}  // namespace principia
