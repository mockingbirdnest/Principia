#include "ksp_plugin/flight_plan.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iterator>
#include <memory>
#include <thread>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "absl/strings/str_cat.h"
#include "base/algebra.hpp"
#include "base/status_utilities.hpp"  // 🧙 For CHECK_OK.
#include "glog/logging.h"
#include "integrators/embedded_explicit_generalized_runge_kutta_nyström_integrator.hpp"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/methods.hpp"
#include "ksp_plugin/integrators.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/make_not_null.hpp"

namespace principia {
namespace ksp_plugin {
namespace _flight_plan {
namespace internal {

using namespace principia::base::_algebra;
using namespace principia::integrators::_embedded_explicit_generalized_runge_kutta_nyström_integrator;  // NOLINT
using namespace principia::integrators::_embedded_explicit_runge_kutta_nyström_integrator;  // NOLINT
using namespace principia::integrators::_methods;
using namespace principia::ksp_plugin::_integrators;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_make_not_null;
using namespace std::chrono_literals;

namespace {

inline absl::Status BadDesiredFinalTime() {
  return absl::Status(FlightPlan::bad_desired_final_time,
                      "Bad desired final time");
}

inline absl::Status DoesNotFit() {
  return absl::Status(FlightPlan::does_not_fit, "Does not fit");
}

inline absl::Status Singular(Square<Speed> const& Δv²) {
  return absl::Status(FlightPlan::singular,
                      absl::StrCat("Singular: ", DebugString(Δv²)));
}

}  // namespace

FlightPlan::FlightPlan(
    Mass const& initial_mass,
    Instant const& initial_time,
    DegreesOfFreedom<Barycentric> initial_degrees_of_freedom,
    Instant const& desired_final_time,
    not_null<Ephemeris<Barycentric>*> const ephemeris,
    Ephemeris<Barycentric>::AdaptiveStepParameters adaptive_step_parameters,
    Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters
        generalized_adaptive_step_parameters)
    : initial_mass_(initial_mass),
      initial_time_(initial_time),
      initial_degrees_of_freedom_(std::move(initial_degrees_of_freedom)),
      ephemeris_(ephemeris),
      desired_final_time_(desired_final_time),
      adaptive_step_parameters_(std::move(adaptive_step_parameters)),
      generalized_adaptive_step_parameters_(
          std::move(generalized_adaptive_step_parameters)) {
  CHECK(desired_final_time_ >= initial_time_);
  MakeProlongator(desired_final_time_);

  // Set the first point of the first coasting trajectory.
  trajectory_.Append(initial_time_, initial_degrees_of_freedom_).IgnoreError();
  segments_.emplace_back(trajectory_.segments().begin());

  coast_analysers_.push_back(make_not_null_unique<OrbitAnalyser>(
      ephemeris_, DefaultHistoryParameters()));
  CHECK(manœuvres_.empty());
  ComputeSegments(manœuvres_.begin(),
                  manœuvres_.end(),
                  max_ephemeris_steps_per_frame).IgnoreError();
}

FlightPlan::FlightPlan(FlightPlan const& other)
    : initial_mass_(other.initial_mass_),
      initial_time_(other.initial_time_),
      initial_degrees_of_freedom_(other.initial_degrees_of_freedom_),
      ephemeris_(other.ephemeris_),
      desired_final_time_(other.desired_final_time_),
      anomalous_segments_(other.anomalous_segments_),
      manœuvres_(other.manœuvres_),
      analysis_is_enabled_(other.analysis_is_enabled_),
      adaptive_step_parameters_(other.adaptive_step_parameters_),
      generalized_adaptive_step_parameters_(
          other.generalized_adaptive_step_parameters_) {
  MakeProlongator(desired_final_time_);
  bool first_segment = true;
  for (auto const& other_segment : other.trajectory_.segments()) {
    if (!first_segment) {
      // The first segment was created by the constructor of the trajectory.
      trajectory_.NewSegment();
    }
    bool first_point = true;
    for (auto const& [time, degrees_of_freedom] : other_segment) {
      // For segments other than the first, `NewSegment` copied the last point
      // of the previous segment.
      if (!first_point || first_segment) {
        CHECK_OK(trajectory_.Append(time, degrees_of_freedom));
      }
      first_point = false;
    }
    first_segment = false;
  }
  for (auto it = trajectory_.segments().begin();
       it != trajectory_.segments().end();
       ++it) {
    segments_.push_back(it);
  }
  for (int i = 0; i < other.coast_analysers_.size(); ++i) {
    coast_analysers_.push_back(make_not_null_unique<OrbitAnalyser>(
        ephemeris_, DefaultHistoryParameters()));
  }
}

Instant FlightPlan::initial_time() const {
  return initial_time_;
}

Instant FlightPlan::actual_final_time() const {
  return segments_.back()->back().time;
}

Instant FlightPlan::desired_final_time() const {
  return desired_final_time_;
}

int FlightPlan::number_of_manœuvres() const {
  return manœuvres_.size();
}

int FlightPlan::number_of_anomalous_manœuvres() const {
  return (anomalous_segments_ - 1) / 2;
}

absl::Status const& FlightPlan::anomalous_status() const {
  return anomalous_status_;
}

NavigationManœuvre const& FlightPlan::GetManœuvre(int const index) const {
  CHECK_LE(0, index);
  CHECK_LT(index, number_of_manœuvres());
  return manœuvres_[index];
}

absl::Status FlightPlan::Insert(NavigationManœuvre::Burn const& burn,
                                int const index) {
  CHECK_GE(index, 0);
  CHECK_LE(index, number_of_manœuvres());
  NavigationManœuvre const manœuvre(
      index == 0 ? initial_mass_ : manœuvres_[index - 1].final_mass(),
      burn);
  if (manœuvre.IsSingular()) {
    return Singular(manœuvre.Δv().Norm²());
  }
  if (!manœuvre.FitsBetween(start_of_previous_coast(index),
                            start_of_burn(index))) {
    return DoesNotFit();
  }
  manœuvres_.insert(manœuvres_.begin() + index, manœuvre);
  coast_analysers_.insert(coast_analysers_.begin() + index + 1,
                          make_not_null_unique<OrbitAnalyser>(
                              ephemeris_, DefaultHistoryParameters()));
  UpdateInitialMassOfManœuvresAfter(index);
  PopSegmentsAffectedByManœuvre(index);
  return ComputeSegments(manœuvres_.begin() + index,
                         manœuvres_.end(),
                         max_ephemeris_steps_per_frame);
}

absl::Status FlightPlan::Remove(int const index) {
  CHECK_GE(index, 0);
  CHECK_LT(index, number_of_manœuvres());
  manœuvres_.erase(manœuvres_.begin() + index);
  coast_analysers_.erase(coast_analysers_.begin() + index + 1);
  UpdateInitialMassOfManœuvresAfter(index);
  PopSegmentsAffectedByManœuvre(index);
  return ComputeSegments(manœuvres_.begin() + index,
                         manœuvres_.end(),
                         max_ephemeris_steps_per_frame);
}

absl::Status FlightPlan::Replace(NavigationManœuvre::Burn const& burn,
                                 int const index) {
  CHECK_LE(0, index);
  CHECK_LT(index, number_of_manœuvres());
  NavigationManœuvre const manœuvre(manœuvres_[index].initial_mass(),
                                    burn);
  if (manœuvre.IsSingular()) {
    return Singular(manœuvre.Δv().Norm²());
  }
  if (index == number_of_manœuvres() - 1) {
    // This is the last manœuvre.  If it doesn't fit just because the flight
    // plan is too short, extend the flight plan.
    if (manœuvre.IsAfter(start_of_previous_coast(index))) {
      desired_final_time_ =
          std::max(desired_final_time_, manœuvre.final_time());
    } else {
      return DoesNotFit();
    }
  } else if (!manœuvre.FitsBetween(start_of_previous_coast(index),
                                   start_of_next_burn(index))) {
    return DoesNotFit();
  }

  // Replace the manœuvre at position `index` and rebuild all the ones that
  // follow as they may have a different initial mass.  Also pop the segments
  // that we'll recompute.
  manœuvres_[index] = manœuvre;
  UpdateInitialMassOfManœuvresAfter(index);

  // TODO(phl): Recompute as late as possible.
  PopSegmentsAffectedByManœuvre(index);
  return ComputeSegments(manœuvres_.begin() + index,
                         manœuvres_.end(),
                         max_ephemeris_steps_per_frame);
}

absl::Status FlightPlan::SetDesiredFinalTime(
    Instant const& desired_final_time) {
  if (desired_final_time < start_of_last_coast()) {
    return BadDesiredFinalTime();
  }

  desired_final_time_ = desired_final_time;
  MakeProlongator(desired_final_time_);

  // Reset the last coast and recompute it.
  ResetLastSegment();
  return ComputeSegments(manœuvres_.end(),
                         manœuvres_.end(),
                         max_ephemeris_steps_per_frame);
}

absl::Status FlightPlan::SetAdaptiveStepParameters(
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        adaptive_step_parameters,
    Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters const&
        generalized_adaptive_step_parameters) {
  adaptive_step_parameters_ = adaptive_step_parameters;
  generalized_adaptive_step_parameters_ = generalized_adaptive_step_parameters;
  return RecomputeAllSegments();
}

Ephemeris<Barycentric>::AdaptiveStepParameters const&
FlightPlan::adaptive_step_parameters() const {
  return adaptive_step_parameters_;
}

Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters const&
FlightPlan::generalized_adaptive_step_parameters() const {
  return generalized_adaptive_step_parameters_;
}

int FlightPlan::number_of_segments() const {
  return segments_.size();
}

DiscreteTrajectorySegmentIterator<Barycentric> FlightPlan::GetSegment(
    int const index) const {
  CHECK_LE(0, index);
  CHECK_LT(index, number_of_segments());
  return segments_[index];
}

DiscreteTrajectory<Barycentric> const& FlightPlan::GetAllSegments() const {
  return trajectory_;
}

DiscreteTrajectorySegmentIterator<Barycentric>
FlightPlan::GetSegmentAvoidingDeadlines(int index) {
  auto const status = RecomputeSegmentsAvoidingDeadlineIfNeeded();
  LOG_IF(INFO, !status.ok()) << "Unable to avoid deadline: " << status;
  return GetSegment(index);
}

DiscreteTrajectory<Barycentric> const&
FlightPlan::GetAllSegmentsAvoidingDeadlines() {
  auto const status = RecomputeSegmentsAvoidingDeadlineIfNeeded();
  LOG_IF(INFO, !status.ok()) << "Unable to avoid deadline: " << status;
  return GetAllSegments();
}

OrbitAnalyser::Analysis* FlightPlan::analysis(int coast_index) {
  if (coast_index > manœuvres_.size() - number_of_anomalous_manœuvres()) {
    // If the coast follows an anomalous manœuvre, no valid initial state was
    // available with which to request an analysis.
    return nullptr;
  }
  coast_analysers_[coast_index]->RefreshAnalysis();
  return coast_analysers_[coast_index]->analysis();
}

double FlightPlan::progress_of_analysis(int coast_index) const {
  if (coast_index > manœuvres_.size() - number_of_anomalous_manœuvres()) {
    return 0.0;
  }
  return coast_analysers_[coast_index]->progress_of_next_analysis();
}

void FlightPlan::EnableAnalysis(bool const enabled) {
  if (enabled != analysis_is_enabled_) {
    if (enabled) {
      // Request analysis of all non-anomalous coasts, and the first anomalous
      // segment if it is a coast.
      for (int index = 0;
           index < segments_.size() - std::max(0, anomalous_segments_ - 1);
           ++index) {
        if (index % 2 == 0) {
          auto const& coast = segments_[index];
          auto const& [first_time, first_degrees_of_freedom] = coast->front();
          coast_analysers_[index / 2]->RequestAnalysis(
              {.first_time = first_time,
               .first_degrees_of_freedom = first_degrees_of_freedom,
               .mission_duration = coast->back().time - first_time,
               .extended_mission_duration = desired_final_time_ - first_time});
        }
      }
    } else {
      // Immediately stop all the analyses.
      for (auto const& coast_analyser : coast_analysers_) {
        coast_analyser->Interrupt();
      }
    }
    analysis_is_enabled_ = enabled;
  }
}

void FlightPlan::WriteToMessage(
    not_null<serialization::FlightPlan*> const message) const {
  initial_mass_.WriteToMessage(message->mutable_initial_mass());
  initial_time_.WriteToMessage(message->mutable_initial_time());
  initial_degrees_of_freedom_.WriteToMessage(
      message->mutable_initial_degrees_of_freedom());
  desired_final_time_.WriteToMessage(message->mutable_desired_final_time());
  adaptive_step_parameters_.WriteToMessage(
      message->mutable_adaptive_step_parameters());
  generalized_adaptive_step_parameters_.WriteToMessage(
      message->mutable_generalized_adaptive_step_parameters());
  for (auto const& manœuvre : manœuvres_) {
    manœuvre.WriteToMessage(message->add_manoeuvre());
  }
}

std::unique_ptr<FlightPlan> FlightPlan::ReadFromMessage(
    serialization::FlightPlan const& message,
    not_null<Ephemeris<Barycentric>*> const ephemeris) {
  Instant initial_time = Instant::ReadFromMessage(message.initial_time());
  std::unique_ptr<DegreesOfFreedom<Barycentric>> initial_degrees_of_freedom;
  CHECK(message.has_adaptive_step_parameters());
  auto const adaptive_step_parameters =
      std::make_unique<Ephemeris<Barycentric>::AdaptiveStepParameters>(
          Ephemeris<Barycentric>::AdaptiveStepParameters::ReadFromMessage(
              message.adaptive_step_parameters()));

  bool const is_pre_erdős =
      !message.has_generalized_adaptive_step_parameters();
  LOG_IF(WARNING, is_pre_erdős) << "Reading pre-Erdős FlightPlan";
  std::unique_ptr<Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters>
      generalized_adaptive_step_parameters;
  if (is_pre_erdős) {
    generalized_adaptive_step_parameters = std::make_unique<
        Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters>(
        DefaultBurnParameters());
  } else {
    generalized_adaptive_step_parameters = std::make_unique<
        Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters>(
        Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters::
            ReadFromMessage(message.generalized_adaptive_step_parameters()));
  }

  initial_degrees_of_freedom =
      std::make_unique<DegreesOfFreedom<Barycentric>>(
          DegreesOfFreedom<Barycentric>::ReadFromMessage(
              message.initial_degrees_of_freedom()));

  auto flight_plan = std::make_unique<FlightPlan>(
      Mass::ReadFromMessage(message.initial_mass()),
      initial_time,
      *initial_degrees_of_freedom,
      Instant::ReadFromMessage(message.desired_final_time()),
      ephemeris,
      *adaptive_step_parameters,
      *generalized_adaptive_step_parameters);

  for (int i = 0; i < message.manoeuvre_size(); ++i) {
    auto const& manoeuvre = message.manoeuvre(i);
    flight_plan->manœuvres_.push_back(
        NavigationManœuvre::ReadFromMessage(manoeuvre, ephemeris));
    flight_plan->coast_analysers_.push_back(make_not_null_unique<OrbitAnalyser>(
        flight_plan->ephemeris_, DefaultHistoryParameters()));
  }
  // We need to forcefully prolong, otherwise we might exceed the ephemeris
  // step limit while recomputing the segments and make the flight plan
  // anomalous for no good reason.
  flight_plan->ephemeris_->Prolong(flight_plan->desired_final_time_)
      .IgnoreError();
  absl::Status const status = flight_plan->RecomputeAllSegments();
  LOG_IF(INFO, flight_plan->anomalous_segments_ > 0)
      << "Loading a flight plan with " << flight_plan->anomalous_segments_
      << " anomalous segments and status " << status << "\n"
      << message.DebugString();

  return flight_plan;
}

FlightPlan::FlightPlan()
    : initial_degrees_of_freedom_(Barycentric::origin, Barycentric::unmoving),
      ephemeris_(make_not_null<Ephemeris<Barycentric>*>()),
      adaptive_step_parameters_(
          EmbeddedExplicitRungeKuttaNyströmIntegrator<
              DormandالمكاوىPrince1986RKN434FM,
              Ephemeris<Barycentric>::NewtonianMotionEquation>(),
          /*max_steps=*/1,
          /*length_integration_tolerance=*/1 * Metre,
          /*speed_integration_tolerance=*/1 * Metre / Second),
      generalized_adaptive_step_parameters_(
          EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<
              Fine1987RKNG34,
              Ephemeris<Barycentric>::GeneralizedNewtonianMotionEquation>(),
          /*max_steps=*/1,
          /*length_integration_tolerance=*/1 * Metre,
          /*speed_integration_tolerance=*/1 * Metre / Second) {}

absl::Status FlightPlan::RecomputeAllSegments() {
  PopLastSegments(segments_.size() - 1);
  ResetLastSegment();
  return ComputeSegments(manœuvres_.begin(),
                         manœuvres_.end(),
                         max_ephemeris_steps_per_frame);
}

absl::Status FlightPlan::RecomputeSegmentsAvoidingDeadlineIfNeeded() {
  if (anomalous_segments_ == 0 ||
      !absl::IsDeadlineExceeded(anomalous_status_)) {
    return absl::OkStatus();
  }

  int const anomalous_manœuvres = anomalous_segments_ / 2;
  auto it = manœuvres_.end();
  for (int i = 0; i < anomalous_manœuvres; ++i) {
    PopLastSegments(2);
    --it;
  }
  ResetLastSegment();

  // Ask for 0 steps because is this called often on the UI thread.
  return ComputeSegments(it, manœuvres_.end(), /*max_ephemeris_steps*/0);
}

absl::Status FlightPlan::BurnSegment(
    NavigationManœuvre const& manœuvre,
    DiscreteTrajectorySegmentIterator<Barycentric> const segment,
    std::int64_t const max_ephemeris_steps) {
  Instant const final_time = manœuvre.final_time();
  if (manœuvre.initial_time() < final_time) {
    // Make sure that the ephemeris covers the entire segment, reanimating and
    // waiting if necessary.
    Instant const starting_time = segment->back().time;
    if (starting_time < ephemeris_->t_min()) {
      ephemeris_->AwaitReanimation(starting_time);
    }

    if (manœuvre.is_inertially_fixed()) {
      return ephemeris_->FlowWithAdaptiveStep(
                             &trajectory_,
                             manœuvre.InertialIntrinsicAcceleration(),
                             final_time,
                             adaptive_step_parameters_,
                             max_ephemeris_steps);
    } else {
      return ephemeris_->FlowWithAdaptiveStep(
                             &trajectory_,
                             manœuvre.FrenetIntrinsicAcceleration(),
                             final_time,
                             generalized_adaptive_step_parameters_,
                             max_ephemeris_steps);
    }
  } else {
    return absl::OkStatus();
  }
}

absl::Status FlightPlan::CoastSegment(
    Instant const& desired_final_time,
    DiscreteTrajectorySegmentIterator<Barycentric> const segment,
    std::int64_t const max_ephemeris_steps) {
  // Make sure that the ephemeris covers the entire segment, reanimating and
  // waiting if necessary.
  Instant const starting_time = segment->back().time;
  if (starting_time < ephemeris_->t_min()) {
    ephemeris_->AwaitReanimation(starting_time);
  }

  return ephemeris_->FlowWithAdaptiveStep(
                         &trajectory_,
                         Ephemeris<Barycentric>::NoIntrinsicAcceleration,
                         desired_final_time,
                         adaptive_step_parameters_,
                         max_ephemeris_steps);
}

absl::Status FlightPlan::ComputeSegments(
    std::vector<NavigationManœuvre>::iterator const begin,
    std::vector<NavigationManœuvre>::iterator const end,
    std::int64_t const max_ephemeris_steps) {
  CHECK(!segments_.empty());
  if (anomalous_segments_ == 0) {
    anomalous_status_ = absl::OkStatus();
  }
  absl::Status overall_status = anomalous_status_;
  for (auto it = begin; it != end; ++it) {
    auto& manœuvre = *it;
    auto& coast = segments_.back();
    manœuvre.clear_coasting_trajectory();

    if (anomalous_segments_ == 0) {
      absl::Status const status = CoastSegment(manœuvre.initial_time(),
                                               coast,
                                               max_ephemeris_steps);
      if (status.ok()) {
        manœuvre.set_coasting_trajectory(coast);
      } else {
        overall_status.Update(status);
        anomalous_segments_ = 1;
        anomalous_status_ = status;
      }
      if (analysis_is_enabled_) {
        auto const& [first_time, first_degrees_of_freedom] = coast->front();
        auto& analyser = coast_analysers_[it - manœuvres_.begin()];
        // Interrupt the analyses for the coasts that are being recomputed as we
        // won't use their result.
        analyser->Interrupt();
        analyser->RequestAnalysis(
            {.first_time = first_time,
             .first_degrees_of_freedom = first_degrees_of_freedom,
             .mission_duration = coast->back().time - first_time,
             .extended_mission_duration = desired_final_time_ - first_time});
      }
    }

    AddLastSegment();

    if (anomalous_segments_ == 0) {
      auto& burn = segments_.back();
      absl::Status const status = BurnSegment(manœuvre,
                                              burn,
                                              max_ephemeris_steps);
      if (!status.ok()) {
        overall_status.Update(status);
        anomalous_segments_ = 1;
        anomalous_status_ = status;
      }
    }

    AddLastSegment();
  }
  if (anomalous_segments_ == 0) {
    // If the desired end time is before the end of the last burn, move it to
    // that point.  Otherwise we might try to integrate towards the past in
    // CoastSegment.  Also, it's the user-friendly thing to do: no point in
    // having to extend the flight plan by hand.
    desired_final_time_ =
        std::max(desired_final_time_, segments_.back()->t_max());
    if (analysis_is_enabled_) {
      auto const& [first_time, first_degrees_of_freedom] =
          segments_.back()->front();
      auto& analyser = coast_analysers_.back();
      // Interrupt the analysis for the last coast as we won't use its result.
      analyser->Interrupt();
      analyser->RequestAnalysis(
          {.first_time = first_time,
           .first_degrees_of_freedom = first_degrees_of_freedom,
           .mission_duration = desired_final_time_ - first_time});
    }
    absl::Status const status = CoastSegment(desired_final_time_,
                                             segments_.back(),
                                             max_ephemeris_steps);
    if (!status.ok()) {
      overall_status.Update(status);
      anomalous_segments_ = 1;
      anomalous_status_ = status;
    }
  }
  return overall_status;
}

void FlightPlan::AddLastSegment() {
  segments_.emplace_back(trajectory_.NewSegment());
  if (anomalous_segments_ > 0) {
    ++anomalous_segments_;
  }
}

void FlightPlan::ResetLastSegment() {
  auto const& last_segment = segments_.back();
  trajectory_.ForgetAfter(std::next(last_segment->begin()));
  if (anomalous_segments_ == 1) {
    anomalous_segments_ = 0;
  }
}

void FlightPlan::PopLastSegments(std::int64_t const count) {
  CHECK_EQ(0, count % 2) << count;
  // It is important that the segments be destroyed in (reverse chronological)
  // order of the forks.
  for (std::int64_t i = 0; i < count; ++i) {
    auto& last_segment = segments_.back();
    trajectory_.DeleteSegments(last_segment);
    segments_.pop_back();
    if (anomalous_segments_ > 0) {
      --anomalous_segments_;
    }
  }
}

void FlightPlan::PopSegmentsAffectedByManœuvre(int const index) {
  // We will keep, for each manœuvre in [0, index[, its burn and the coast
  // preceding it, as well as the coast preceding manœuvre `index`.
  int const segments_kept = 2 * index + 1;
  PopLastSegments(number_of_segments() - segments_kept);
  ResetLastSegment();
}

void FlightPlan::UpdateInitialMassOfManœuvresAfter(int const index) {
  if (index >= manœuvres_.size()) {
    return;
  }
  Mass initial_mass = manœuvres_[index].final_mass();
  for (int i = index + 1; i < manœuvres_.size(); ++i) {
    manœuvres_[i] = NavigationManœuvre(initial_mass, manœuvres_[i].burn());
    initial_mass = manœuvres_[i].final_mass();
  }
}

void FlightPlan::MakeProlongator(Instant const& prolongation_time) {
  // A helper lambda to swallow the status of RETURN_IF_STOPPED.
  auto const prolong_with_status = [this](Instant const& prolongation_time) {
    // The loop makes sure that we give the main thread a chance to call
    // methods of the ephemeris.  The call to `Prolong` below is expected
    // to take about 40 ms.
    do {
      RETURN_IF_STOPPED;
      std::this_thread::sleep_for(20ms);
      RETURN_IF_STOPPED;
      // Do not return on an error here because it could be a deadline exceeded.
      ephemeris_->Prolong(prolongation_time, max_ephemeris_steps_per_frame)
          .IgnoreError();
    } while (ephemeris_->t_max() < prolongation_time);
    return absl::OkStatus();
  };

  if (prolongation_time < last_prolongation_time_) {
    // The desired prolongation became shorter, just kill the prolongator
    // thread.  We may recreate it below, but shorter.
    prolongator_ = jthread();
  }
  if (ephemeris_->t_max() < prolongation_time) {
    // The ephemeris is too short, start a thread to prolong it.  Note that we
    // must copy `prolong_with_status` since it's called on another thread.
    last_prolongation_time_ = prolongation_time;
    prolongator_ =
        MakeStoppableThread([prolong_with_status, prolongation_time]() {
          prolong_with_status(prolongation_time).IgnoreError();
        });
  }
}

Instant FlightPlan::start_of_last_coast() const {
  return manœuvres_.empty() ? initial_time_ : manœuvres_.back().final_time();
}

Instant FlightPlan::start_of_burn(int const index) const {
  return index == manœuvres_.size() ? desired_final_time_
                                    : manœuvres_[index].initial_time();
}

Instant FlightPlan::start_of_next_burn(int const index) const {
  return start_of_burn(index + 1);
}

Instant FlightPlan::start_of_previous_coast(int const index) const {
  return index == 0 ? initial_time_ : manœuvres_[index - 1].final_time();
}

}  // namespace internal
}  // namespace _flight_plan
}  // namespace ksp_plugin
}  // namespace principia
