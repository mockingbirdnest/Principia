
#include "ksp_plugin/flight_plan.hpp"

#include <algorithm>
#include <optional>
#include <utility>
#include <vector>

#include "integrators/embedded_explicit_generalized_runge_kutta_nyström_integrator.hpp"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/methods.hpp"
#include "ksp_plugin/integrators.hpp"
#include "testing_utilities/make_not_null.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_flight_plan {

using base::make_not_null_unique;
using geometry::Position;
using geometry::Vector;
using geometry::Velocity;
using integrators::EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator;
using integrators::EmbeddedExplicitRungeKuttaNyströmIntegrator;
using integrators::methods::DormandالمكاوىPrince1986RKN434FM;
using integrators::methods::Fine1987RKNG34;
using quantities::Acceleration;
using quantities::si::Metre;
using quantities::si::Second;

inline absl::Status BadDesiredFinalTime() {
  return absl::Status(FlightPlan::bad_desired_final_time,
                      "Bad desired final time");
}

inline absl::Status DoesNotFit() {
  return absl::Status(FlightPlan::does_not_fit, "Does not fit");
}

inline absl::Status Singular() {
  return absl::Status(FlightPlan::singular, "Singular");
}

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
      desired_final_time_(desired_final_time),
      root_(make_not_null_unique<DiscreteTrajectory<Barycentric>>()),
      ephemeris_(ephemeris),
      adaptive_step_parameters_(std::move(adaptive_step_parameters)),
      generalized_adaptive_step_parameters_(
          std::move(generalized_adaptive_step_parameters)) {
  CHECK(desired_final_time_ >= initial_time_);

  // Set the (single) point of the root.
  root_->Append(initial_time_, initial_degrees_of_freedom_);

  // Create a fork for the first coasting trajectory.
  segments_.emplace_back(root_->NewForkWithoutCopy(initial_time_));
  coast_analysers_.push_back(make_not_null_unique<OrbitAnalyser>(
      ephemeris_, DefaultHistoryParameters()));
  CHECK(manœuvres_.empty());
  ComputeSegments(manœuvres_.begin(), manœuvres_.end());
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
    return Singular();
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
  return ComputeSegments(manœuvres_.begin() + index, manœuvres_.end());
}

absl::Status FlightPlan::Remove(int index) {
  CHECK_GE(index, 0);
  CHECK_LT(index, number_of_manœuvres());
  manœuvres_.erase(manœuvres_.begin() + index);
  coast_analysers_.erase(coast_analysers_.begin() + index + 1);
  UpdateInitialMassOfManœuvresAfter(index);
  PopSegmentsAffectedByManœuvre(index);
  return ComputeSegments(manœuvres_.begin() + index, manœuvres_.end());
}

absl::Status FlightPlan::Replace(NavigationManœuvre::Burn const& burn,
                                 int const index) {
  CHECK_LE(0, index);
  CHECK_LT(index, number_of_manœuvres());
  NavigationManœuvre const manœuvre(manœuvres_[index].initial_mass(),
                                    burn);
  if (manœuvre.IsSingular()) {
    return Singular();
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

  // Replace the manœuvre at position |index| and rebuild all the ones that
  // follow as they may have a different initial mass.  Also pop the segments
  // that we'll recompute.
  manœuvres_[index] = manœuvre;
  UpdateInitialMassOfManœuvresAfter(index);

  // TODO(phl): Recompute as late as possible.
  PopSegmentsAffectedByManœuvre(index);
  return ComputeSegments(manœuvres_.begin() + index, manœuvres_.end());
}

absl::Status FlightPlan::SetDesiredFinalTime(
    Instant const& desired_final_time) {
  if (desired_final_time < start_of_last_coast()) {
    return BadDesiredFinalTime();
  }
  desired_final_time_ = desired_final_time;
  // Reset the last coast and recompute it.
  ResetLastSegment();
  return ComputeSegments(manœuvres_.end(), manœuvres_.end());
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

void FlightPlan::GetSegment(
    int const index,
    DiscreteTrajectory<Barycentric>::Iterator& begin,
    DiscreteTrajectory<Barycentric>::Iterator& end) const {
  CHECK_LE(0, index);
  CHECK_LT(index, number_of_segments());
  begin = segments_[index]->Fork();
  end = segments_[index]->end();
}

void FlightPlan::GetAllSegments(
    DiscreteTrajectory<Barycentric>::Iterator& begin,
    DiscreteTrajectory<Barycentric>::Iterator& end) const {
  begin = segments_.back()->Find(segments_.front()->Fork()->time);
  end = segments_.back()->end();
  CHECK(begin != end);
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
  flight_plan->ephemeris_->Prolong(flight_plan->desired_final_time_);
  absl::Status const status = flight_plan->RecomputeAllSegments();
  LOG_IF(INFO, flight_plan->anomalous_segments_ > 0)
      << "Loading a flight plan with " << flight_plan->anomalous_segments_
      << " anomalous segments and status " << status << "\n"
      << message.DebugString();

  return flight_plan;
}

FlightPlan::FlightPlan()
    : initial_degrees_of_freedom_(Barycentric::origin, Barycentric::unmoving),
      root_(make_not_null_unique<DiscreteTrajectory<Barycentric>>()),
      ephemeris_(testing_utilities::make_not_null<Ephemeris<Barycentric>*>()),
      adaptive_step_parameters_(
          EmbeddedExplicitRungeKuttaNyströmIntegrator<
              DormandالمكاوىPrince1986RKN434FM,
              Position<Barycentric>>(),
          /*max_steps=*/1,
          /*length_integration_tolerance=*/1 * Metre,
          /*speed_integration_tolerance=*/1 * Metre / Second),
      generalized_adaptive_step_parameters_(
          EmbeddedExplicitGeneralizedRungeKuttaNyströmIntegrator<
              Fine1987RKNG34,
              Position<Barycentric>>(),
          /*max_steps=*/1,
          /*length_integration_tolerance=*/1 * Metre,
          /*speed_integration_tolerance=*/1 * Metre / Second) {}

absl::Status FlightPlan::RecomputeAllSegments() {
  // It is important that the segments be destroyed in (reverse chronological)
  // order of the forks.
  while (segments_.size() > 1) {
    PopLastSegment();
  }
  ResetLastSegment();
  return ComputeSegments(manœuvres_.begin(), manœuvres_.end());
}

absl::Status FlightPlan::BurnSegment(
    NavigationManœuvre const& manœuvre,
    not_null<DiscreteTrajectory<Barycentric>*> const segment) {
  Instant const final_time = manœuvre.final_time();
  if (manœuvre.initial_time() < final_time) {
    // Make sure that the ephemeris covers the entire segment, reanimating and
    // waiting if necessary.
    Instant const starting_time = segment->back().time;
    if (starting_time < ephemeris_->t_min()) {
      ephemeris_->RequestReanimation(starting_time);
      ephemeris_->WaitForReanimation(starting_time);
    }

    if (manœuvre.is_inertially_fixed()) {
      return ephemeris_->FlowWithAdaptiveStep(
                             segment,
                             manœuvre.InertialIntrinsicAcceleration(),
                             final_time,
                             adaptive_step_parameters_,
                             max_ephemeris_steps_per_frame);
    } else {
      return ephemeris_->FlowWithAdaptiveStep(
                             segment,
                             manœuvre.FrenetIntrinsicAcceleration(),
                             final_time,
                             generalized_adaptive_step_parameters_,
                             max_ephemeris_steps_per_frame);
    }
  } else {
    return absl::OkStatus();
  }
}

absl::Status FlightPlan::CoastSegment(
    Instant const& desired_final_time,
    not_null<DiscreteTrajectory<Barycentric>*> const segment) {
  // Make sure that the ephemeris covers the entire segment, reanimating and
  // waiting if necessary.
  Instant const starting_time = segment->back().time;
  if (starting_time < ephemeris_->t_min()) {
    ephemeris_->RequestReanimation(starting_time);
    ephemeris_->WaitForReanimation(starting_time);
  }

  return ephemeris_->FlowWithAdaptiveStep(
                         segment,
                         Ephemeris<Barycentric>::NoIntrinsicAcceleration,
                         desired_final_time,
                         adaptive_step_parameters_,
                         max_ephemeris_steps_per_frame);
}

absl::Status FlightPlan::ComputeSegments(
    std::vector<NavigationManœuvre>::iterator const begin,
    std::vector<NavigationManœuvre>::iterator const end) {
  CHECK(!segments_.empty());
  if (anomalous_segments_ == 0) {
    anomalous_status_ = absl::OkStatus();
  }
  absl::Status overall_status = anomalous_status_;
  for (auto it = begin; it != end; ++it) {
    auto& manœuvre = *it;
    auto& coast = segments_.back();
    manœuvre.set_coasting_trajectory(coast);

    if (anomalous_segments_ == 0) {
      absl::Status const status = CoastSegment(manœuvre.initial_time(), coast);
      if (!status.ok()) {
        overall_status.Update(status);
        anomalous_segments_ = 1;
        anomalous_status_ = status;
      }
      coast_analysers_[it - manœuvres_.begin()]->RequestAnalysis(
          {.first_time = coast->Fork()->time,
           .first_degrees_of_freedom = coast->Fork()->degrees_of_freedom,
           .mission_duration = coast->back().time - coast->Fork()->time,
           .extended_mission_duration =
               desired_final_time_ - coast->Fork()->time});
    }

    AddLastSegment();

    if (anomalous_segments_ == 0) {
      auto& burn = segments_.back();
      absl::Status const status = BurnSegment(manœuvre, burn);
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
    coast_analysers_.back()->RequestAnalysis(
        {.first_time = segments_.back()->Fork()->time,
         .first_degrees_of_freedom =
             segments_.back()->Fork()->degrees_of_freedom,
         .mission_duration =
             desired_final_time_ - segments_.back()->Fork()->time});
    absl::Status const status =
        CoastSegment(desired_final_time_, segments_.back());
    if (!status.ok()) {
      overall_status.Update(status);
      anomalous_segments_ = 1;
      anomalous_status_ = status;
    }
  }
  return overall_status;
}

void FlightPlan::AddLastSegment() {
  segments_.emplace_back(segments_.back()->NewForkAtLast());
  if (anomalous_segments_ > 0) {
    ++anomalous_segments_;
  }
}

void FlightPlan::ResetLastSegment() {
  segments_.back()->ForgetAfter(segments_.back()->Fork()->time);
  if (anomalous_segments_ == 1) {
    anomalous_segments_ = 0;
  }
}

void FlightPlan::PopLastSegment() {
  DiscreteTrajectory<Barycentric>* trajectory = segments_.back();
  CHECK(!trajectory->is_root());
  trajectory->parent()->DeleteFork(trajectory);
  segments_.pop_back();
  if (anomalous_segments_ > 0) {
    --anomalous_segments_;
  }
}

void FlightPlan::PopSegmentsAffectedByManœuvre(int const index) {
  // We will keep, for each manœuvre in [0, index[, its burn and the coast
  // preceding it, as well as the coast preceding manœuvre |index|.
  int const segments_kept = 2 * index + 1;
  while (number_of_segments() > segments_kept) {
    PopLastSegment();
  }
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

}  // namespace internal_flight_plan
}  // namespace ksp_plugin
}  // namespace principia
