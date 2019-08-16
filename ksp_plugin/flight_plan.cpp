
#include "ksp_plugin/flight_plan.hpp"

#include <algorithm>
#include <optional>
#include <vector>

#include "integrators/embedded_explicit_generalized_runge_kutta_nyström_integrator.hpp"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/methods.hpp"
#include "ksp_plugin/integrators.hpp"
#include "testing_utilities/make_not_null.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_flight_plan {

using base::Error;
using base::make_not_null_unique;
using base::Status;
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

inline Status const BadDesiredFinalTime() {
  return Status(FlightPlan::bad_desired_final_time, "Bad desired final time");
}

inline Status const DoesNotFit() {
  return Status(FlightPlan::does_not_fit, "Does not fit");
}

inline Status const Singular() {
  return Status(FlightPlan::singular, "Singular");
}

FlightPlan::FlightPlan(
    Mass const& initial_mass,
    Instant const& initial_time,
    DegreesOfFreedom<Barycentric> const& initial_degrees_of_freedom,
    Instant const& desired_final_time,
    not_null<Ephemeris<Barycentric>*> const ephemeris,
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        adaptive_step_parameters,
    Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters const&
        generalized_adaptive_step_parameters)
    : initial_mass_(initial_mass),
      initial_time_(initial_time),
      initial_degrees_of_freedom_(initial_degrees_of_freedom),
      desired_final_time_(desired_final_time),
      root_(make_not_null_unique<DiscreteTrajectory<Barycentric>>()),
      ephemeris_(ephemeris),
      adaptive_step_parameters_(adaptive_step_parameters),
      generalized_adaptive_step_parameters_(
          generalized_adaptive_step_parameters) {
  CHECK(desired_final_time_ >= initial_time_);

  // Set the (single) point of the root.
  root_->Append(initial_time_, initial_degrees_of_freedom_);

  // Create a fork for the first coasting trajectory.
  segments_.emplace_back(root_->NewForkWithoutCopy(initial_time_));
  CHECK(manœuvres_.empty());
  ComputeSegments(manœuvres_.begin(), manœuvres_.end());
}

Instant FlightPlan::initial_time() const {
  return initial_time_;
}

Instant FlightPlan::actual_final_time() const {
  return segments_.back()->last().time();
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

Status FlightPlan::Append(NavigationManœuvre::Burn const& burn) {
  NavigationManœuvre const manœuvre(
      manœuvres_.empty() ? initial_mass_ : manœuvres_.back().final_mass(),
      burn);
  if (manœuvre.IsSingular()) {
    return Singular();
  }
  if (!manœuvre.FitsBetween(start_of_last_coast(), desired_final_time_)) {
    return DoesNotFit();
  }
  // Reset the last coast and integrate it followed by a burn followed by a new
  // last coast;
  manœuvres_.push_back(manœuvre);
  ResetLastSegment();
  return ComputeSegments(manœuvres_.begin() + manœuvres_.size() - 1,
                         manœuvres_.end());
}

void FlightPlan::ForgetBefore(Instant const& time,
                              std::function<void()> const& on_empty) {
  // Find the first segment to keep.  Note that incrementing by 2 ensures that
  // we only look at coasts.
  std::optional<int> first_to_keep;
  for (int i = 0; i < segments_.size(); i += 2) {
    if (time <= segments_[i]->last().time()) {
      first_to_keep = i;
      break;
    }
  }
  if (!first_to_keep) {
    // The entire flight plan needs to go away.
    on_empty();
    return;
  }

  // Detach the first coast to keep, truncate its beginning, and reattach it
  // to a new root.
  std::unique_ptr<DiscreteTrajectory<Barycentric>> new_first_coast =
      segments_[*first_to_keep]->DetachFork();
  new_first_coast->ForgetBefore(time);
  root_ = make_not_null_unique<DiscreteTrajectory<Barycentric>>();
  root_->Append(new_first_coast->Begin().time(),
                new_first_coast->Begin().degrees_of_freedom());
  root_->AttachFork(std::move(new_first_coast));

  // Remove from the vectors the trajectories and manœuvres that we don't want
  // to keep.
  segments_.erase(segments_.cbegin(),
                  segments_.cbegin() + *first_to_keep);
  manœuvres_.erase(manœuvres_.cbegin(),
                   manœuvres_.cbegin() + *first_to_keep / 2);

  auto const root_begin = root_->Begin();
  initial_time_ = root_begin.time();
  initial_degrees_of_freedom_ = root_begin.degrees_of_freedom();
}

Status FlightPlan::RemoveLast() {
  CHECK(!manœuvres_.empty());
  manœuvres_.pop_back();
  PopLastSegment();  // Last coast.
  PopLastSegment();  // Last burn.
  // Clear and recompute the last coast.
  ResetLastSegment();
  return ComputeSegments(manœuvres_.end(), manœuvres_.end());
}

Status FlightPlan::Replace(NavigationManœuvre::Burn const& burn,
                           int const index) {
  CHECK_LE(0, index);
  CHECK_LT(index, number_of_manœuvres());
  NavigationManœuvre const manœuvre(manœuvres_[index].initial_mass(),
                                    burn);
  if (manœuvre.IsSingular()) {
    return Singular();
  }
  if (!manœuvre.FitsBetween(start_of_previous_coast(index),
                            start_of_next_burn(index))) {
    return DoesNotFit();
  }

  // Replace the manœuvre at position |index| and rebuild all the ones that
  // follow as they may have a different initial mass.  Also pop the segments
  // that we'll recompute.
  manœuvres_[index] = manœuvre;
  PopLastSegment();  // Last coast.
  PopLastSegment();  // Last burn.

  Mass initial_mass = manœuvre.final_mass();
  for (int i = index + 1; i < manœuvres_.size(); ++i) {
    manœuvres_[i] = NavigationManœuvre(initial_mass, manœuvres_[i].burn());
    initial_mass = manœuvres_[i].final_mass();
    PopLastSegment();  // Last coast.
    PopLastSegment();  // Last burn.
  }

  // TODO(phl): Recompute as late as possible.
  // At this point the last coast is the one to which |manœuvre| gets attached.
  // Clear it and recompute everything that follows.
  ResetLastSegment();
  return ComputeSegments(manœuvres_.begin() + index, manœuvres_.end());
}

Status FlightPlan::SetDesiredFinalTime(Instant const& desired_final_time) {
  if (desired_final_time < start_of_last_coast()) {
    return BadDesiredFinalTime();
  }
  desired_final_time_ = desired_final_time;
  // Reset the last coast and recompute it.
  ResetLastSegment();
  return ComputeSegments(manœuvres_.end(), manœuvres_.end());
}

Status FlightPlan::SetAdaptiveStepParameters(
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
  end = segments_[index]->End();
}

void FlightPlan::GetAllSegments(
    DiscreteTrajectory<Barycentric>::Iterator& begin,
    DiscreteTrajectory<Barycentric>::Iterator& end) const {
  begin = segments_.back()->Find(segments_.front()->Fork().time());
  end = segments_.back()->End();
  CHECK(begin != end);
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
  }
  // We need to forcefully prolong, otherwise we might exceed the ephemeris
  // step limit while recomputing the segments and make the flight plan
  // anomalous for no good reason.
  flight_plan->ephemeris_->Prolong(flight_plan->desired_final_time_);
  Status const status = flight_plan->RecomputeAllSegments();
  LOG_IF(INFO, flight_plan->anomalous_segments_ > 0)
      << "Loading a flight plan with " << flight_plan->anomalous_segments_
      << " anomalous segments and status " << status << "\n"
      << message.DebugString();

  return flight_plan;
}

FlightPlan::FlightPlan()
    : initial_degrees_of_freedom_(Barycentric::origin, Velocity<Barycentric>()),
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

Status FlightPlan::RecomputeAllSegments() {
  // It is important that the segments be destroyed in (reverse chronological)
  // order of the forks.
  while (segments_.size() > 1) {
    PopLastSegment();
  }
  ResetLastSegment();
  return ComputeSegments(manœuvres_.begin(), manœuvres_.end());
}

Status FlightPlan::BurnSegment(
    NavigationManœuvre const& manœuvre,
    not_null<DiscreteTrajectory<Barycentric>*> const segment) {
  Instant const final_time = manœuvre.final_time();
  if (manœuvre.initial_time() < final_time) {
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
    return Status::OK;
  }
}

Status FlightPlan::CoastSegment(
    Instant const& desired_final_time,
    not_null<DiscreteTrajectory<Barycentric>*> const segment) {
  return ephemeris_->FlowWithAdaptiveStep(
                         segment,
                         Ephemeris<Barycentric>::NoIntrinsicAcceleration,
                         desired_final_time,
                         adaptive_step_parameters_,
                         max_ephemeris_steps_per_frame);
}

Status FlightPlan::ComputeSegments(
    std::vector<NavigationManœuvre>::iterator const begin,
    std::vector<NavigationManœuvre>::iterator const end) {
  CHECK(!segments_.empty());
  if (anomalous_segments_ == 0) {
    anomalous_status_ = Status::OK;
  }
  Status overall_status = anomalous_status_;
  for (auto it = begin; it != end; ++it) {
    auto& manœuvre = *it;
    auto& coast = segments_.back();
    manœuvre.set_coasting_trajectory(coast);

    if (anomalous_segments_ == 0) {
      Status const status = CoastSegment(manœuvre.initial_time(), coast);
      if (!status.ok()) {
        overall_status.Update(status);
        anomalous_segments_ = 1;
        anomalous_status_ = status;
      }
    }

    AddLastSegment();

    if (anomalous_segments_ == 0) {
      auto& burn = segments_.back();
      Status const status = BurnSegment(manœuvre, burn);
      if (!status.ok()) {
        overall_status.Update(status);
        anomalous_segments_ = 1;
        anomalous_status_ = status;
      }
    }

    AddLastSegment();
  }
  if (anomalous_segments_ == 0) {
    Status const status = CoastSegment(desired_final_time_, segments_.back());
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
  segments_.back()->ForgetAfter(segments_.back()->Fork().time());
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

Instant FlightPlan::start_of_last_coast() const {
  return manœuvres_.empty() ? initial_time_ : manœuvres_.back().final_time();
}

Instant FlightPlan::start_of_next_burn(int const index) const {
  return index == manœuvres_.size() - 1
             ? desired_final_time_
             : manœuvres_[index + 1].initial_time();
}

Instant FlightPlan::start_of_previous_coast(int const index) const {
  return index == 0 ? initial_time_ : manœuvres_[index - 1].final_time();
}

}  // namespace internal_flight_plan
}  // namespace ksp_plugin
}  // namespace principia
