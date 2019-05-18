
#include "ksp_plugin/flight_plan.hpp"

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
  RecomputeSegments(manœuvres_, segments_);
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

NavigationManœuvre const& FlightPlan::GetManœuvre(int const index) const {
  CHECK_LE(0, index);
  CHECK_LT(index, number_of_manœuvres());
  return manœuvres_[index];
}

bool FlightPlan::Append(NavigationManœuvre::Burn const& burn) {
  NavigationManœuvre const manœuvre(
      manœuvres_.empty() ? initial_mass_ : manœuvres_.back().final_mass(),
      burn);
  if (manœuvre.FitsBetween(start_of_last_coast(), desired_final_time_) &&
      !manœuvre.IsSingular()) {
    DiscreteTrajectory<Barycentric>* recomputed_last_coast =
        CoastIfReachesManœuvreInitialTime(last_coast(), manœuvre);
    if (recomputed_last_coast != nullptr) {
      ReplaceLastSegment(recomputed_last_coast);
      Append(manœuvre);
      return true;
    }
  }
  return false;
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

void FlightPlan::RemoveLast() {
  CHECK(!manœuvres_.empty());
  manœuvres_.pop_back();
  PopLastSegment();  // Last coast.
  PopLastSegment();  // Last burn.
  ResetLastSegment();
  RecomputeSegments(manœuvres_, segments_);
}

Status FlightPlan::Replace(NavigationManœuvre::Burn const& burn,
                           int const index) {
  CHECK_LE(0, index);
  CHECK_LT(index, number_of_manœuvres());
  NavigationManœuvre const manœuvre(manœuvres_[index].initial_mass(),
                                    burn);
  if (manœuvre.IsSingular()) {
    return Status(Error::OUT_OF_RANGE, "Singular");
  }
  if (!manœuvre.FitsBetween(start_of_previous_coast(index),
                            start_of_next_burn(index))) {
    return Status(Error::INVALID_ARGUMENT, "Doesn't fit");
  }

  DiscreteTrajectory<Barycentric>* replaced_coast = previous_coast(index);
  DiscreteTrajectory<Barycentric>& replaced_coast_parent =
      *replaced_coast->parent();
  // TODO(phl): Fork as late as possible.
  DiscreteTrajectory<Barycentric>* tentative_coast =
      replaced_coast_parent.NewForkWithoutCopy(replaced_coast->Fork().time());
  std::vector<not_null<DiscreteTrajectory<Barycentric>*>> tentative_segments = {
      tentative_coast};
  std::vector<NavigationManœuvre> tentative_manœuvres = {manœuvre};
  Mass initial_mass = manœuvre.final_mass();
  for (int i = index + 1; i < manœuvres_.size(); ++i) {
    NavigationManœuvre const tentative_manœuvre(initial_mass,
                                                manœuvres_[i].burn());
    tentative_manœuvres.push_back(tentative_manœuvre);
    initial_mass = tentative_manœuvre.final_mass();
  }
  Status const status = RecomputeSegments(tentative_manœuvres,
                                          tentative_segments);
  if (status.ok()) {
    // Keep the tentative manœuvres and segments: copy them into the vectors.
    manœuvres_.erase(manœuvres_.cbegin() + index,
                     manœuvres_.cend());
    std::copy(tentative_manœuvres.cbegin(),
              tentative_manœuvres.cend(),
              std::back_inserter(manœuvres_));
    segments_.erase(segments_.cbegin() + 2 * index,
                    segments_.cend());
    std::copy(tentative_segments.cbegin(),
              tentative_segments.cend(),
              std::back_inserter(segments_));
    replaced_coast_parent.DeleteFork(replaced_coast);
  } else {
    // Drop the tentative manœuvres and segments.
    replaced_coast_parent.DeleteFork(tentative_coast);
  }
  return status;
}

bool FlightPlan::ReplaceLast(NavigationManœuvre::Burn const& burn) {
#if 1
  return Replace(burn, manœuvres_.size() - 1).ok();
#else
  CHECK(!manœuvres_.empty());
  NavigationManœuvre const manœuvre(manœuvres_.back().initial_mass(),
                                    burn);
  if (manœuvre.FitsBetween(start_of_penultimate_coast(), desired_final_time_) &&
      !manœuvre.IsSingular()) {
    DiscreteTrajectory<Barycentric>* recomputed_penultimate_coast =
        CoastIfReachesManœuvreInitialTime(penultimate_coast(), manœuvre);
    if (recomputed_penultimate_coast != nullptr) {
      manœuvres_.pop_back();
      PopLastSegment();  // Last coast.
      PopLastSegment();  // Last burn.
      ReplaceLastSegment(recomputed_penultimate_coast);
      Append(manœuvre);
      return true;
    }
  }
  return false;
#endif
}

bool FlightPlan::SetDesiredFinalTime(Instant const& desired_final_time) {
  if (start_of_last_coast() > desired_final_time) {
    return false;
  } else {
    desired_final_time_ = desired_final_time;
    ResetLastSegment();
    RecomputeSegments(manœuvres_, segments_);
    return true;//TODO(phl)OK?
  }
}

Ephemeris<Barycentric>::AdaptiveStepParameters const&
FlightPlan::adaptive_step_parameters() const {
  return adaptive_step_parameters_;
}

Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters const&
FlightPlan::generalized_adaptive_step_parameters() const {
  return generalized_adaptive_step_parameters_;
}

bool FlightPlan::SetAdaptiveStepParameters(
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        adaptive_step_parameters,
    Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters const&
        generalized_adaptive_step_parameters) {
  auto const original_adaptive_step_parameters = adaptive_step_parameters_;
  auto const original_generalized_adaptive_step_parameters =
      generalized_adaptive_step_parameters_;
  adaptive_step_parameters_ = adaptive_step_parameters;
  generalized_adaptive_step_parameters_ = generalized_adaptive_step_parameters;
  if (RecomputeSegments()) {
    return true;
  } else {
    // If the recomputation fails, leave this place as clean as we found it.
    adaptive_step_parameters_ = original_adaptive_step_parameters;
    generalized_adaptive_step_parameters_ =
        original_generalized_adaptive_step_parameters;
    CHECK(RecomputeSegments());
    return false;
  }
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
  // step limit while recomputing the segments and fail the check.
  flight_plan->ephemeris_->Prolong(flight_plan->desired_final_time_);
  CHECK(flight_plan->RecomputeSegments()) << message.DebugString();

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

void FlightPlan::Append(NavigationManœuvre const& manœuvre) {
  manœuvres_.push_back(manœuvre);
  CHECK_EQ(manœuvre.initial_time(), segments_.back()->last().time());
  RecomputeSegments(manœuvres_, segments_);
}

bool FlightPlan::RecomputeSegments() {
  // It is important that the segments be destroyed in (reverse chronological)
  // order of the forks.
  while (segments_.size() > 1) {
    PopLastSegment();
  }
  Status const status = RecomputeSegments(manœuvres_, segments_);
  return status.ok() && anomalous_segments_ <= 2;
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
                             max_ephemeris_steps_per_frame,
                             /*last_point_only=*/false);
    } else {
      return ephemeris_->FlowWithAdaptiveStep(
                             segment,
                             manœuvre.FrenetIntrinsicAcceleration(),
                             final_time,
                             generalized_adaptive_step_parameters_,
                             max_ephemeris_steps_per_frame,
                             /*last_point_only=*/false);
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
                         max_ephemeris_steps_per_frame,
                         /*last_point_only=*/false);
}

Status FlightPlan::RecomputeSegments(
    std::vector<NavigationManœuvre>& manœuvres,
    std::vector<not_null<DiscreteTrajectory<Barycentric>*>>& segments) {
  CHECK(!segments.empty());
  Status overall_status;
  int anomalous_segments = 0;
  for (auto& manœuvre : manœuvres) {
    auto& coast = segments.back();
    manœuvre.set_coasting_trajectory(coast);

    if (anomalous_segments == 0) {
      Status const status = CoastSegment(manœuvre.initial_time(), coast);
      if (!status.ok()) {
        overall_status.Update(status);
        anomalous_segments = 1;
      }
    }

    AddSegment(segments);
    if (anomalous_segments > 0) {
      ++anomalous_segments;
    }

    if (anomalous_segments == 0) {
      auto& burn = segments.back();
      Status const status = BurnSegment(manœuvre, burn);
      if (!status.ok()) {
        overall_status.Update(status);
        anomalous_segments = 1;
      }
    }

    AddSegment(segments);
    if (anomalous_segments > 0) {
      ++anomalous_segments;
    }
  }
  if (anomalous_segments == 0) {
    Status const status = CoastSegment(desired_final_time_, segments.back());
    if (!status.ok()) {
      overall_status.Update(status);
      anomalous_segments = 1;
    }
  }
  return overall_status;
}

void FlightPlan::ReplaceLastSegment(
    not_null<DiscreteTrajectory<Barycentric>*> const segment) {
  CHECK_EQ(segment->parent(), segments_.back()->parent());
  CHECK_EQ(segment->Fork().time(), segments_.back()->Fork().time());
  PopLastSegment();
  // |segment| must not be anomalous, so it cannot not follow an anomalous
  // segment.
  CHECK_EQ(0, anomalous_segments_);
  segments_.emplace_back(segment);
}

void FlightPlan::AddSegment(
    std::vector<not_null<DiscreteTrajectory<Barycentric>*>>& segments) {
  segments.emplace_back(segments.back()->NewForkAtLast());
}

void FlightPlan::ResetLastSegment() {
  segments_.back()->ForgetAfter(segments_.back()->Fork().time());
  }
}

void FlightPlan::PopLastSegment() {
  DiscreteTrajectory<Barycentric>* trajectory = segments_.back();
  CHECK(!trajectory->is_root());
  trajectory->parent()->DeleteFork(trajectory);
  segments_.pop_back();
}

DiscreteTrajectory<Barycentric>* FlightPlan::CoastIfReachesManœuvreInitialTime(
    DiscreteTrajectory<Barycentric>& coast,
    NavigationManœuvre const& manœuvre) {
  DiscreteTrajectory<Barycentric>* recomputed_coast =
      coast.parent()->NewForkWithoutCopy(coast.Fork().time());
  bool const reached_manœuvre_initial_time =
      ephemeris_->FlowWithAdaptiveStep(
          recomputed_coast,
          Ephemeris<Barycentric>::NoIntrinsicAcceleration,
          manœuvre.initial_time(),
          adaptive_step_parameters_,
          max_ephemeris_steps_per_frame,
          /*last_point_only=*/false).ok();
  if (!reached_manœuvre_initial_time) {
    recomputed_coast->parent()->DeleteFork(recomputed_coast);
  }
  return recomputed_coast;
}

Instant FlightPlan::start_of_last_coast() const {
  return manœuvres_.empty() ? initial_time_ : manœuvres_.back().final_time();
}

Instant FlightPlan::start_of_penultimate_coast() const {
  CHECK(!manœuvres_.empty());
  return manœuvres_.size() == 1
             ? initial_time_
             : manœuvres_[manœuvres_.size() - 2].final_time();
}

Instant FlightPlan::start_of_next_burn(int const index) const {
  return index == manœuvres_.size() - 1
             ? desired_final_time_
             : manœuvres_[index + 1].initial_time();
}

Instant FlightPlan::start_of_previous_coast(int const index) const {
  return index == 0 ? initial_time_ : manœuvres_[index - 1].final_time();
}

DiscreteTrajectory<Barycentric>& FlightPlan::last_coast() {
  return *segments_.back();
}

DiscreteTrajectory<Barycentric>& FlightPlan::penultimate_coast() {
  // The penultimate coast is the antepenultimate segment.
  return *segments_[segments_.size() - 3];
}

not_null<DiscreteTrajectory<Barycentric>*> FlightPlan::previous_coast(
    int const index) {
  return segments_[2 * index];
}

}  // namespace internal_flight_plan
}  // namespace ksp_plugin
}  // namespace principia
