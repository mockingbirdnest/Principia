
#include "ksp_plugin/flight_plan.hpp"

#include <optional>
#include <vector>

#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/methods.hpp"
#include "testing_utilities/make_not_null.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_flight_plan {

using base::make_not_null_unique;
using geometry::Position;
using geometry::Velocity;
using integrators::EmbeddedExplicitRungeKuttaNyströmIntegrator;
using integrators::methods::DormandالمكاوىPrince1986RKN434FM;
using quantities::si::Metre;
using quantities::si::Second;

FlightPlan::FlightPlan(
    Mass const& initial_mass,
    Instant const& initial_time,
    DegreesOfFreedom<Barycentric> const& initial_degrees_of_freedom,
    Instant const& desired_final_time,
    not_null<Ephemeris<Barycentric>*> const ephemeris,
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        adaptive_step_parameters)
    : initial_mass_(initial_mass),
      initial_time_(initial_time),
      initial_degrees_of_freedom_(initial_degrees_of_freedom),
      desired_final_time_(desired_final_time),
      root_(make_not_null_unique<DiscreteTrajectory<Barycentric>>()),
      ephemeris_(ephemeris),
      adaptive_step_parameters_(adaptive_step_parameters) {
  CHECK(desired_final_time_ >= initial_time_);

  // Set the (single) point of the root.
  root_->Append(initial_time_, initial_degrees_of_freedom_);

  // Create a fork for the first coasting trajectory.
  segments_.emplace_back(root_->NewForkWithoutCopy(initial_time_));
  CoastLastSegment(desired_final_time_);
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

bool FlightPlan::Append(Burn burn) {
  auto manœuvre =
      MakeNavigationManœuvre(
          std::move(burn),
          manœuvres_.empty() ? initial_mass_ : manœuvres_.back().final_mass());
  if (manœuvre.FitsBetween(start_of_last_coast(), desired_final_time_) &&
      !manœuvre.IsSingular()) {
    DiscreteTrajectory<Barycentric>* recomputed_last_coast =
        CoastIfReachesManœuvreInitialTime(last_coast(), manœuvre);
    if (recomputed_last_coast != nullptr) {
      ReplaceLastSegment(recomputed_last_coast);
      Append(std::move(manœuvre));
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
  root_->AttachFork(std::move(new_first_coast));

  // Remove from the vectors the trajectories and manœuvres that we don't want
  // to keep.
  segments_.erase(segments_.cbegin(),
                  segments_.cbegin() + *first_to_keep);
  // For some reason manœuvres_.erase() doesn't work because it wants to
  // copy, hence this dance.
  std::vector<NavigationManœuvre> m;
  std::move(manœuvres_.begin() + *first_to_keep / 2,
            manœuvres_.end(),
            std::back_inserter(m));
  manœuvres_.swap(m);

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
  CoastLastSegment(desired_final_time_);
}

bool FlightPlan::ReplaceLast(Burn burn) {
  CHECK(!manœuvres_.empty());
  auto manœuvre = MakeNavigationManœuvre(std::move(burn),
                                         manœuvres_.back().initial_mass());
  if (manœuvre.FitsBetween(start_of_penultimate_coast(), desired_final_time_) &&
      !manœuvre.IsSingular()) {
    DiscreteTrajectory<Barycentric>* recomputed_penultimate_coast =
        CoastIfReachesManœuvreInitialTime(penultimate_coast(), manœuvre);
    if (recomputed_penultimate_coast != nullptr) {
      manœuvres_.pop_back();
      PopLastSegment();  // Last coast.
      PopLastSegment();  // Last burn.
      ReplaceLastSegment(recomputed_penultimate_coast);
      Append(std::move(manœuvre));
      return true;
    }
  }
  return false;
}

bool FlightPlan::SetDesiredFinalTime(Instant const& desired_final_time) {
  if (start_of_last_coast() > desired_final_time) {
    return false;
  } else {
    desired_final_time_ = desired_final_time;
    ResetLastSegment();
    CoastLastSegment(desired_final_time_);
    return true;
  }
}

Ephemeris<Barycentric>::AdaptiveStepParameters const&
FlightPlan::adaptive_step_parameters() const {
  return adaptive_step_parameters_;
}

bool FlightPlan::SetAdaptiveStepParameters(
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        adaptive_step_parameters) {
  auto const original_adaptive_step_parameters = adaptive_step_parameters_;
  adaptive_step_parameters_ = adaptive_step_parameters;
  if (RecomputeSegments()) {
    return true;
  } else {
    // If the recomputation fails, leave this place as clean as we found it.
    adaptive_step_parameters_ = original_adaptive_step_parameters;
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
  for (auto const& manœuvre : manœuvres_) {
    manœuvre.WriteToMessage(message->add_manoeuvre());
  }
}

std::unique_ptr<FlightPlan> FlightPlan::ReadFromMessage(
    serialization::FlightPlan const& message,
    not_null<Ephemeris<Barycentric>*> const ephemeris) {
  std::unique_ptr<Ephemeris<Barycentric>::AdaptiveStepParameters>
      adaptive_step_parameters;
  Instant initial_time = Instant::ReadFromMessage(message.initial_time());
  std::unique_ptr<DegreesOfFreedom<Barycentric>> initial_degrees_of_freedom;
  CHECK(message.has_adaptive_step_parameters());
  adaptive_step_parameters =
      std::make_unique<Ephemeris<Barycentric>::AdaptiveStepParameters>(
          Ephemeris<Barycentric>::AdaptiveStepParameters::ReadFromMessage(
              message.adaptive_step_parameters()));
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
      *adaptive_step_parameters);

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
          /*speed_integration_tolerance=*/1 * Metre / Second) {}

void FlightPlan::Append(NavigationManœuvre manœuvre) {
  manœuvres_.emplace_back(std::move(manœuvre));
  {
    // Hide the moved-from |manœuvre|.
    NavigationManœuvre& manœuvre = manœuvres_.back();
    CHECK_EQ(manœuvre.initial_time(), segments_.back()->last().time());
    manœuvre.set_coasting_trajectory(segments_.back());
    AddSegment();
    BurnLastSegment(manœuvre);
    AddSegment();
    CoastLastSegment(desired_final_time_);
  }
}

bool FlightPlan::RecomputeSegments() {
  // It is important that the segments be destroyed in (reverse chronological)
  // order of the forks.
  while (segments_.size() > 1) {
    PopLastSegment();
  }
  ResetLastSegment();
  for (auto& manœuvre : manœuvres_) {
    CoastLastSegment(manœuvre.initial_time());
    manœuvre.set_coasting_trajectory(segments_.back());
    AddSegment();
    BurnLastSegment(manœuvre);
    AddSegment();
  }
  CoastLastSegment(desired_final_time_);
  return anomalous_segments_ <= 2;
}

void FlightPlan::BurnLastSegment(NavigationManœuvre const& manœuvre) {
  if (anomalous_segments_ > 0) {
    return;
  } else if (manœuvre.initial_time() < manœuvre.final_time()) {
    if (manœuvre.is_inertially_fixed()) {
      bool const reached_desired_final_time =
          ephemeris_->FlowWithAdaptiveStep(segments_.back(),
                                           manœuvre.IntrinsicAcceleration(),
                                           manœuvre.final_time(),
                                           adaptive_step_parameters_,
                                           max_ephemeris_steps_per_frame,
                                           /*last_point_only=*/false).ok();
      if (!reached_desired_final_time) {
        anomalous_segments_ = 1;
      }
    } else {
      // We decompose the manœuvre in smaller unguided manœuvres (movements),
      // which are unguided.
      // TODO(egg): eventually we should just compute the direction in the
      // right-hand-side of the integrator, but for that we need the velocity,
      // which means we need general Runge-Kutta integrators.
      constexpr int movements = 100;
      Mass remaining_mass = manœuvre.initial_mass();
      serialization::DynamicFrame serialized_manœuvre_frame;
      manœuvre.frame()->WriteToMessage(&serialized_manœuvre_frame);
      for (int i = 0; i < movements; ++i) {
        NavigationManœuvre movement(manœuvre.thrust(),
                                    remaining_mass,
                                    manœuvre.specific_impulse(),
                                    manœuvre.direction(),
                                    NavigationFrame::ReadFromMessage(
                                        serialized_manœuvre_frame, ephemeris_),
                                    /*is_inertially_fixed=*/true);
        movement.set_duration(manœuvre.duration() / movements);
        movement.set_initial_time(segments_.back()->last().time());
        movement.set_coasting_trajectory(segments_.back());
        remaining_mass = movement.final_mass();
        bool const reached_desired_final_time =
            ephemeris_->FlowWithAdaptiveStep(segments_.back(),
                                             movement.IntrinsicAcceleration(),
                                             movement.final_time(),
                                             adaptive_step_parameters_,
                                             max_ephemeris_steps_per_frame,
                                             /*last_point_only=*/false).ok();
        if (!reached_desired_final_time) {
          anomalous_segments_ = 1;
          return;
        }
      }
    }
  }
}

void FlightPlan::CoastLastSegment(Instant const& desired_final_time) {
  if (anomalous_segments_ > 0) {
    return;
  } else {
    bool const reached_desired_final_time =
        ephemeris_->FlowWithAdaptiveStep(
                        segments_.back(),
                        Ephemeris<Barycentric>::NoIntrinsicAcceleration,
                        desired_final_time,
                        adaptive_step_parameters_,
                        max_ephemeris_steps_per_frame,
                        /*last_point_only=*/false).ok();
    if (!reached_desired_final_time) {
      anomalous_segments_ = 1;
    }
  }
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

void FlightPlan::AddSegment() {
  segments_.emplace_back(segments_.back()->NewForkAtLast());
  if (anomalous_segments_ > 0) {
    ++anomalous_segments_;
  }
}

void FlightPlan::ResetLastSegment() {
  segments_.back()->ForgetAfter(segments_.back()->Fork().time());
  if (anomalous_segments_ == 1) {
    // If there was one anomalous segment, it was the last one, which was
    // anomalous because it ended early.  It is no longer anomalous.
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

DiscreteTrajectory<Barycentric>& FlightPlan::last_coast() {
  return *segments_.back();
}

DiscreteTrajectory<Barycentric>& FlightPlan::penultimate_coast() {
  // The penultimate coast is the antepenultimate segment.
  return *segments_[segments_.size() - 3];
}

}  // namespace internal_flight_plan
}  // namespace ksp_plugin
}  // namespace principia
