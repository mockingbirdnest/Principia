
#include "ksp_plugin/flight_plan.hpp"

#include "testing_utilities/make_not_null.hpp"

namespace principia {
namespace ksp_plugin {

FlightPlan::FlightPlan(
    not_null<DiscreteTrajectory<Barycentric>*> const root,
    Instant const& initial_time,
    Instant const& final_time,
    Mass const& initial_mass,
    not_null<Ephemeris<Barycentric>*> const ephemeris,
    AdaptiveStepSizeIntegrator<
        Ephemeris<Barycentric>::NewtonianMotionEquation> const& integrator,
    Length const& length_integration_tolerance,
    Speed const& speed_integration_tolerance)
    : initial_time_(initial_time),
      final_time_(final_time),
      initial_mass_(initial_mass),
      ephemeris_(ephemeris),
      integrator_(integrator),
      length_integration_tolerance_(length_integration_tolerance),
      speed_integration_tolerance_(speed_integration_tolerance) {
  CHECK(final_time_ >= initial_time_);
  auto it = root->LowerBound(initial_time_);
  if (it.time() != initial_time_) {
    --it;
    initial_time_ = it.time();
  }

  segments_.emplace_back(make_not_null_unique<DiscreteTrajectory<Barycentric>>());
  segments_.back()->Append(it.time(), it.degrees_of_freedom());

  CoastLastSegment(final_time_);
}

Instant FlightPlan::initial_time() const {
  return initial_time_;
}

Instant FlightPlan::final_time() const {
  return final_time_;
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
  if (manœuvre.FitsBetween(start_of_last_coast(), final_time_)) {
    auto recomputed_last_coast =
        CoastIfReachesManœuvreInitialTime(last_coast(), manœuvre);
    if (recomputed_last_coast != nullptr) {
      ReplaceLastSegment(std::move(recomputed_last_coast));
      Append(std::move(manœuvre));
      return true;
    }
  }
  return false;
}

void FlightPlan::RemoveLast() {
  CHECK(!manœuvres_.empty());
  manœuvres_.pop_back();
  PopLastSegment();  // Last coast.
  PopLastSegment();  // Last burn.
  ResetLastSegment();
  CoastLastSegment(final_time_);
}

bool FlightPlan::ReplaceLast(Burn burn) {
  CHECK(!manœuvres_.empty());
  auto manœuvre =
      MakeNavigationManœuvre(std::move(burn), manœuvres_.back().initial_mass());
  if (manœuvre.FitsBetween(start_of_penultimate_coast(), final_time_)) {
    auto recomputed_penultimate_coast =
        CoastIfReachesManœuvreInitialTime(penultimate_coast(), manœuvre);
    if (recomputed_penultimate_coast != nullptr) {
      manœuvres_.pop_back();
      PopLastSegment();  // Last coast.
      PopLastSegment();  // Last burn.
      ReplaceLastSegment(std::move(recomputed_penultimate_coast));
      Append(std::move(manœuvre));
      return true;
    }
  }
  return false;
}

bool FlightPlan::SetFinalTime(Instant const& final_time) {
  if (start_of_last_coast() > final_time) {
    return false;
  } else {
    final_time_ = final_time;
    ResetLastSegment();
    CoastLastSegment(final_time_);
    return true;
  }
}

void FlightPlan::SetTolerances(
    Length const& length_integration_tolerance,
    Speed const& speed_integration_tolerance) {
  length_integration_tolerance_ = length_integration_tolerance;
  speed_integration_tolerance_ = speed_integration_tolerance;
  RecomputeSegments();
}

int FlightPlan::number_of_segments() const {
  return segments_.size();
}

void FlightPlan::GetSegment(
    int const index,
    not_null<DiscreteTrajectory<Barycentric>::Iterator*> begin,
    not_null<DiscreteTrajectory<Barycentric>::Iterator*> end) const {
  CHECK_LE(0, index);
  CHECK_LT(index, number_of_segments());
  *begin = segments_[index]->Begin();
  *end = segments_[index]->End();
}

void FlightPlan::WriteToMessage(
    not_null<serialization::FlightPlan*> const message) const {
  initial_mass_.WriteToMessage(message->mutable_initial_mass());
  initial_time_.WriteToMessage(message->mutable_initial_time());
  final_time_.WriteToMessage(message->mutable_final_time());
  length_integration_tolerance_.WriteToMessage(
      message->mutable_length_integration_tolerance());
  speed_integration_tolerance_.WriteToMessage(
      message->mutable_speed_integration_tolerance());
  integrator_.WriteToMessage(message->mutable_integrator());
  for (auto const& manœuvre : manœuvres_) {
    manœuvre.WriteToMessage(message->add_manoeuvre());
  }
}

std::unique_ptr<FlightPlan> FlightPlan::ReadFromMessage(
    serialization::FlightPlan const& message,
    not_null<DiscreteTrajectory<Barycentric>*> const root,
    not_null<Ephemeris<Barycentric>*> const ephemeris) {
  auto flight_plan = std::make_unique<FlightPlan>(
      root,
      Instant::ReadFromMessage(message.initial_time()),
      Instant::ReadFromMessage(message.final_time()),
      Mass::ReadFromMessage(message.initial_mass()),
      ephemeris,
      AdaptiveStepSizeIntegrator<
          Ephemeris<Barycentric>::NewtonianMotionEquation>::
          ReadFromMessage(message.integrator()),
      Length::ReadFromMessage(message.length_integration_tolerance()),
      Speed::ReadFromMessage(message.speed_integration_tolerance()));
  for (int i = 0; i < message.manoeuvre_size(); ++i) {
    flight_plan->manœuvres_.emplace_back(
        NavigationManœuvre::ReadFromMessage(message.manoeuvre(i), ephemeris));
  }
  flight_plan->RecomputeSegments();
  return std::move(flight_plan);
}

FlightPlan::FlightPlan()
    : ephemeris_(testing_utilities::make_not_null<Ephemeris<Barycentric>*>()),
      integrator_(
          *testing_utilities::make_not_null<
              AdaptiveStepSizeIntegrator<
                  Ephemeris<Barycentric>::NewtonianMotionEquation>*>()) {}

void FlightPlan::Append(NavigationManœuvre manœuvre) {
  manœuvres_.emplace_back(std::move(manœuvre));
  {
    // Hide the moved-from |manœuvre|.
    NavigationManœuvre& manœuvre = manœuvres_.back();
    CHECK_EQ(manœuvre.initial_time(), segments_.back()->last().time());
    manœuvre.set_coasting_trajectory(segments_.back().get());
    AddSegment();
    BurnLastSegment(manœuvre);
    AddSegment();
    CoastLastSegment(final_time_);
  }
}

void FlightPlan::RecomputeSegments() {
  // It is important that the segments be destroyed in (reverse chronological)
  // order of the forks.
  while (segments_.size() > 1) {
    PopLastSegment();
  }
  ResetLastSegment();
  for (auto& manœuvre : manœuvres_) {
    CoastLastSegment(manœuvre.initial_time());
    manœuvre.set_coasting_trajectory(segments_.back().get());
    AddSegment();
    BurnLastSegment(manœuvre);
    AddSegment();
  }
  CoastLastSegment(final_time_);
}

void FlightPlan::BurnLastSegment(NavigationManœuvre const& manœuvre) {
  if (anomalous_segments_ > 0) {
    return;
  } else if (manœuvre.initial_time() < manœuvre.final_time()) {
    bool const reached_final_time =
        ephemeris_->FlowWithAdaptiveStep(segments_.back().get(),
                                         manœuvre.IntrinsicAcceleration(),
                                         length_integration_tolerance_,
                                         speed_integration_tolerance_,
                                         integrator_,
                                         manœuvre.final_time());
    if (!reached_final_time) {
      ++anomalous_segments_;
    }
  }
}

void FlightPlan::CoastLastSegment(Instant const& final_time) {
  if (anomalous_segments_ > 0) {
    return;
  } else {
    bool const reached_final_time =
          ephemeris_->FlowWithAdaptiveStep(
                          segments_.back().get(),
                          Ephemeris<Barycentric>::kNoIntrinsicAcceleration,
                          length_integration_tolerance_,
                          speed_integration_tolerance_,
                          integrator_,
                          final_time);
    if (!reached_final_time) {
      ++anomalous_segments_;
    }
  }
}

void FlightPlan::ReplaceLastSegment(
    not_null<std::unique_ptr<DiscreteTrajectory<Barycentric>>> segment) {
  CHECK_EQ(segment->Begin().time(), segments_.back()->Begin().time());
  PopLastSegment();
  // |segment| must not be anomalous, so it cannot not follow an anomalous
  // segment.
  CHECK_EQ(0, anomalous_segments_);
  segments_.emplace_back(std::move(segment));
}

void FlightPlan::AddSegment() {
  auto const penultimate_last = segments_.back()->last();
  segments_.emplace_back(
      make_not_null_unique<DiscreteTrajectory<Barycentric>>());
  segments_.back()->Append(penultimate_last.time(),
                           penultimate_last.degrees_of_freedom());
  if (anomalous_segments_ > 0) {
    ++anomalous_segments_;
  }
}

void FlightPlan::ResetLastSegment() {
  segments_.back()->ForgetAfter(segments_.back()->Begin().time());
  if (anomalous_segments_ == 1) {
    // If there was one anomalous segment, it was the last one, which was
    // anomalous because it ended early.  It is no longer anomalous.
    anomalous_segments_ = 0;
  }
}

void FlightPlan::PopLastSegment() {
  segments_.pop_back();
  if (anomalous_segments_ > 0) {
    --anomalous_segments_;
  }
}

std::unique_ptr<DiscreteTrajectory<Barycentric>>
FlightPlan::CoastIfReachesManœuvreInitialTime(
    DiscreteTrajectory<Barycentric> const& coast,
    NavigationManœuvre const& manœuvre) const {
  auto recomputed_coast = std::make_unique<DiscreteTrajectory<Barycentric>>();
  recomputed_coast->Append(coast.Begin().time(),
                           coast.Begin().degrees_of_freedom());
  bool const reached_manœuvre_initial_time =
      ephemeris_->FlowWithAdaptiveStep(
          recomputed_coast.get(),
          Ephemeris<Barycentric>::kNoIntrinsicAcceleration,
          length_integration_tolerance_,
          speed_integration_tolerance_,
          integrator_,
          manœuvre.initial_time());
  if (!reached_manœuvre_initial_time) {
    recomputed_coast.reset();
  }
  return std::move(recomputed_coast);
}

Instant FlightPlan::start_of_last_coast() const {
  return manœuvres_.empty() ? initial_time_ : manœuvres_.back().final_time();
}

Instant FlightPlan::start_of_penultimate_coast() const {
  return manœuvres_.size() == 1
             ? initial_time_
             : manœuvres_[manœuvres_.size() - 2].final_time();
}

DiscreteTrajectory<Barycentric> const& FlightPlan::last_coast() const {
  return *segments_.back();
}

DiscreteTrajectory<Barycentric> const& FlightPlan::penultimate_coast() const {
  // The penultimate coast is the antepenultimate segment.
  return *segments_[segments_.size() - 3];
}

}  // namespace ksp_plugin
}  // namespace principia
