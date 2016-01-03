#include "ksp_plugin/flight_plan.hpp"

namespace principia {
namespace ksp_plugin {

FlightPlan::FlightPlan(
    not_null<DiscreteTrajectory<Barycentric>*> root,
    Instant const& initial_time,
    Instant const& final_time,
    Mass const& initial_mass,
    not_null<Ephemeris<Barycentric>*> ephemeris,
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
  //ephemeris_->Prolong(final_time_);
  auto it = root->LowerBound(initial_time);
  if (it.time() != initial_time_) {
    --it;
  }
  segments_.emplace(root->NewForkWithCopy(it.time()));
  CoastLastSegment(final_time_);
}

int FlightPlan::size() const {
  return manœuvres_.size();
}

NavigationManœuvre const& FlightPlan::Get(int index) {
  CHECK_LE(0, index);
  CHECK_LT(index, size());
  return manœuvres_[index];
}

bool FlightPlan::Append(Burn burn) {
  auto manœuvre =
      MakeManœuvre(
          std::move(burn),
          manœuvres_.empty() ? initial_mass_ : manœuvres_.back().final_mass());
  if (manœuvre.FitsBetween(
          manœuvres_.empty() ? initial_time_ : manœuvres_.back().final_time(),
          final_time_)) {
    Append(std::move(manœuvre));
    return true;
  } else {
    return false;
  }
}

void FlightPlan::RemoveLast() {
  CHECK(!manœuvres_.empty());
  manœuvres_.pop_back();
  segments_.pop();  // Last coast.
  segments_.pop();  // Last burn.
  ResetLastSegment();
  CoastLastSegment(final_time_);
}

bool FlightPlan::ReplaceLast(Burn burn) {
  CHECK(!manœuvres_.empty());
  auto manœuvre =
      MakeManœuvre(std::move(burn), manœuvres_.back().initial_mass());
  if (manœuvre.FitsBetween(manœuvres_.size() == 1
                               ? initial_time_
                               : manœuvres_[manœuvres_.size() - 2].final_time(),
                           final_time_)) {
    manœuvres_.pop_back();
    segments_.pop();  // Last coast.
    segments_.pop();  // Last burn.
    Append(std::move(manœuvre));
    return true;
  } else {
    return false;
  }
}

bool FlightPlan::SetFinalTime(Instant const& final_time) {
  if (!manœuvres_.empty() && manœuvres_.back().final_time() > final_time ||
      initial_time_ > final_time) {
    return false;
  } else {
    final_time_ = final_time;
    //ephemeris_->Prolong(final_time_);
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

void FlightPlan::Append(NavigationManœuvre manœuvre) {
  manœuvres_.emplace_back(std::move(manœuvre));
  {
    // Hide the moved-from |manœuvre|.
    NavigationManœuvre const& manœuvre = manœuvres_.back();
    ResetLastSegment();
    CoastLastSegment(manœuvre.initial_time());
    AddSegment();
    BurnLastSegment(manœuvre);
    AddSegment();
    CoastLastSegment(final_time_);
  }
}

void FlightPlan::RecomputeSegments() {
  // It is important that the segments be destroyed in (reverse chronological)
  // order, since the destructor of each segment references the previous
  // segment.
  while(segments_.size() > 1) {
    segments_.pop();
  }
  ResetLastSegment();
  for (auto const& manœuvre : manœuvres_) {
    CoastLastSegment(manœuvre.initial_time());
    AddSegment();
    BurnLastSegment(manœuvre);
    AddSegment();
  }
  CoastLastSegment(final_time_);
}

void FlightPlan::BurnLastSegment(NavigationManœuvre const& manœuvre) {
  ephemeris_->FlowWithAdaptiveStep(
      segments_.top().get(),
      manœuvre.acceleration(*segments_.top()),
      length_integration_tolerance_,
      speed_integration_tolerance_,
      integrator_,
      manœuvre.final_time());
}

void FlightPlan::CoastLastSegment(Instant const& final_time) {
  ephemeris_->FlowWithAdaptiveStep(
      segments_.top().get(),
      Ephemeris<Barycentric>::kNoIntrinsicAcceleration,
      length_integration_tolerance_,
      speed_integration_tolerance_,
      integrator_,
      final_time);
}

void FlightPlan::AddSegment() {
  segments_.emplace(segments_.top()->NewForkAtLast());
}

void FlightPlan::ResetLastSegment() {
  segments_.top()->ForgetAfter(segments_.top()->Fork().time());
}

NavigationManœuvre FlightPlan::MakeManœuvre(Burn burn,
                                            Mass const& initial_mass) {
  NavigationManœuvre manœuvre(burn.thrust,
      initial_mass,
      burn.specific_impulse,
      Normalize(burn.Δv),
      std::move(burn.frame));
  manœuvre.set_initial_time(burn.initial_time);
  manœuvre.set_Δv(burn.Δv.Norm());
  return std::move(manœuvre);
}

}  // namespace ksp_plugin
}  // namespace principia
