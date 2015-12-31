#include "ksp_plugin/flight_plan.hpp"

namespace principia {
namespace ksp_plugin {

FlightPlan::FlightPlan(
    not_null<DiscreteTrajectory<Barycentric>*> root,
    Instant const& initial_time,
    Mass const& initial_mass,
    not_null<Ephemeris<Barycentric>*> ephemeris,
    AdaptiveStepSizeIntegrator<
        Ephemeris<Barycentric>::NewtonianMotionEquation> const& integrator)
    : initial_time_(initial_time),
      initial_mass_(initial_mass),
      ephemeris_(ephemeris),
      integrator_(integrator) {
  segments_.emplace(
      root->NewForkWithCopy(root->LowerBound(initial_time).time()));
}

int FlightPlan::size() const {
  return burns_.size();
}

Manœuvre<Barycentric, Navigation> const& FlightPlan::Get(int index) {
  CHECK_LE(0, index);
  CHECK_LT(index, size());
  return manœuvres_[index];
}

bool FlightPlan::Append(BurnDefinition burn) {
  auto manœuvre = 
    Manœuvre<Barycentric, Navigation>(
        burn.thrust,
        manœuvres_.empty() ? initial_mass_ : manœuvres_.back().final_mass(),
        burn.specific_impulse,
        Normalize(burn.Δv),
        burn.frame.get());
  manœuvre.set_initial_time(burn.initial_time);
  manœuvre.set_Δv(burn.Δv.Norm());
  if (manœuvre.FitsBetween(
          manœuvres_.empty() ? initial_time_ : manœuvres_.back().final_time(),
          final_time_)) {
    return false;
  } else {
    burns_.emplace(std::move(burn));
    {
      // Hide the moved-from |burn|.
      BurnDefinition& burn = burns_.top();
      manœuvres_.emplace_back(manœuvre);
      // Reset the last coast.
      segments_.top()->ForgetBefore(segments_.top()->Fork().time());
      // Prolong the last coast until the start of the new burn.
      ephemeris_->FlowWithAdaptiveStep(
          segments_.top().get(),
          ephemeris_->kNoIntrinsicAcceleration,
          length_integration_tolerance_,
          speed_integration_tolerance_,
          integrator_,
          burn.initial_time);
      // Create the burn segment.
      segments_.emplace(segments_.top()->NewForkAtLast());
      ephemeris_->FlowWithAdaptiveStep(
          segments_.top().get(),
          manœuvres_.back().acceleration(*segments_.top()),
          length_integration_tolerance_,
          speed_integration_tolerance_,
          integrator_,
          manœuvres_.back().final_time());
      // Create the final coast segment.
      segments_.emplace(segments_.top()->NewForkAtLast());
      ephemeris_->FlowWithAdaptiveStep(
          segments_.top().get(),
          ephemeris_->kNoIntrinsicAcceleration,
          length_integration_tolerance_,
          speed_integration_tolerance_,
          integrator_,
          final_time_);
      return true;
    }
  }
}

void FlightPlan::RemoveLast() {
  CHECK(!burns_.empty());
  burns_.pop();
  manœuvres_.pop_back();
  segments_.pop();  // Last coast.
  segments_.pop();  // Last burn.
  // Reset the last remaining coast.
  segments_.top()->ForgetBefore(segments_.top()->Fork().time());
  // Prolong the last remaining coast until |final_time_|.
  ephemeris_->FlowWithAdaptiveStep(
      segments_.top().get(),
      ephemeris_->kNoIntrinsicAcceleration,
      length_integration_tolerance_,
      speed_integration_tolerance_,
      integrator_,
      final_time_);
}

bool FlightPlan::ReplaceLast(BurnDefinition burn) {
  CHECK(!burns_.empty());

  auto manœuvre =
    Manœuvre<Barycentric, Navigation>(burn.thrust,
                                      manœuvres_.back().initial_mass(),
                                      burn.specific_impulse,
                                      Normalize(burn.Δv),
                                      burn.frame.get());
  manœuvre.set_initial_time(burn.initial_time);
  manœuvre.set_Δv(burn.Δv.Norm());
  if (manœuvre.FitsBetween(manœuvres_.size() == 1
                               ? initial_time_
                               : manœuvres_[manœuvres_.size() - 2].final_time(),
                           final_time_)) {
    return false;
  } else {
    burns_.pop();
    manœuvres_.pop_back();
    segments_.pop();  // Last coast.
    segments_.pop();  // Last burn.
    CHECK(Append(std::move(burn)));
    return true;
  }
}

bool FlightPlan::set_final_time(Instant const& final_time) {
  if (!manœuvres_.empty() && manœuvres_.back().final_time() > final_time ||
      initial_time_ > final_time) {
    return false;
  } else {
    final_time_ = final_time;
  segments_.top()->ForgetBefore(segments_.top()->Fork().time());
  ephemeris_->FlowWithAdaptiveStep(
      segments_.top().get(),
      ephemeris_->kNoIntrinsicAcceleration,
      length_integration_tolerance_,
      speed_integration_tolerance_,
      integrator_,
      final_time_);
    return true;
  }
}

}  // namespace ksp_plugin
}  // namespace principia
