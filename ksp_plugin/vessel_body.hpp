
#pragma once

#include "ksp_plugin/vessel.hpp"

#include <vector>

#include "quantities/si.hpp"
#include "testing_utilities/make_not_null.hpp"

namespace principia {

using quantities::si::Kilogram;

namespace ksp_plugin {

inline Vessel::Vessel(not_null<Celestial const*> const parent)
    : body_(),
      parent_(parent) {}

inline not_null<MasslessBody const*> Vessel::body() const {
  return &body_;
}

inline bool Vessel::is_synchronized() const {
  bool const synchronized = history_ != nullptr;
  if (synchronized) {
    CHECK(owned_prolongation_ == nullptr);
  }
  return synchronized;
}

inline bool Vessel::is_initialized() const {
  bool const initialized = prolongation_ != nullptr;
  if (!initialized) {
    CHECK(owned_prolongation_ == nullptr);
  }
  return initialized;
}

inline not_null<Celestial const*> Vessel::parent() const {
  return parent_;
}

inline void Vessel::set_parent(not_null<Celestial const*> const parent) {
  parent_ = parent;
}

inline DiscreteTrajectory<Barycentric> const& Vessel::history() const {
  CHECK(is_synchronized());
  return *history_;
}

inline not_null<DiscreteTrajectory<Barycentric>*> Vessel::mutable_history() {
  CHECK(is_synchronized());
  return history_.get();
}

inline DiscreteTrajectory<Barycentric> const& Vessel::prolongation() const {
  CHECK(is_initialized());
  return *prolongation_;
}

inline not_null<DiscreteTrajectory<Barycentric>*>
Vessel::mutable_prolongation() {
  CHECK(is_initialized());
  return prolongation_;
}

inline not_null<FlightPlan*> Vessel::flight_plan() const {
  CHECK(is_initialized());
  return flight_plan_.get();
}

inline bool Vessel::has_flight_plan() const {
  return flight_plan_ != nullptr;
}

inline DiscreteTrajectory<Barycentric> const& Vessel::prediction() const {
  CHECK(has_prediction());
  return *prediction_;
}

inline bool Vessel::has_prediction() const {
  return prediction_ != nullptr;
}

inline void Vessel::CreateProlongation(
    Instant const& time,
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
  CHECK(!is_synchronized());
  CHECK(!is_initialized());
  CHECK(owned_prolongation_ == nullptr);
  owned_prolongation_ = std::make_unique<DiscreteTrajectory<Barycentric>>();
  owned_prolongation_->Append(time, degrees_of_freedom);
  prolongation_ = owned_prolongation_.get();
}

inline void Vessel::CreateHistoryAndForkProlongation(
    Instant const& time,
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
  CHECK(!is_synchronized());
  history_ = std::make_unique<DiscreteTrajectory<Barycentric>>();
  history_->Append(time, degrees_of_freedom);
  prolongation_ = history_->NewForkAtLast();
  owned_prolongation_.reset();
}

inline void Vessel::ResetProlongation(Instant const& time) {
  CHECK(is_initialized());
  CHECK(is_synchronized());
  CHECK(owned_prolongation_ == nullptr);
  history_->DeleteFork(&prolongation_);
  prolongation_ = history_->NewForkWithCopy(time);
}

inline void Vessel::CreateFlightPlan(
    Instant const& final_time,
    Mass const& initial_mass,
    not_null<Ephemeris<Barycentric>*> ephemeris,
    AdaptiveStepSizeIntegrator<
        Ephemeris<Barycentric>::NewtonianMotionEquation> const& integrator,
    Length const& length_integration_tolerance,
    Speed const& speed_integration_tolerance) {
  if (!is_synchronized()) {
    return;
  }
  flight_plan_ = std::make_unique<FlightPlan>(
                     mutable_history(),
                     /*initial_time=*/history().last().time(),
                     /*final_time=*/final_time,
                     initial_mass,
                     ephemeris,
                     Ephemeris<Barycentric>::AdaptiveStepParameters(
                         integrator,
                         length_integration_tolerance,
                         speed_integration_tolerance));
}

inline void Vessel::DeleteFlightPlan() {
  flight_plan_.reset();
}

inline void Vessel::UpdatePrediction(
    not_null<Ephemeris<Barycentric>*> ephemeris,
    AdaptiveStepSizeIntegrator<
        Ephemeris<Barycentric>::NewtonianMotionEquation> const& integrator,
    Instant const& last_time,
    Length const& prediction_length_tolerance,
    Speed const& prediction_speed_tolerance) {
  if (!is_synchronized()) {
    return;
  }
  DeletePrediction();
  prediction_ = mutable_history()->NewForkAtLast();
  if (history().last().time() != prolongation().last().time()) {
    prediction_->Append(prolongation().last().time(),
                        prolongation().last().degrees_of_freedom());
  }
  ephemeris->FlowWithAdaptiveStep(
      prediction_,
      Ephemeris<Barycentric>::kNoIntrinsicAcceleration,
      last_time,
      Ephemeris<Barycentric>::AdaptiveStepParameters(
          integrator,
          prediction_length_tolerance,
          prediction_speed_tolerance));
}

inline void Vessel::DeletePrediction() {
  if (has_prediction()) {
    mutable_history()->DeleteFork(&prediction_);
  }
}

inline void Vessel::WriteToMessage(
    not_null<serialization::Vessel*> const message) const {
  CHECK(is_initialized());
  body_.WriteToMessage(message->mutable_body());
  if (is_synchronized()) {
    history_->WriteToMessage(
        message->mutable_history_and_prolongation()->mutable_history());
    prolongation_->WritePointerToMessage(
        message->mutable_history_and_prolongation()->mutable_prolongation());
  } else {
    owned_prolongation_->WriteToMessage(message->mutable_owned_prolongation());
  }
  if (prediction_ != nullptr) {
    prediction_->WritePointerToMessage(message->mutable_prediction());
  }
  if (flight_plan_ != nullptr) {
    flight_plan_->WriteToMessage(message->mutable_flight_plan());
  }
}

inline not_null<std::unique_ptr<Vessel>> Vessel::ReadFromMessage(
    serialization::Vessel const& message,
    not_null<Ephemeris<Barycentric>*> const ephemeris,
    not_null<Celestial const*> const parent) {
  auto vessel = make_not_null_unique<Vessel>(parent);
  // NOTE(egg): for now we do not read the |MasslessBody| as it can contain no
  // information.
  if (message.has_history_and_prolongation()) {
    vessel->history_ =
        DiscreteTrajectory<Barycentric>::ReadFromMessage(
            message.history_and_prolongation().history());
    vessel->prolongation_ =
        DiscreteTrajectory<Barycentric>::ReadPointerFromMessage(
            message.history_and_prolongation().prolongation(),
            vessel->history_.get());
    if (message.has_prediction()) {
      vessel->prediction_ =
          DiscreteTrajectory<Barycentric>::ReadPointerFromMessage(
              message.prediction(),
              vessel->history_.get());
    }
    if (message.has_flight_plan()) {
      vessel->flight_plan_ = FlightPlan::ReadFromMessage(
          message.flight_plan(), vessel->history_.get(), ephemeris);
    }
  } else if (message.has_owned_prolongation()) {
    vessel->owned_prolongation_ =
        DiscreteTrajectory<Barycentric>::ReadFromMessage(
            message.owned_prolongation());
    vessel->prolongation_ = vessel->owned_prolongation_.get();
    CHECK(!message.has_prediction());
    CHECK(!message.has_flight_plan());
  } else {
    LOG(FATAL) << "message does not represent an initialized Vessel";
    base::noreturn();
  }
  return vessel;
}

inline Vessel::Vessel()
    : body_(),
      parent_(testing_utilities::make_not_null<Celestial const*>()) {}

}  // namespace ksp_plugin
}  // namespace principia
