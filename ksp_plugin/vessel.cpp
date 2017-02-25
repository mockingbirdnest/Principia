
#include "ksp_plugin/vessel.hpp"

#include <algorithm>
#include <limits>
#include <list>
#include <vector>

#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "ksp_plugin/pile_up.hpp"
#include "ksp_plugin/vessel_subsets.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/make_not_null.hpp"
#include "vessel.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_vessel {

using base::make_not_null_unique;
using geometry::BarycentreCalculator;
using geometry::Position;
using integrators::DormandElMikkawyPrince1986RKN434FM;
using integrators::McLachlanAtela1992Order5Optimal;
using quantities::IsFinite;
using quantities::Time;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Second;

Vessel::~Vessel() {
  CHECK(!is_piled_up());
}

Vessel::Vessel(not_null<Celestial const*> const parent,
               not_null<Ephemeris<Barycentric>*> const ephemeris,
               Ephemeris<Barycentric>::AdaptiveStepParameters const&
                   prediction_adaptive_step_parameters)
    : body_(),
      prediction_adaptive_step_parameters_(prediction_adaptive_step_parameters),
      parent_(parent),
      ephemeris_(ephemeris) {}

not_null<MasslessBody const*> Vessel::body() const {
  return &body_;
}

bool Vessel::is_initialized() const {
  CHECK_EQ(history_ == nullptr, prolongation_ == nullptr);
  CHECK_EQ(history_ == nullptr, prediction_ == nullptr);
  return history_ != nullptr;
}

not_null<Celestial const*> Vessel::parent() const {
  return parent_;
}

void Vessel::set_parent(not_null<Celestial const*> const parent) {
  parent_ = parent;
}

void Vessel::clear_parts() {
  parts_.clear();
}

void Vessel::add_part(not_null<Part const*> part) {
  parts.push_back(part);
}

DiscreteTrajectory<Barycentric> const& Vessel::history() const {
  CHECK(is_initialized());
  return *history_;
}

DiscreteTrajectory<Barycentric> const& Vessel::prolongation() const {
  CHECK(is_initialized());
  return *prolongation_;
}

DiscreteTrajectory<Barycentric> const& Vessel::prediction() const {
  CHECK(is_initialized());
  return *prediction_;
}

FlightPlan& Vessel::flight_plan() const {
  CHECK(has_flight_plan());
  return *flight_plan_;
}

bool Vessel::has_flight_plan() const {
  return flight_plan_ != nullptr;
}

void Vessel::AdvanceTime(Instant const& time) {
  std::vector<DiscreteTrajectory<Barycentric>::Iterator> its;
  for (auto const part : parts_) {
    its.push_back(part->tail().Begin());
  }
  bool at_end = false;
  while (!at_end) {
    Instant const time = its[0].time();
    BarycentreCalculator<DegreesOfFreedom<Barycentric>, Mass> calculator;
    for (int i = 0; i < parts_.size(); ++i) {
      auto const part = parts_[i];
      auto const it = its[i];
      CHECK_EQ(time, it.time());
      calculator.Add(it.degrees_of_freedom(), part->mass());
    }
    DegreesOfFreedom<Barycentric> const vessel_degrees_of_freedom =
        calculator.Get();
    AppendToPsychohistory(
        common_time,
        vessel_degrees_of_freedom,
        /*authoritative=*/its[0] != parts_[0]->tail().last() ||
            parts_[0]->tail_is_authoritative);
  }
  FlowProlongation(time);
}

void Vessel::ForgetBefore(Instant const& time) {
  CHECK(is_initialized());
  if (prediction_->Fork().time() < time) {
    history_->DeleteFork(prediction_);
    prediction_ = history_->NewForkAtLast();
  }
  if (flight_plan_ != nullptr) {
    flight_plan_->ForgetBefore(time, [this]() { flight_plan_.reset(); });
  }
  history_->ForgetBefore(time);
}

void Vessel::CreateFlightPlan(
    Instant const& final_time,
    Mass const& initial_mass,
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        flight_plan_adaptive_step_parameters) {
  auto const history_last = history().last();
  flight_plan_ = std::make_unique<FlightPlan>(
      initial_mass,
      /*initial_time=*/history_last.time(),
      /*initial_degrees_of_freedom=*/history_last.degrees_of_freedom(),
      final_time,
      ephemeris_,
      flight_plan_adaptive_step_parameters);
}

void Vessel::DeleteFlightPlan() {
  flight_plan_.reset();
}

void Vessel::UpdatePrediction(Instant const& last_time) {
  CHECK(is_initialized());
  history_->DeleteFork(prediction_);
  prediction_ = history_->NewForkAtLast();
  auto const prolongation_last = prolongation_->last();
  if (history_->last().time() != prolongation_last.time()) {
    prediction_->Append(prolongation_last.time(),
                        prolongation_last.degrees_of_freedom());
  }
  FlowPrediction(last_time);
}

DiscreteTrajectory<Barycentric> const& Vessel::psychohistory() const {
  return psychohistory_;
}

bool Vessel::psychohistory_is_authoritative() const {
  return psychohistory_is_authoritative_;
}

void Vessel::WriteToMessage(
    not_null<serialization::Vessel*> const message) const {
  // TODO(phl): Implement.
}

not_null<std::unique_ptr<Vessel>> Vessel::ReadFromMessage(
    serialization::Vessel const& message,
    not_null<Ephemeris<Barycentric>*> const ephemeris,
    not_null<Celestial const*> const parent) {
  // TODO(phl): Implement.
  return std::unique_ptr<Vessel>{};
}

Vessel::Vessel()
    : body_(),
      prediction_adaptive_step_parameters_(DefaultPredictionParameters()),
      parent_(testing_utilities::make_not_null<Celestial const*>()),
      ephemeris_(testing_utilities::make_not_null<Ephemeris<Barycentric>*>()) {}

void Vessel::FlowPrediction(Instant const& time) {
  if (time > prediction_->last().time()) {
    bool const finite_time = IsFinite(time - prediction_->last().time());
    Instant const t = finite_time ? time : ephemeris_->t_max();
    // This will not prolong the ephemeris if |time| is infinite (but it may do
    // so if it is finite).
    bool const reached_t = ephemeris_->FlowWithAdaptiveStep(
        prediction_,
        Ephemeris<Barycentric>::NoIntrinsicAcceleration,
        t,
        prediction_adaptive_step_parameters_,
        FlightPlan::max_ephemeris_steps_per_frame,
        /*last_point_only=*/false);
    if (!finite_time && reached_t) {
      // This will prolong the ephemeris by |max_ephemeris_steps_per_frame|.
      ephemeris_->FlowWithAdaptiveStep(
        prediction_,
        Ephemeris<Barycentric>::NoIntrinsicAcceleration,
        time,
        prediction_adaptive_step_parameters_,
        FlightPlan::max_ephemeris_steps_per_frame,
        /*last_point_only=*/false);
    }
  }
}

Ephemeris<Barycentric>::FixedStepParameters DefaultHistoryParameters() {
  return Ephemeris<Barycentric>::FixedStepParameters(
             McLachlanAtela1992Order5Optimal<Position<Barycentric>>(),
             /*step=*/10 * Second);
}

Ephemeris<Barycentric>::AdaptiveStepParameters DefaultProlongationParameters() {
  return Ephemeris<Barycentric>::AdaptiveStepParameters(
             DormandElMikkawyPrince1986RKN434FM<Position<Barycentric>>(),
             /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
             /*length_integration_tolerance=*/1 * Milli(Metre),
             /*speed_integration_tolerance=*/1 * Milli(Metre) / Second);
}

Ephemeris<Barycentric>::AdaptiveStepParameters DefaultPredictionParameters() {
  return Ephemeris<Barycentric>::AdaptiveStepParameters(
             DormandElMikkawyPrince1986RKN434FM<Position<Barycentric>>(),
             /*max_steps=*/1000,
             /*length_integration_tolerance=*/1 * Metre,
             /*speed_integration_tolerance=*/1 * Metre / Second);
}

}  // namespace internal_vessel
}  // namespace ksp_plugin
}  // namespace principia
