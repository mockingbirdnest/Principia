
#include "ksp_plugin/vessel.hpp"

#include <algorithm>
#include <list>
#include <vector>

#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "ksp_plugin/pile_up.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/make_not_null.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_vessel {

using base::FindOrDie;
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

Vessel::Vessel(not_null<Celestial const*> const parent,
               not_null<Ephemeris<Barycentric>*> const ephemeris,
               Ephemeris<Barycentric>::AdaptiveStepParameters const&
                   prediction_adaptive_step_parameters)
    : body_(),
      prediction_adaptive_step_parameters_(prediction_adaptive_step_parameters),
      parent_(parent),
      ephemeris_(ephemeris),
      prediction_(make_not_null_unique<DiscreteTrajectory<Barycentric>>()) {}

not_null<MasslessBody const*> Vessel::body() const {
  return &body_;
}

not_null<Celestial const*> Vessel::parent() const {
  return parent_;
}

void Vessel::set_parent(not_null<Celestial const*> const parent) {
  parent_ = parent;
}

void Vessel::AddPart(not_null<std::unique_ptr<Part>> part) {
  PartId const id = part->part_id();
  parts_.emplace(id, std::move(part));
  kept_parts_.insert(part.get());
}

not_null<std::unique_ptr<Part>> Vessel::ExtractPart(PartId const id) {
  auto const it = parts_.find(id);
  CHECK(it != parts_.end()) << id;
  auto result = std::move(it->second);
  parts_.erase(it);
  return result;
}

void Vessel::KeepPart(PartId const id) {
  kept_parts_.insert(FindOrDie(parts_, id).get());
}

void Vessel::FreeParts() {
  dummy_part_.reset();
  for (auto it = parts_.begin(); it != parts_.end();) {
    not_null<Part const*> part = it->second.get();
    if (kept_parts_.count(part) == 0) {
      it = parts_.erase(it);
    } else {
      ++it;
    }
  }
  CHECK(!parts_.empty());
  kept_parts_.clear();
}

not_null<Part*> Vessel::part(PartId const id) const {
  return FindOrDie(parts_, id).get();
}

void Vessel::ForSomePart(std::function<void(Part&)> action) const {
  CHECK(!parts_.empty());
  action(*parts_.begin()->second);
}

void Vessel::ForAllParts(std::function<void(Part&)> action) const {
  for (auto const& pair : parts_) {
    action(*pair.second);
  }
}

DiscreteTrajectory<Barycentric> const& Vessel::prediction() const {
  return *prediction_;
}

void Vessel::set_prediction_adaptive_step_parameters(
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        prediction_adaptive_step_parameters) {
  prediction_adaptive_step_parameters_ = prediction_adaptive_step_parameters;
}

Ephemeris<Barycentric>::AdaptiveStepParameters const&
Vessel::prediction_adaptive_step_parameters() const {
  return prediction_adaptive_step_parameters_;
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
  for (auto const& pair : parts_) {
    Part const& part = *pair.second;
    its.push_back(part.tail().Begin());
  }
  for (;;) {
    Part const& first_part = *parts_.begin()->second;
    Instant const time = its[0].time();
    bool const at_end_of_tail = its[0] == first_part.tail().last();
    bool const tail_is_authoritative = first_part.tail_is_authoritative();

    BarycentreCalculator<DegreesOfFreedom<Barycentric>, Mass> calculator;
    int i = 0;
    for (auto const& pair : parts_) {
      Part& part = *pair.second;
      auto& it = its[i];
      calculator.Add(it.degrees_of_freedom(), part.mass());
      CHECK_EQ(time, it.time());
      CHECK_EQ(at_end_of_tail, it == part.tail().last());
      CHECK_EQ(tail_is_authoritative, part.tail_is_authoritative());
      if (at_end_of_tail) {
        part.tail().ForgetBefore(astronomy::InfiniteFuture);
      } else {
        ++it;
      }
      ++i;
    }
    DegreesOfFreedom<Barycentric> const vessel_degrees_of_freedom =
        calculator.Get();
    AppendToPsychohistory(
        time,
        vessel_degrees_of_freedom,
        /*authoritative=*/!at_end_of_tail || tail_is_authoritative);

    if (at_end_of_tail) {
      return;
    }
  }
}

void Vessel::ForgetBefore(Instant const& time) {
  prediction_->ForgetBefore(time);
  if (flight_plan_ != nullptr) {
    flight_plan_->ForgetBefore(time, [this]() { flight_plan_.reset(); });
  }
}

void Vessel::CreateFlightPlan(
    Instant const& final_time,
    Mass const& initial_mass,
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        flight_plan_adaptive_step_parameters) {
  auto const last = last_authoritative();
  flight_plan_ = std::make_unique<FlightPlan>(
      initial_mass,
      /*initial_time=*/last.time(),
      /*initial_degrees_of_freedom=*/last.degrees_of_freedom(),
      final_time,
      ephemeris_,
      flight_plan_adaptive_step_parameters);
}

void Vessel::DeleteFlightPlan() {
  flight_plan_.reset();
}

void Vessel::UpdatePrediction(Instant const& last_time) {
  prediction_ = make_not_null_unique<DiscreteTrajectory<Barycentric>>();
  auto const last = psychohistory_.last();
  prediction_->Append(last.time(), last.degrees_of_freedom());
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
      ephemeris_(testing_utilities::make_not_null<Ephemeris<Barycentric>*>()),
      prediction_(make_not_null_unique<DiscreteTrajectory<Barycentric>>()) {}

void Vessel::AppendToPsychohistory(
    Instant const& time,
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom,
    bool const authoritative) {
  if (!psychohistory_is_authoritative_) {
    psychohistory_.ForgetAfter((--psychohistory_.last()).time());
  }
  psychohistory_.Append(time, degrees_of_freedom);
  psychohistory_is_authoritative_ = authoritative;
}

void Vessel::FlowPrediction(Instant const& time) {
  if (time > prediction_->last().time()) {
    bool const finite_time = IsFinite(time - prediction_->last().time());
    Instant const t = finite_time ? time : ephemeris_->t_max();
    // This will not prolong the ephemeris if |time| is infinite (but it may do
    // so if it is finite).
    bool const reached_t = ephemeris_->FlowWithAdaptiveStep(
        prediction_.get(),
        Ephemeris<Barycentric>::NoIntrinsicAcceleration,
        t,
        prediction_adaptive_step_parameters_,
        FlightPlan::max_ephemeris_steps_per_frame,
        /*last_point_only=*/false);
    if (!finite_time && reached_t) {
      // This will prolong the ephemeris by |max_ephemeris_steps_per_frame|.
      ephemeris_->FlowWithAdaptiveStep(
        prediction_.get(),
        Ephemeris<Barycentric>::NoIntrinsicAcceleration,
        time,
        prediction_adaptive_step_parameters_,
        FlightPlan::max_ephemeris_steps_per_frame,
        /*last_point_only=*/false);
    }
  }
}

DiscreteTrajectory<Barycentric>::Iterator Vessel::last_authoritative() const {
  auto it = psychohistory_.last();
  if (!psychohistory_is_authoritative_) {
    --it;
  }
  return it;
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
