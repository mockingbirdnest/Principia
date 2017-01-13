
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
               Ephemeris<Barycentric>::FixedStepParameters const&
                   history_fixed_step_parameters,
               Ephemeris<Barycentric>::AdaptiveStepParameters const&
                   prolongation_adaptive_step_parameters,
               Ephemeris<Barycentric>::AdaptiveStepParameters const&
                   prediction_adaptive_step_parameters)
    : body_(),
      history_fixed_step_parameters_(history_fixed_step_parameters),
      prolongation_adaptive_step_parameters_(
          prolongation_adaptive_step_parameters),
      prediction_adaptive_step_parameters_(prediction_adaptive_step_parameters),
      parent_(parent),
      ephemeris_(ephemeris),
      subset_node_(make_not_null_unique<Subset<Vessel>::Node>()) {}

not_null<MasslessBody const*> Vessel::body() const {
  return &body_;
}

bool Vessel::is_initialized() const {
  CHECK_EQ(history_ == nullptr, prediction_ == nullptr);
  return history_ != nullptr;
}

not_null<Celestial const*> Vessel::parent() const {
  return parent_;
}

void Vessel::set_parent(not_null<Celestial const*> const parent) {
  parent_ = parent;
}

DiscreteTrajectory<Barycentric> const& Vessel::history() const {
  CHECK(is_initialized());
  return *history_;
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

void Vessel::set_dirty() {
  is_dirty_ = true;
}

bool Vessel::is_dirty() const {
  return is_dirty_;
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

void Vessel::CreateHistoryAndForkProlongation(
    Instant const& time,
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
  CHECK(!is_initialized());
  history_ = std::make_unique<DiscreteTrajectory<Barycentric>>();
  history_->Append(time, degrees_of_freedom);
  prolongation_ = history_->NewForkAtLast();
  prediction_ = history_->NewForkAtLast();
}

void Vessel::AdvanceTimeNotInBubble(Instant const& time) {
  CHECK(is_initialized());
  AdvanceHistoryIfNeeded(time);
  FlowProlongation(time);
}

void Vessel::AdvanceTimeInBubble(
    Instant const& time,
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom) {
  CHECK(is_initialized());
  AdvanceHistoryIfNeeded(time);
  prolongation_->Append(time, degrees_of_freedom);
  is_dirty_ = true;
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

void Vessel::clear_mass() {
  mass_ = Mass();
}

void Vessel::increment_mass(Mass const& mass) {
  mass_ += mass;
}

Mass const & Vessel::mass() const {
  return mass_;
}

void Vessel::clear_intrinsic_force() {
  intrinsic_force_ = Vector<Force, Barycentric>();
}

void Vessel::increment_intrinsic_force(
    Vector<Force, Barycentric> const& intrinsic_force) {
  intrinsic_force_ += intrinsic_force;
}

Vector<Force, Barycentric> const & Vessel::intrinsic_force() const {
  return intrinsic_force_;
}

void Vessel::set_containing_pile_up(IteratorOn<std::list<PileUp>> pile_up) {
  CHECK(!is_piled_up());
  containing_pile_up_ = pile_up;
}

std::experimental::optional<IteratorOn<std::list<PileUp>>>
Vessel::containing_pile_up() const {
  return containing_pile_up_;
}

bool Vessel::is_piled_up() const {
  // TODO(egg): |has_value()| once we have a standard |optional|.
  return static_cast<bool>(containing_pile_up_);
}

void Vessel::clear_pile_up() {
  if (is_piled_up()) {
    IteratorOn<std::list<PileUp>> pile_up = *containing_pile_up_;
    for (not_null<Vessel*> const vessel : pile_up.iterator()->vessels()) {
      vessel->containing_pile_up_ = std::experimental::nullopt;
    }
    CHECK(!is_piled_up());
    pile_up.Erase();
  }
}

void Vessel::WriteToMessage(
    not_null<serialization::Vessel*> const message) const {
  CHECK(is_initialized());
  body_.WriteToMessage(message->mutable_body());
  prolongation_adaptive_step_parameters_.WriteToMessage(
      message->mutable_prolongation_adaptive_step_parameters());
  history_fixed_step_parameters_.WriteToMessage(
      message->mutable_history_fixed_step_parameters());
  history_->WriteToMessage(message->mutable_history(), {prolongation_});
  prediction_->Fork().time().WriteToMessage(
      message->mutable_prediction_fork_time());
  prediction_->last().time().WriteToMessage(
      message->mutable_prediction_last_time());
  prediction_adaptive_step_parameters_.WriteToMessage(
      message->mutable_prediction_adaptive_step_parameters());
  if (flight_plan_ != nullptr) {
    flight_plan_->WriteToMessage(message->mutable_flight_plan());
  }
  message->set_is_dirty(is_dirty_);
}

not_null<std::unique_ptr<Vessel>> Vessel::ReadFromMessage(
    serialization::Vessel const& message,
    not_null<Ephemeris<Barycentric>*> const ephemeris,
    not_null<Celestial const*> const parent) {
  // NOTE(egg): for now we do not read the |MasslessBody| as it can contain no
  // information.
  std::unique_ptr<Vessel> vessel;
  bool const is_pre_буняковский = message.has_history_and_prolongation() ||
                                  message.has_owned_prolongation();

  if (is_pre_буняковский) {
    vessel = make_not_null_unique<Vessel>(parent,
                                          ephemeris,
                                          DefaultHistoryParameters(),
                                          DefaultProlongationParameters(),
                                          DefaultPredictionParameters());

    if (message.has_history_and_prolongation()) {
      vessel->history_ =
          DiscreteTrajectory<Barycentric>::ReadFromMessage(
              message.history_and_prolongation().history(), /*forks=*/{});
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
    } else {
      vessel->history_ = DiscreteTrajectory<Barycentric>::ReadFromMessage(
                             message.owned_prolongation(), /*forks=*/{});
      vessel->prolongation_ = vessel->history_->NewForkAtLast();
      CHECK(!message.has_prediction());
      CHECK(!message.has_flight_plan());
    }
    if (vessel->prediction_ == nullptr) {
      vessel->prediction_ = vessel->history_->NewForkAtLast();
    }
  } else {
    CHECK(message.has_history() &&
          message.has_history_fixed_step_parameters() &&
          message.has_prolongation_adaptive_step_parameters() &&
          message.has_prediction_fork_time() &&
          message.has_prediction_last_time() &&
          message.has_prediction_adaptive_step_parameters())
        << message.DebugString();
    vessel = make_not_null_unique<Vessel>(
        parent,
        ephemeris,
        Ephemeris<Barycentric>::FixedStepParameters::ReadFromMessage(
            message.history_fixed_step_parameters()),
        Ephemeris<Barycentric>::AdaptiveStepParameters::ReadFromMessage(
            message.prolongation_adaptive_step_parameters()),
        Ephemeris<Barycentric>::AdaptiveStepParameters::ReadFromMessage(
            message.prediction_adaptive_step_parameters()));
    vessel->history_ = DiscreteTrajectory<Barycentric>::ReadFromMessage(
        message.history(), {&vessel->prolongation_});
    vessel->prediction_ = vessel->history_->NewForkWithoutCopy(
        Instant::ReadFromMessage(message.prediction_fork_time()));
    vessel->FlowPrediction(
        Instant::ReadFromMessage(message.prediction_last_time()));
    if (message.has_flight_plan()) {
      vessel->flight_plan_ = FlightPlan::ReadFromMessage(
          message.flight_plan(), vessel->history_.get(), ephemeris);
    }
    vessel->is_dirty_ = message.is_dirty();
  }
  return std::move(vessel);
}

Vessel::Vessel()
    : body_(),
      history_fixed_step_parameters_(DefaultHistoryParameters()),
      prolongation_adaptive_step_parameters_(DefaultProlongationParameters()),
      prediction_adaptive_step_parameters_(DefaultPredictionParameters()),
      parent_(testing_utilities::make_not_null<Celestial const*>()),
      ephemeris_(testing_utilities::make_not_null<Ephemeris<Barycentric>*>()),
      subset_node_(make_not_null_unique<Subset<Vessel>::Node>()) {}

void Vessel::AdvanceHistoryIfNeeded(Instant const& time) {
  Instant const& history_last_time = history_->last().time();
  Time const& Δt = history_fixed_step_parameters_.step();

  if (history_last_time + Δt < time) {
    if (is_dirty_) {
      FlowProlongation(history_last_time + Δt);
      history_->Append(history_last_time + Δt,
                       prolongation_->last().degrees_of_freedom());
      is_dirty_ = false;
    }
    FlowHistory(time);
    history_->DeleteFork(prolongation_);
    prolongation_ = history_->NewForkAtLast();
  }
}

void Vessel::FlowHistory(Instant const& time) {
  ephemeris_->FlowWithFixedStep(
      {history_.get()},
      Ephemeris<Barycentric>::NoIntrinsicAccelerations,
      time,
      history_fixed_step_parameters_);
}

void Vessel::FlowProlongation(Instant const& time) {
  Instant const& prolongation_last_time = prolongation_->last().time();
  CHECK_LE(prolongation_last_time, time);
  if (prolongation_last_time == time) {
    return;
  }
  ephemeris_->FlowWithAdaptiveStep(
      prolongation_,
      Ephemeris<Barycentric>::NoIntrinsicAcceleration,
      time,
      prolongation_adaptive_step_parameters_,
      Ephemeris<Barycentric>::unlimited_max_ephemeris_steps);
}

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
        FlightPlan::max_ephemeris_steps_per_frame);
    if (!finite_time && reached_t) {
      // This will prolong the ephemeris by |max_ephemeris_steps_per_frame|.
      ephemeris_->FlowWithAdaptiveStep(
        prediction_,
        Ephemeris<Barycentric>::NoIntrinsicAcceleration,
        time,
        prediction_adaptive_step_parameters_,
        FlightPlan::max_ephemeris_steps_per_frame);
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
