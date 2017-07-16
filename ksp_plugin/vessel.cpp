
#include "ksp_plugin/vessel.hpp"

#include <algorithm>
#include <limits>
#include <list>
#include <string>
#include <vector>

#include "ksp_plugin/integrators.hpp"
#include "ksp_plugin/pile_up.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/make_not_null.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_vessel {

using base::Contains;
using base::FindOrDie;
using base::make_not_null_unique;
using geometry::BarycentreCalculator;
using geometry::Position;
using quantities::IsFinite;
using quantities::Time;

Vessel::Vessel(GUID const& guid,
               std::string const& name,
               not_null<Celestial const*> const parent,
               not_null<Ephemeris<Barycentric>*> const ephemeris,
               Ephemeris<Barycentric>::AdaptiveStepParameters const&
                   prediction_adaptive_step_parameters)
    : guid_(guid),
      name_(name),
      body_(),
      prediction_adaptive_step_parameters_(prediction_adaptive_step_parameters),
      parent_(parent),
      ephemeris_(ephemeris),
      psychohistory_(make_not_null_unique<DiscreteTrajectory<Barycentric>>()),
      prediction_(make_not_null_unique<DiscreteTrajectory<Barycentric>>()) {}

Vessel::~Vessel() {
  LOG(INFO) << "Destroying vessel " << ShortDebugString();
  // The parts must remove themselves from their pile-ups *before* any of them
  // starts to destroy, otherwise |clear_pile_up| might access destroyed parts.
  for (auto const& pair : parts_) {
    auto const& part = pair.second;
    part->clear_pile_up();
  }
}

GUID const& Vessel::guid() const {
  return guid_;
}

std::string const& Vessel::name() const {
  return name_;
}

void Vessel::set_name(std::string const& new_name) {
  LOG(INFO) << "Vessel " << ShortDebugString() << " is now known as "
            << new_name;
  name_ = new_name;
}

not_null<MasslessBody const*> Vessel::body() const {
  return &body_;
}

not_null<Celestial const*> Vessel::parent() const {
  return parent_;
}

void Vessel::set_parent(not_null<Celestial const*> const parent) {
  LOG(INFO) << "Vessel " << ShortDebugString() << " switches parent from "
            << parent_->body()->name() << " to " << parent->body()->name();
  parent_ = parent;
}

void Vessel::AddPart(not_null<std::unique_ptr<Part>> part) {
  LOG(INFO) << "Adding part " << part->ShortDebugString() << " to vessel "
            << ShortDebugString();
  parts_.emplace(part->part_id(), std::move(part));
}

not_null<std::unique_ptr<Part>> Vessel::ExtractPart(PartId const id) {
  CHECK_LE(kept_parts_.size(), parts_.size());
  auto const it = parts_.find(id);
  CHECK(it != parts_.end()) << id;
  auto result = std::move(it->second);
  LOG(INFO) << "Extracting part " << result->ShortDebugString()
            << " from vessel " << ShortDebugString();
  parts_.erase(it);
  kept_parts_.erase(id);
  return result;
}

void Vessel::KeepPart(PartId const id) {
  CHECK_LE(kept_parts_.size(), parts_.size());
  CHECK(Contains(parts_, id)) << id;
  kept_parts_.insert(id);
}

void Vessel::FreeParts() {
  CHECK_LE(kept_parts_.size(), parts_.size());
  for (auto it = parts_.begin(); it != parts_.end();) {
    not_null<Part*> const part = it->second.get();
    if (Contains(kept_parts_, part->part_id())) {
      ++it;
    } else {
      part->clear_pile_up();
      it = parts_.erase(it);
    }
  }
  CHECK(!parts_.empty());
  kept_parts_.clear();
}

void Vessel::ClearAllIntrinsicForces() {
  for (auto const& pair : parts_) {
    auto const& part = pair.second;
    part->clear_intrinsic_force();
  }
}

void Vessel::PreparePsychohistory(Instant const& t) {
  CHECK(!parts_.empty());
  if (psychohistory_->Empty()) {
    LOG(INFO) << "Preparing psychohistory of vessel " << ShortDebugString()
              << " at " << t;
    BarycentreCalculator<DegreesOfFreedom<Barycentric>, Mass> calculator;
    ForAllParts([&calculator](Part& part) {
      calculator.Add(part.degrees_of_freedom(), part.mass());
    });
    psychohistory_->Append(t, calculator.Get());
    psychohistory_is_authoritative_ = true;
  }
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

void Vessel::AdvanceTime() {
  CHECK(!parts_.empty());
  std::vector<DiscreteTrajectory<Barycentric>::Iterator> its;
  std::vector<DiscreteTrajectory<Barycentric>::Iterator> tails;
  its.reserve(parts_.size());
  tails.reserve(parts_.size());
  for (auto const& pair : parts_) {
    Part const& part = *pair.second;
    CHECK(!part.tail().Empty()) << part.ShortDebugString()
                                << " " << ShortDebugString();
    its.push_back(part.tail().Begin());
    tails.push_back(part.tail().last());
  }

  Part const& first_part = *parts_.begin()->second;
  bool const tail_is_authoritative = first_part.tail_is_authoritative();

  // Loop over the times of the tails.
  for (;;) {
    Instant const first_time = its[0].time();
    bool const at_end_of_tail = its[0] == tails[0];

    // Loop over the parts at a given time.
    BarycentreCalculator<DegreesOfFreedom<Barycentric>, Mass> calculator;
    int i = 0;
    for (auto const& pair : parts_) {
      Part& part = *pair.second;
      auto& it = its[i];
      calculator.Add(it.degrees_of_freedom(), part.mass());
      CHECK_EQ(first_time, it.time());
      CHECK_EQ(at_end_of_tail, it == tails[i]);
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
        first_time,
        vessel_degrees_of_freedom,
        /*authoritative=*/!at_end_of_tail || tail_is_authoritative);

    if (at_end_of_tail) {
      return;
    }
  }
}

void Vessel::ForgetBefore(Instant const& time) {
  // Make sure that the psychohistory keep at least an authoritative point (and
  // possibly a non-authoritative one).  We cannot use the parts because they
  // may have been moved to the future already.
  psychohistory_->ForgetBefore(std::min(time, last_authoritative().time()));
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
  CHECK(!psychohistory_->Empty());
  auto const last = psychohistory_->last();
  prediction_->Append(last.time(), last.degrees_of_freedom());
  FlowPrediction(last_time);
}

DiscreteTrajectory<Barycentric>::Iterator Vessel::last_authoritative() const {
  auto it = psychohistory_->last();
  if (!psychohistory_is_authoritative_) {
    --it;
  }
  return it;
}

DiscreteTrajectory<Barycentric> const& Vessel::psychohistory() const {
  return *psychohistory_;
}

bool Vessel::psychohistory_is_authoritative() const {
  return psychohistory_is_authoritative_;
}

void Vessel::WriteToMessage(
    not_null<serialization::Vessel*> const message) const {
  message->set_guid(guid_);
  message->set_name(name_);
  body_.WriteToMessage(message->mutable_body());
  prediction_adaptive_step_parameters_.WriteToMessage(
      message->mutable_prediction_adaptive_step_parameters());
  for (auto const& pair : parts_) {
    auto const& part = pair.second;
    part->WriteToMessage(message->add_parts());
  }
  for (auto const& part_id : kept_parts_) {
    CHECK(Contains(parts_, part_id));
    message->add_kept_parts(part_id);
  }
  psychohistory_->WriteToMessage(message->mutable_psychohistory(),
                                 /*forks=*/{});
  message->set_psychohistory_is_authoritative(psychohistory_is_authoritative_);
  prediction_->WriteToMessage(message->mutable_prediction(),
                              /*forks=*/{});
  if (flight_plan_ != nullptr) {
    flight_plan_->WriteToMessage(message->mutable_flight_plan());
  }
}

not_null<std::unique_ptr<Vessel>> Vessel::ReadFromMessage(
    serialization::Vessel const& message,
    not_null<Celestial const*> const parent,
    not_null<Ephemeris<Barycentric>*> const ephemeris,
    std::function<void(PartId)> const& deletion_callback) {
  // NOTE(egg): for now we do not read the |MasslessBody| as it can contain no
  // information.
  auto vessel = make_not_null_unique<Vessel>(
      message.guid(),
      message.name(),
      parent,
      ephemeris,
      Ephemeris<Barycentric>::AdaptiveStepParameters::ReadFromMessage(
          message.prediction_adaptive_step_parameters()));
  for (auto const& serialized_part : message.parts()) {
    PartId const part_id = serialized_part.part_id();
    auto part =
        Part::ReadFromMessage(serialized_part, [deletion_callback, part_id]() {
          if (deletion_callback != nullptr) {
            deletion_callback(part_id);
          }
        });
    vessel->parts_.emplace(part_id, std::move(part));
  }
  for (PartId const part_id : message.kept_parts()) {
    CHECK(Contains(vessel->parts_, part_id));
    vessel->kept_parts_.insert(part_id);
  }
  vessel->psychohistory_ =
      DiscreteTrajectory<Barycentric>::ReadFromMessage(message.psychohistory(),
                                                       /*forks=*/{});
  vessel->psychohistory_is_authoritative_ =
      message.psychohistory_is_authoritative();
  vessel->prediction_ =
      DiscreteTrajectory<Barycentric>::ReadFromMessage(message.prediction(),
                                                       /*forks=*/{});
  if (message.has_flight_plan()) {
    vessel->flight_plan_ = FlightPlan::ReadFromMessage(message.flight_plan(),
                                                       ephemeris);
  }
  return std::move(vessel);
}

void Vessel::FillContainingPileUpsFromMessage(
    serialization::Vessel const& message,
    not_null<std::list<PileUp>*> const pile_ups) {
  for (auto const& part_message : message.parts()) {
    auto const& part = FindOrDie(parts_, part_message.part_id());
    part->FillContainingPileUpFromMessage(part_message, pile_ups);
  }
}

std::string Vessel::ShortDebugString() const {
  return name_ + " (" + guid_ + ")";
}

Vessel::Vessel()
    : body_(),
      prediction_adaptive_step_parameters_(DefaultPredictionParameters()),
      parent_(testing_utilities::make_not_null<Celestial const*>()),
      ephemeris_(testing_utilities::make_not_null<Ephemeris<Barycentric>*>()),
      psychohistory_(make_not_null_unique<DiscreteTrajectory<Barycentric>>()),
      prediction_(make_not_null_unique<DiscreteTrajectory<Barycentric>>()) {}

void Vessel::AppendToPsychohistory(
    Instant const& time,
    DegreesOfFreedom<Barycentric> const& degrees_of_freedom,
    bool const authoritative) {
  if (!psychohistory_is_authoritative_) {
    CHECK(!psychohistory_->Empty());
    psychohistory_->ForgetAfter((--psychohistory_->last()).time());
  }
  psychohistory_->Append(time, degrees_of_freedom);
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

}  // namespace internal_vessel
}  // namespace ksp_plugin
}  // namespace principia
