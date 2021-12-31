
#include "ksp_plugin/vessel.hpp"

#include <algorithm>
#include <functional>
#include <limits>
#include <list>
#include <string>
#include <vector>

#include "base/map_util.hpp"
#include "geometry/named_quantities.hpp"
#include "ksp_plugin/integrators.hpp"
#include "ksp_plugin/pile_up.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/make_not_null.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_vessel {

using base::Contains;
using base::FindOrDie;
using base::jthread;
using base::make_not_null_unique;
using base::MakeStoppableThread;
using geometry::BarycentreCalculator;
using geometry::InfiniteFuture;
using geometry::InfinitePast;
using geometry::Position;
using quantities::IsFinite;
using quantities::Length;
using quantities::Time;
using quantities::si::Metre;
using ::std::placeholders::_1;

using namespace std::chrono_literals;

bool operator!=(Vessel::PrognosticatorParameters const& left,
                Vessel::PrognosticatorParameters const& right) {
  return left.first_time != right.first_time ||
         left.first_degrees_of_freedom != right.first_degrees_of_freedom ||
         &left.adaptive_step_parameters.integrator() !=
             &right.adaptive_step_parameters.integrator() ||
         left.adaptive_step_parameters.max_steps() !=
             right.adaptive_step_parameters.max_steps() ||
         left.adaptive_step_parameters.length_integration_tolerance() !=
             right.adaptive_step_parameters.length_integration_tolerance() ||
         left.adaptive_step_parameters.speed_integration_tolerance() !=
             right.adaptive_step_parameters.speed_integration_tolerance();
}

Vessel::Vessel(GUID guid,
               std::string name,
               not_null<Celestial const*> const parent,
               not_null<Ephemeris<Barycentric>*> const ephemeris,
               Ephemeris<Barycentric>::AdaptiveStepParameters
                   prediction_adaptive_step_parameters)
    : guid_(std::move(guid)),
      name_(std::move(name)),
      body_(),
      prediction_adaptive_step_parameters_(
          std::move(prediction_adaptive_step_parameters)),
      parent_(parent),
      ephemeris_(ephemeris),
      history_(trajectory_.segments().begin()),
      psychohistory_(trajectory_.segments().end()),
      prediction_(trajectory_.segments().end()),
      prognosticator_(
          [this](PrognosticatorParameters const& parameters) {
            return FlowPrognostication(parameters);
          },
          20ms)  // 50 Hz.
{}

Vessel::~Vessel() {
  LOG(INFO) << "Destroying vessel " << ShortDebugString();
  // Ask the prognosticator to shut down.  This may take a while.
  StopPrognosticator();
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

bool Vessel::WillKeepPart(PartId const id) const {
  return Contains(kept_parts_, id);
}

void Vessel::FreeParts() {
  CHECK_LE(kept_parts_.size(), parts_.size());
  for (auto it = parts_.begin(); it != parts_.end();) {
    not_null<Part*> const part = it->second.get();
    if (Contains(kept_parts_, part->part_id())) {
      ++it;
    } else {
      part->reset_containing_pile_up();
      it = parts_.erase(it);
    }
  }
  CHECK(!parts_.empty());
  kept_parts_.clear();
}

void Vessel::ClearAllIntrinsicForcesAndTorques() {
  for (auto const& [_, part] : parts_) {
    part->clear_intrinsic_force();
    part->clear_intrinsic_torque();
  }
}

void Vessel::PrepareHistory(
    Instant const& t,
    DiscreteTrajectorySegment<Barycentric>::DownsamplingParameters const&
        downsampling_parameters) {
  CHECK(!parts_.empty());
  if (trajectory_.empty()) {
    LOG(INFO) << "Preparing history of vessel " << ShortDebugString()
              << " at " << t;
    BarycentreCalculator<DegreesOfFreedom<Barycentric>, Mass> calculator;
    ForAllParts([&calculator](Part& part) {
      calculator.Add(
          part.rigid_motion()({RigidPart::origin, RigidPart::unmoving}),
          part.mass());
    });
    CHECK(psychohistory_ == trajectory_.segments().end());
    history_->SetDownsampling(downsampling_parameters);
    trajectory_.Append(t, calculator.Get()).IgnoreError();
    psychohistory_ = trajectory_.NewSegment();
    prediction_ = trajectory_.NewSegment();
  }
}

void Vessel::DisableDownsampling() {
  history_->ClearDownsampling();
}

not_null<Part*> Vessel::part(PartId const id) const {
  return FindOrDie(parts_, id).get();
}

void Vessel::ForSomePart(std::function<void(Part&)> action) const {
  CHECK(!parts_.empty());
  action(*parts_.begin()->second);
}

void Vessel::ForAllParts(std::function<void(Part&)> action) const {
  for (auto const& [_, part] : parts_) {
    action(*part);
  }
}

DiscreteTrajectory<Barycentric> const& Vessel::trajectory() const {
  return trajectory_;
}

DiscreteTrajectorySegmentIterator<Barycentric> Vessel::history() const {
  return history_;
}

DiscreteTrajectorySegmentIterator<Barycentric> Vessel::psychohistory() const {
  return psychohistory_;
}

DiscreteTrajectorySegmentIterator<Barycentric> Vessel::prediction() const {
  return prediction_;
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

bool Vessel::has_flight_plan() const {
  return has_deserialized_flight_plan() ||
         std::holds_alternative<serialization::FlightPlan>(flight_plan_);
}

FlightPlan& Vessel::flight_plan() const {
  CHECK(has_deserialized_flight_plan());
  auto& flight_plan =
      *std::get<std::unique_ptr<FlightPlan>>(flight_plan_);
  return flight_plan;
}

void Vessel::ReadFlightPlanFromMessage() {
  if (std::holds_alternative<serialization::FlightPlan>(flight_plan_)) {
    auto const& message = std::get<serialization::FlightPlan>(flight_plan_);
    flight_plan_ = FlightPlan::ReadFromMessage(message, ephemeris_);
  }
}

void Vessel::AdvanceTime() {
  // Squirrel away the prediction so that we can reattach it if we don't have a
  // prognostication.
  auto prediction = trajectory_.DetachSegments(prediction_);
  prediction_ = trajectory_.segments().end();

  // Read the wall of text below and realize that this can happen for the
  // history as well as the psychohistory, if the history of the part was
  // obtained using an adaptive step integrator, which is the case during a
  // burn.  See #2931.
  trajectory_.DeleteSegments(psychohistory_);
  AppendToVesselTrajectory(&Part::history_begin,
                           &Part::history_end,
                           *history_);
  psychohistory_ = trajectory_.NewSegment();

  // The reason why we may want to skip the start of the psychohistory is
  // subtle.  Say that we have a vessel A with points at t₀, t₀ + 10 s,
  // t₀ + 20 s in its history.  Say that a vessel B is created at t₀ + 23 s,
  // maybe because of an undocking or a staging, and (some) of the parts of A
  // are transfered to B.  Now time moves a bit and at t₀ = 24 s we want to
  // attach the psychohistory of these parts (their centre of mass, really) to
  // vessel B.  Most of the time the psychohistory will have a *single* point at
  // t₀ + 24 s and everything will be fine.  However, because the psychohistory
  // is integrated using an adaptive step, it is possible that it would have
  // multiple points, say one at t₀ + 21 s and one at t₀ + 24 s.  In this case
  // trying to insert the point at t₀ + 21 s would put us before the last point
  // of the history of B and would fail a check.  Therefore, we just ignore that
  // point.  See #2507 and the |last_time| in AppendToVesselTrajectory.
  AppendToVesselTrajectory(&Part::psychohistory_begin,
                           &Part::psychohistory_end,
                           *psychohistory_);

  // Attach the prognostication, if there is one.  Otherwise fall back to the
  // pre-existing prediction.
  auto optional_prognostication = prognosticator_.Get();
  if (optional_prognostication.has_value()) {
    AttachPrediction(std::move(optional_prognostication.value()));
  } else {
    AttachPrediction(std::move(prediction));
  }

  for (auto const& [_, part] : parts_) {
    part->ClearHistory();
  }
}

void Vessel::CreateFlightPlan(
    Instant const& final_time,
    Mass const& initial_mass,
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        flight_plan_adaptive_step_parameters,
    Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters const&
        flight_plan_generalized_adaptive_step_parameters) {
  auto const history_back = history_->back();
  flight_plan_ = std::make_unique<FlightPlan>(
      initial_mass,
      /*initial_time=*/history_back.time,
      /*initial_degrees_of_freedom=*/history_back.degrees_of_freedom,
      final_time,
      ephemeris_,
      flight_plan_adaptive_step_parameters,
      flight_plan_generalized_adaptive_step_parameters);
}

void Vessel::DeleteFlightPlan() {
  flight_plan_.emplace<std::unique_ptr<FlightPlan>>();
}

absl::Status Vessel::RebaseFlightPlan(Mass const& initial_mass) {
  CHECK(has_deserialized_flight_plan());
  auto& flight_plan =
      std::get<std::unique_ptr<FlightPlan>>(flight_plan_);
  Instant const new_initial_time = history_->back().time;
  int first_manœuvre_kept = 0;
  for (int i = 0; i < flight_plan->number_of_manœuvres(); ++i) {
    auto const& manœuvre = flight_plan->GetManœuvre(i);
    if (manœuvre.initial_time() < new_initial_time) {
      first_manœuvre_kept = i + 1;
      if (new_initial_time < manœuvre.final_time()) {
        return absl::UnavailableError(
            u8"Cannot rebase during planned manœuvre execution");
      }
    }
  }
  not_null<std::unique_ptr<FlightPlan>> const original_flight_plan =
      std::move(flight_plan);
  Instant const new_desired_final_time =
      new_initial_time >= original_flight_plan->desired_final_time()
          ? new_initial_time + (original_flight_plan->desired_final_time() -
                                original_flight_plan->initial_time())
          : original_flight_plan->desired_final_time();
  CreateFlightPlan(
      new_desired_final_time,
      initial_mass,
      original_flight_plan->adaptive_step_parameters(),
      original_flight_plan->generalized_adaptive_step_parameters());
  for (int i = first_manœuvre_kept;
       i < original_flight_plan->number_of_manœuvres();
       ++i) {
    auto const& manœuvre = original_flight_plan->GetManœuvre(i);
    flight_plan->Insert(manœuvre.burn(), i - first_manœuvre_kept)
        .IgnoreError();
  }
  return absl::OkStatus();
}

void Vessel::RefreshPrediction() {
  // The |prognostication| is a trajectory which is computed asynchronously and
  // may be used as a prediction;
  std::optional<DiscreteTrajectory<Barycentric>> prognostication;

  // Note that we know that |RefreshPrediction| is called on the main thread,
  // therefore the ephemeris currently covers the last time of the
  // psychohistory.  Were this to change, this code might have to change.
  PrognosticatorParameters prognosticator_parameters{
      psychohistory_->back().time,
      psychohistory_->back().degrees_of_freedom,
      prediction_adaptive_step_parameters_};
  if (synchronous_) {
    auto status_or_prognostication =
        FlowPrognostication(std::move(prognosticator_parameters));
    if (status_or_prognostication.ok()) {
      prognostication = std::move(status_or_prognostication).value();
    }
  } else {
    prognosticator_.Put(std::move(prognosticator_parameters));
    prognosticator_.Start();
    prognostication = prognosticator_.Get();
  }
  if (prognostication.has_value()) {
    AttachPrediction(std::move(prognostication).value());
  }
}

void Vessel::RefreshPrediction(Instant const& time) {
  RefreshPrediction();
  trajectory_.ForgetAfter(trajectory_.upper_bound(time));
}

void Vessel::StopPrognosticator() {
  prognosticator_.Stop();
}

std::string Vessel::ShortDebugString() const {
  return name_ + " (" + guid_ + ")";
}

void Vessel::WriteToMessage(not_null<serialization::Vessel*> const message,
                            PileUp::SerializationIndexForPileUp const&
                                serialization_index_for_pile_up) const {
  message->set_guid(guid_);
  message->set_name(name_);
  body_.WriteToMessage(message->mutable_body());
  prediction_adaptive_step_parameters_.WriteToMessage(
      message->mutable_prediction_adaptive_step_parameters());
  for (auto const& [_, part] : parts_) {
    part->WriteToMessage(message->add_parts(), serialization_index_for_pile_up);
  }
  for (auto const& part_id : kept_parts_) {
    CHECK(Contains(parts_, part_id));
    message->add_kept_parts(part_id);
  }
  // Starting with Gateaux we don't save the prediction, see #2685.  Instead we
  // just save its first point and re-read as if it was the whole prediction.
  trajectory_.WriteToMessage(
      message->mutable_history(),
      /*begin=*/trajectory_.begin(),
      /*end=*/std::next(prediction_->begin()),
      /*tracked=*/{history_, psychohistory_, prediction_},
      /*exact=*/{});
  if (std::holds_alternative<serialization::FlightPlan>(flight_plan_)) {
    *message->mutable_flight_plan() =
        std::get<serialization::FlightPlan>(flight_plan_);
  } else if (std::holds_alternative<std::unique_ptr<FlightPlan>>(
                 flight_plan_)) {
    auto& flight_plan = std::get<std::unique_ptr<FlightPlan>>(flight_plan_);
    if (flight_plan != nullptr) {
      flight_plan->WriteToMessage(message->mutable_flight_plan());
    }
  } else {
    LOG(FATAL) << "Unexpected flight plan variant " << flight_plan_.index();
  }
}

not_null<std::unique_ptr<Vessel>> Vessel::ReadFromMessage(
    serialization::Vessel const& message,
    not_null<Celestial const*> const parent,
    not_null<Ephemeris<Barycentric>*> const ephemeris,
    std::function<void(PartId)> const& deletion_callback) {
  bool const is_pre_cesàro = message.has_psychohistory_is_authoritative();
  bool const is_pre_chasles = message.has_prediction();
  bool const is_pre_陈景润 = !message.history().has_downsampling() &&
                             message.history().segment_size() == 0;
  bool const is_pre_hamilton = message.history().segment_size() == 0;
  LOG_IF(WARNING, is_pre_hamilton)
      << "Reading pre-"
      << (is_pre_cesàro    ? u8"Cesàro"
          : is_pre_chasles ? "Chasles"
          : is_pre_陈景润   ? u8"陈景润"
                           : "Hamilton") << " Vessel";

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

  if (is_pre_cesàro) {
    auto const psychohistory =
        DiscreteTrajectory<Barycentric>::ReadFromMessage(message.history(),
                                                         /*forks=*/{});
    // The |history_| has been created by the constructor above.  Reconstruct
    // it from the |psychohistory|.
    for (auto it = psychohistory.begin(); it != psychohistory.end();) {
      auto const& [time, degrees_of_freedom] = *it;
      ++it;
      if (it == psychohistory.end() &&
          !message.psychohistory_is_authoritative()) {
        vessel->psychohistory_ = vessel->trajectory_.NewSegment();
      }
      vessel->trajectory_.Append(time, degrees_of_freedom).IgnoreError();
    }
    if (message.psychohistory_is_authoritative()) {
      vessel->psychohistory_ = vessel->trajectory_.NewSegment();
    }
    vessel->prediction_ = vessel->trajectory_.NewSegment();
  } else if (is_pre_chasles) {
    vessel->trajectory_ = DiscreteTrajectory<Barycentric>::ReadFromMessage(
        message.history(),
        /*tracked=*/{&vessel->psychohistory_});
    vessel->history_ = vessel->trajectory_.segments().begin();
    vessel->prediction_ = vessel->trajectory_.NewSegment();
  } else if (is_pre_hamilton) {
    vessel->trajectory_ = DiscreteTrajectory<Barycentric>::ReadFromMessage(
        message.history(),
        /*tracked=*/{&vessel->psychohistory_, &vessel->prediction_});
    vessel->history_ = vessel->trajectory_.segments().begin();
    // Necessary after Εὔδοξος because the ephemeris has not been prolonged
    // during deserialization.  Doesn't hurt prior to Εὔδοξος.
    ephemeris->Prolong(vessel->prediction_->back().time).IgnoreError();
  } else {
    vessel->trajectory_ = DiscreteTrajectory<Barycentric>::ReadFromMessage(
        message.history(),
        /*tracked=*/{&vessel->history_,
                     &vessel->psychohistory_,
                     &vessel->prediction_});
    // Necessary after Εὔδοξος because the ephemeris has not been prolonged
    // during deserialization.
    ephemeris->Prolong(vessel->prediction_->back().time).IgnoreError();
  }

  if (is_pre_陈景润) {
    vessel->history_->SetDownsamplingUnconditionally(
        DefaultDownsamplingParameters());
  }

  if (message.has_flight_plan()) {
    // Starting with हरीश चंद्र we deserialize the flight plan lazily.
    vessel->flight_plan_
        .emplace<serialization::FlightPlan>(message.flight_plan());
  }
  return vessel;
}

void Vessel::FillContainingPileUpsFromMessage(
    serialization::Vessel const& message,
    PileUp::PileUpForSerializationIndex const&
        pile_up_for_serialization_index) {
  for (auto const& part_message : message.parts()) {
    auto const& part = FindOrDie(parts_, part_message.part_id());
    part->FillContainingPileUpFromMessage(part_message,
                                          pile_up_for_serialization_index);
  }
}

void Vessel::RequestOrbitAnalysis(Time const& mission_duration) {
  if (!orbit_analyser_.has_value()) {
    // TODO(egg): perhaps we should get the history parameters from the plugin;
    // on the other hand, these are probably overkill for high orbits anyway,
    // and given that we know many things about our trajectory in the analyser,
    // perhaps we should pick something appropriate automatically instead.  The
    // default will do in the meantime.
    orbit_analyser_.emplace(ephemeris_, DefaultHistoryParameters());
  }
  if (orbit_analyser_->last_parameters().has_value() &&
      orbit_analyser_->last_parameters()->mission_duration !=
          mission_duration) {
    orbit_analyser_->Interrupt();
  }
  orbit_analyser_->RequestAnalysis(
      {.first_time = psychohistory_->back().time,
       .first_degrees_of_freedom = psychohistory_->back().degrees_of_freedom,
       .mission_duration = mission_duration});
}

void Vessel::ClearOrbitAnalyser() {
  orbit_analyser_.reset();
}

double Vessel::progress_of_orbit_analysis() const {
  if (!orbit_analyser_.has_value()) {
    return 0;
  }
  return orbit_analyser_->progress_of_next_analysis();
}

void Vessel::RefreshOrbitAnalysis() {
  if (orbit_analyser_.has_value()) {
    orbit_analyser_->RefreshAnalysis();
  }
}

OrbitAnalyser::Analysis* Vessel::orbit_analysis() {
  return orbit_analyser_.has_value() ? orbit_analyser_->analysis() : nullptr;
}

void Vessel::MakeAsynchronous() {
  synchronous_ = false;
}

void Vessel::MakeSynchronous() {
  synchronous_ = true;
}

Vessel::Vessel()
    : body_(),
      prediction_adaptive_step_parameters_(DefaultPredictionParameters()),
      parent_(testing_utilities::make_not_null<Celestial const*>()),
      ephemeris_(testing_utilities::make_not_null<Ephemeris<Barycentric>*>()),
      history_(trajectory_.segments().begin()),
      psychohistory_(trajectory_.segments().end()),
      prediction_(trajectory_.segments().end()),
      prognosticator_(nullptr, 20ms) {}

absl::StatusOr<DiscreteTrajectory<Barycentric>>
Vessel::FlowPrognostication(
    PrognosticatorParameters prognosticator_parameters) {
  DiscreteTrajectory<Barycentric> prognostication;
  prognostication.Append(
      prognosticator_parameters.first_time,
      prognosticator_parameters.first_degrees_of_freedom).IgnoreError();
  absl::Status status;
  status = ephemeris_->FlowWithAdaptiveStep(
      &prognostication,
      Ephemeris<Barycentric>::NoIntrinsicAcceleration,
      ephemeris_->t_max(),
      prognosticator_parameters.adaptive_step_parameters,
      FlightPlan::max_ephemeris_steps_per_frame);
  bool const reached_t_max = status.ok();
  if (reached_t_max) {
    // This will prolong the ephemeris by |max_ephemeris_steps_per_frame|.
    status = ephemeris_->FlowWithAdaptiveStep(
        &prognostication,
        Ephemeris<Barycentric>::NoIntrinsicAcceleration,
        InfiniteFuture,
        prognosticator_parameters.adaptive_step_parameters,
        FlightPlan::max_ephemeris_steps_per_frame);
  }
  LOG_IF_EVERY_N(INFO, !status.ok(), 50)
      << "Prognostication from " << prognosticator_parameters.first_time
      << " finished at " << prognostication.back().time << " with "
      << status.ToString() << " for " << ShortDebugString();
  if (absl::IsCancelled(status)) {
    return status;
  } else {
    // Unless we were stopped, ignore the status, which indicates a failure to
    // reach |t_max|, and provide a short prognostication.
    return std::move(prognostication);
  }
}

void Vessel::AppendToVesselTrajectory(
    TrajectoryIterator const part_trajectory_begin,
    TrajectoryIterator const part_trajectory_end,
    DiscreteTrajectorySegment<Barycentric> const& segment) {
  CHECK(!parts_.empty());
  std::vector<DiscreteTrajectory<Barycentric>::iterator> its;
  std::vector<DiscreteTrajectory<Barycentric>::iterator> ends;
  its.reserve(parts_.size());
  ends.reserve(parts_.size());
  for (auto const& [_, part] : parts_) {
    its.push_back((*part.*part_trajectory_begin)());
    ends.push_back((*part.*part_trajectory_end)());
  }

  // We cannot append a point before this time, see the comments in AdvanceTime.
  Instant const last_time =
      segment.empty() ? InfinitePast : segment.back().time;

  // Loop over the times of the trajectory.
  for (;;) {
    auto const& it0 = its[0];
    bool const at_end_of_part_trajectory = it0 == ends[0];
    Instant const first_time = at_end_of_part_trajectory ? Instant()
                                                         : it0->time;
    bool const can_be_appended = !at_end_of_part_trajectory &&
                                 first_time > last_time;

    // Loop over the parts at a given time.
    BarycentreCalculator<DegreesOfFreedom<Barycentric>, Mass> calculator;
    int i = 0;
    for (auto const& [_, part] : parts_) {
      auto& it = its[i];
      CHECK_EQ(at_end_of_part_trajectory, it == ends[i]);
      if (!at_end_of_part_trajectory) {
        calculator.Add(it->degrees_of_freedom, part->mass());
        CHECK_EQ(first_time, it->time);
        ++it;
      }
      ++i;
    }

    if (at_end_of_part_trajectory) {
      return;
    }

    // Append the parts' barycentre to the trajectory.
    if (can_be_appended) {
      DegreesOfFreedom<Barycentric> const vessel_degrees_of_freedom =
          calculator.Get();
      trajectory_.Append(first_time, vessel_degrees_of_freedom).IgnoreError();
    }
  }
}

void Vessel::AttachPrediction(DiscreteTrajectory<Barycentric>&& trajectory) {
  trajectory.ForgetBefore(psychohistory_->back().time);
  if (trajectory.empty()) {
    prediction_ = trajectory_.NewSegment();
  } else {
    if (prediction_ != trajectory_.segments().end()) {
      trajectory_.DeleteSegments(prediction_);
    }
    prediction_ = trajectory_.AttachSegments(std::move(trajectory));
  }
}

bool Vessel::has_deserialized_flight_plan() const {
  return std::holds_alternative<std::unique_ptr<FlightPlan>>(flight_plan_) &&
         std::get<std::unique_ptr<FlightPlan>>(flight_plan_) != nullptr;
}

// Run the prognostication in both synchronous and asynchronous mode in tests to
// avoid code rot.
#if defined(_DEBUG)
std::atomic_bool Vessel::synchronous_(true);
#else
std::atomic_bool Vessel::synchronous_(false);
#endif

}  // namespace internal_vessel
}  // namespace ksp_plugin
}  // namespace principia
