#include "ksp_plugin/vessel.hpp"

#include <algorithm>
#include <functional>
#include <limits>
#include <list>
#include <set>
#include <string>
#include <vector>

#include "absl/container/btree_set.h"
#include "base/macros.hpp"
#include "base/map_util.hpp"
#include "geometry/named_quantities.hpp"
#include "ksp_plugin/integrators.hpp"
#include "ksp_plugin/pile_up.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/make_not_null.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_vessel {

using physics::Client;
using quantities::IsFinite;
using quantities::Length;
using quantities::Time;
using quantities::Torque;
using quantities::si::Metre;
using ::std::placeholders::_1;
using namespace principia::base::_jthread;
using namespace principia::base::_map_util;
using namespace principia::base::_not_null;
using namespace principia::geometry::_barycentre_calculator;
using namespace principia::geometry::_named_quantities;

using namespace std::chrono_literals;

// TODO(phl): Move this to some kind of parameters.
constexpr std::int64_t max_points_to_serialize = 20'000;

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

Vessel::Vessel(
    GUID guid,
    std::string name,
    not_null<Celestial const*> const parent,
    not_null<Ephemeris<Barycentric>*> const ephemeris,
    Ephemeris<Barycentric>::AdaptiveStepParameters
        prediction_adaptive_step_parameters,
    DiscreteTrajectorySegment<Barycentric>::DownsamplingParameters const&
        downsampling_parameters)
    : guid_(std::move(guid)),
      name_(std::move(name)),
      body_(),
      prediction_adaptive_step_parameters_(
          std::move(prediction_adaptive_step_parameters)),
      parent_(parent),
      ephemeris_(ephemeris),
      downsampling_parameters_(downsampling_parameters),
      checkpointer_(make_not_null_unique<Checkpointer<serialization::Vessel>>(
          MakeCheckpointerWriter(),
          MakeCheckpointerReader())),
      reanimator_(
          [this](Instant const& desired_t_min) {
            return Reanimate(desired_t_min);
          },
          20ms),  // 50 Hz.
      reanimator_clientele_(/*default_value=*/InfiniteFuture),
      backstory_(trajectory_.segments().begin()),
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
  reanimator_.Stop();
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

void Vessel::DetectCollapsibilityChange() {
  bool const will_be_collapsible = IsCollapsible();

  // It is always correct to mark as non-collapsible a collapsible segment or to
  // append collapsible points to a non-collapsible segment (but not
  // vice-versa).  If a non-collapsible segment is being closed but is very
  // short, we don't actually close it but keep appending points to it until it
  // is long enough to have been downsampled.  If downsampling is disabled,
  // surely this is not going to happen so no point in waiting for Godot.
  bool const collapsibility_changes = is_collapsible_ != will_be_collapsible;
  bool const becomes_non_collapsible = collapsibility_changes &&
                                       !will_be_collapsible;
  bool const awaits_first_downsampling = downsampling_parameters_.has_value() &&
                                         !backstory_->was_downsampled();

  if (collapsibility_changes &&
      (becomes_non_collapsible || !awaits_first_downsampling)) {
    // If collapsibility changes, we create a new history segment.  This ensures
    // that downsampling does not change collapsibility boundaries.

    // In normal situations we create a new segment with the collapsibility
    // given by |will_be_collapsible|.  In one cornercase we delete the current
    // segment.
    enum {
      Create,
      Delete,
    } segment_action = Create;

    if (!is_collapsible_) {
      // If the segment that is being closed is not collapsible, we have no way
      // to reconstruct it, so we must serialize it in a checkpoint.  Note that
      // the last point of the backstory specifies the initial conditions of the
      // next (collapsible) segment.
      Instant const checkpoint = backstory_->back().time;

      // In some cornercases we might try to create multiple checkpoints at the
      // same time, see #3280.  The checkpointer doesn't support that.
      bool const create_checkpoint =
          checkpointer_->newest_checkpoint() < checkpoint;

      if (create_checkpoint) {
        LOG(INFO) << "Writing " << ShortDebugString()
                  << " to checkpoint at: " << checkpoint;
        checkpointer_->WriteToCheckpoint(checkpoint);

        // If there are no checkpoints in the current trajectory (this would
        // happen if we restored the last part of trajectory and it didn't
        // overlap with a checkpoint and no reanimation happened) then the
        // |oldest_reanimated_checkpoint_| need to be updated to reflect the
        // newly created checkpoint.
        absl::MutexLock l(&lock_);
        if (oldest_reanimated_checkpoint_ == InfiniteFuture) {
          oldest_reanimated_checkpoint_ = checkpoint;
        } else {
          CHECK_LT(oldest_reanimated_checkpoint_, checkpoint);
        }
      } else {
        // Not only don't we create a new checkpoint and a new segment, but we
        // also delete the current, non-collapsible, 1-point segment, so that we
        // keep appending to the previous collapsible, 1-point segment.  See
        // #3332.
        LOG(INFO) << "Not writing " << ShortDebugString()
                  << " to duplicate checkpoint at: " << checkpoint;
        segment_action = Delete;
      }
    }

    auto psychohistory = trajectory_.DetachSegments(psychohistory_);
    switch (segment_action) {
      case Create: {
        backstory_ = trajectory_.NewSegment();
        if (downsampling_parameters_.has_value()) {
          backstory_->SetDownsampling(downsampling_parameters_.value());
        }
        break;
      }
      case Delete: {
        // Let's hope that no-one has kept an iterator to the deleted backstory.
        trajectory_.DeleteSegments(backstory_);
        CHECK(!trajectory_.segments().empty());
        backstory_ = std::prev(trajectory_.segments().end());
        break;
      }
    };
    psychohistory_ = trajectory_.AttachSegments(std::move(psychohistory));

    // Not updated if we chose to append to the current segment.
    is_collapsible_ = will_be_collapsible;
  }
}

void Vessel::CreateTrajectoryIfNeeded(Instant const& t) {
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
    if (downsampling_parameters_.has_value()) {
      backstory_->SetDownsampling(downsampling_parameters_.value());
    }
    trajectory_.Append(t, calculator.Get()).IgnoreError();
    psychohistory_ = trajectory_.NewSegment();
    prediction_ = trajectory_.NewSegment();
  }
}

void Vessel::DisableDownsampling() {
  backstory_->ClearDownsampling();
  // From now on, no downsampling will happen.
  downsampling_parameters_ = std::nullopt;
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
  return !flight_plans_.empty();
}

int Vessel::flight_plan_count() const {
  return flight_plans_.size();
}

int Vessel::selected_flight_plan_index() const {
  return selected_flight_plan_index_;
}

void Vessel::SelectFlightPlan(int index) {
  CHECK_GE(index, 0);
  CHECK_LT(index, flight_plan_count());
  selected_flight_plan_index_ = index;
}

FlightPlan& Vessel::flight_plan() const {
  CHECK(has_deserialized_flight_plan());
  auto& flight_plan =
      *std::get<not_null<std::unique_ptr<FlightPlan>>>(selected_flight_plan());
  return flight_plan;
}

void Vessel::ReadFlightPlanFromMessage() {
  if (!flight_plans_.empty() &&
      std::holds_alternative<serialization::FlightPlan>(
          selected_flight_plan())) {
    auto const& message =
        std::get<serialization::FlightPlan>(selected_flight_plan());
    selected_flight_plan() = FlightPlan::ReadFromMessage(message, ephemeris_);
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
                           *backstory_);
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

void Vessel::RequestReanimation(Instant const& desired_t_min) {
  reanimator_.Start();

  // No locking here because vessel reanimation is only invoked from the main
  // thread.
  bool must_restart;
  {
    // If the reanimator is asked to do significantly less work (in terms of
    // checkpoints to reanimate) than it is currently doing, interrupt it.  Note
    // that this is fundamentally racy: for instance the reanimator may not have
    // picked the last input given by Put.  But it helps if the user was doing a
    // very long reanimation and wants to shorten it.  Note however that we must
    // not move the desired t_min beyond the point where there are clients
    // waiting for reanimation as they would never succeed.
    Instant const allowable_desired_t_min =
        std::min(desired_t_min, reanimator_clientele_.first());
    must_restart =
        last_desired_t_min_.has_value() &&
        checkpointer_->checkpoint_at_or_before(last_desired_t_min_.value()) <
            checkpointer_->checkpoint_at_or_before(allowable_desired_t_min);
    LOG_IF(WARNING, must_restart)
        << "Restarting reanimator because desired t_min went from "
        << last_desired_t_min_.value() << " to " << allowable_desired_t_min;
    last_desired_t_min_ = allowable_desired_t_min;
  }

  if (must_restart) {
    reanimator_.Restart();
  }

  {
    absl::MutexLock l(&lock_);
    if (DesiredTMinReachedOrFullyReanimated(last_desired_t_min_.value())) {
      return;
    }
  }

  reanimator_.Put(last_desired_t_min_.value());
}

void Vessel::AwaitReanimation(Instant const& desired_t_min) {
  auto desired_t_min_reached_or_fully_reanimated = [this, desired_t_min]() {
    return DesiredTMinReachedOrFullyReanimated(desired_t_min);
  };

  Client me(desired_t_min, reanimator_clientele_);
  RequestReanimation(desired_t_min);
  absl::ReaderMutexLock l(&lock_);
  lock_.Await(absl::Condition(&desired_t_min_reached_or_fully_reanimated));
}

void Vessel::CreateFlightPlan(
    Instant const& final_time,
    Mass const& initial_mass,
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        flight_plan_adaptive_step_parameters,
    Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters const&
        flight_plan_generalized_adaptive_step_parameters) {
  auto const flight_plan_start = backstory_->back();
  flight_plans_.emplace_back(make_not_null_unique<FlightPlan>(
      initial_mass,
      /*initial_time=*/flight_plan_start.time,
      /*initial_degrees_of_freedom=*/flight_plan_start.degrees_of_freedom,
      final_time,
      ephemeris_,
      flight_plan_adaptive_step_parameters,
      flight_plan_generalized_adaptive_step_parameters));
  selected_flight_plan_index_ = flight_plans_.size() - 1;
}

void Vessel::DuplicateFlightPlan() {
  auto const& original = selected_flight_plan();
  auto const it = flight_plans_.begin() + selected_flight_plan_index_;
  // There is no need to adjust either flight plan after the copy; in particular
  // flight plans currently do not have a name, so there is no need to append a
  // (Copy).
  // If we needed to adjust something, it would probably be the original, to
  // mimic the behaviour of adding a copy after the current tab and switching to
  // that copy, even though behind the scenes we are inserting a copy before for
  // the sake of laziness.
  if (std::holds_alternative<serialization::FlightPlan>(original)) {
    flight_plans_.emplace(it, std::get<serialization::FlightPlan>(original));
  } else if (std::holds_alternative<not_null<std::unique_ptr<FlightPlan>>>(
                 original)) {
    std::get<not_null<std::unique_ptr<FlightPlan>>>(original)->WriteToMessage(
        &std::get<serialization::FlightPlan>(*flight_plans_.emplace(
            it, std::in_place_type<serialization::FlightPlan>)));
  } else {
    LOG(FATAL) << "Unexpected flight plan variant " << original.index();
  }
  ++selected_flight_plan_index_;
}

void Vessel::DeleteFlightPlan() {
  flight_plans_.erase(flight_plans_.begin() + selected_flight_plan_index_);
  if (selected_flight_plan_index_ == flight_plans_.size()) {
    --selected_flight_plan_index_;
  }
}

absl::Status Vessel::RebaseFlightPlan(Mass const& initial_mass) {
  CHECK(has_deserialized_flight_plan());
  auto& flight_plan =
      std::get<not_null<std::unique_ptr<FlightPlan>>>(selected_flight_plan());
  Instant const new_initial_time = backstory_->back().time;
  int first_manœuvre_kept = 0;
  for (int i = 0; i < flight_plan->number_of_manœuvres(); ++i) {
    auto const& manœuvre = flight_plan->GetManœuvre(i);
    if (manœuvre.initial_time() < new_initial_time) {
      first_manœuvre_kept = i + 1;
      if (new_initial_time < manœuvre.final_time()) {
        return absl::UnavailableError(
            "Cannot rebase during planned manœuvre execution");
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
  flight_plan = make_not_null_unique<FlightPlan>(
      initial_mass,
      /*initial_time=*/new_initial_time,
      /*initial_degrees_of_freedom=*/backstory_->back().degrees_of_freedom,
      new_desired_final_time,
      ephemeris_,
      original_flight_plan->adaptive_step_parameters(),
      original_flight_plan->generalized_adaptive_step_parameters());
  for (int i = first_manœuvre_kept;
       i < original_flight_plan->number_of_manœuvres();
       ++i) {
    auto const& manœuvre = original_flight_plan->GetManœuvre(i);
    flight_plan->Insert(manœuvre.burn(), i - first_manœuvre_kept).IgnoreError();
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
  if (downsampling_parameters_.has_value()) {
    auto* const serialized_downsampling_parameters =
        message->mutable_downsampling_parameters();
    serialized_downsampling_parameters->set_max_dense_intervals(
        downsampling_parameters_->max_dense_intervals);
    downsampling_parameters_->tolerance.WriteToMessage(
        serialized_downsampling_parameters->mutable_tolerance());
  }
  for (auto const& [_, part] : parts_) {
    part->WriteToMessage(message->add_parts(), serialization_index_for_pile_up);
  }
  for (auto const& part_id : kept_parts_) {
    CHECK(Contains(parts_, part_id));
    message->add_kept_parts(part_id);
  }

  // If the vessel is collapsible, we serialize at most the last
  // |max_points_to_serialize| of the part of the trajectory that ends at the
  // |backstory_|.  If it is not, however, we must serialize at least the entire
  // |backstory_| otherwise we'd lose the beginning of a non-collapsible
  // segment.
  std::int64_t const history_size = backstory_->end() - trajectory_.begin();
  std::int64_t const max_points_to_serialize_present_in_history =
      std::min(max_points_to_serialize, history_size);
  std::int64_t const serialized_points =
      is_collapsible_ ? max_points_to_serialize_present_in_history
                      : std::max(max_points_to_serialize_present_in_history,
                                 backstory_->size());

  // Starting with Gateaux we don't save the prediction, see #2685.  Instead we
  // just save its first point and re-read as if it was the whole prediction.
  trajectory_.WriteToMessage(
      message->mutable_history(),
      /*begin=*/backstory_->end() - serialized_points,
      /*end=*/std::next(prediction_->begin()),
      /*tracked=*/{backstory_, psychohistory_, prediction_},
      /*exact=*/{});
  for (auto const& flight_plan : flight_plans_) {
    if (std::holds_alternative<serialization::FlightPlan>(flight_plan)) {
      *message->add_flight_plans() =
          std::get<serialization::FlightPlan>(flight_plan);
    } else if (std::holds_alternative<not_null<std::unique_ptr<FlightPlan>>>(
                   flight_plan)) {
      auto& deserialized_flight_plan =
          std::get<not_null<std::unique_ptr<FlightPlan>>>(flight_plan);
      deserialized_flight_plan->WriteToMessage(message->add_flight_plans());
    } else {
      LOG(FATAL) << "Unexpected flight plan variant " << flight_plan.index();
    }
  }
  message->set_selected_flight_plan_index(selected_flight_plan_index_);
  message->set_is_collapsible(is_collapsible_);
  checkpointer_->WriteToMessage(message->mutable_checkpoint());
  LOG(INFO) << name_ << " " << NAMED(message->SpaceUsed()) << " "
            << NAMED(message->ByteSize());
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
  bool const is_pre_हरीश_चंद्र = !message.has_is_collapsible();
  bool const is_pre_hilbert = !message.has_selected_flight_plan_index();
  LOG_IF(WARNING, is_pre_hilbert)
      << "Reading pre-"
      << (is_pre_cesàro     ? "Cesàro"
          : is_pre_chasles  ? "Chasles"
          : is_pre_陈景润    ? "陈景润"
          : is_pre_hamilton ? "Hamilton"
          : is_pre_हरीश_चंद्र  ? "हरीश चंद्र"
                            : "Hilbert") << " Vessel";

  // NOTE(egg): for now we do not read the |MasslessBody| as it can contain no
  // information.
  auto vessel = make_not_null_unique<Vessel>(
      message.guid(),
      message.name(),
      parent,
      ephemeris,
      Ephemeris<Barycentric>::AdaptiveStepParameters::ReadFromMessage(
          message.prediction_adaptive_step_parameters()),
      DefaultDownsamplingParameters());
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
                                                         /*tracked=*/{});
    // The |backstory_| has been created by the constructor above.  Reconstruct
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
    vessel->backstory_ = std::prev(vessel->psychohistory_);
    vessel->prediction_ = vessel->trajectory_.NewSegment();
    vessel->downsampling_parameters_ = DefaultDownsamplingParameters();
  } else if (is_pre_chasles) {
    vessel->trajectory_ = DiscreteTrajectory<Barycentric>::ReadFromMessage(
        message.history(),
        /*tracked=*/{&vessel->psychohistory_});
    vessel->backstory_ = vessel->trajectory_.segments().begin();
    CHECK(vessel->backstory_ == std::prev(vessel->psychohistory_));
    vessel->prediction_ = vessel->trajectory_.NewSegment();
    vessel->downsampling_parameters_ = DefaultDownsamplingParameters();
  } else if (is_pre_hamilton) {
    vessel->trajectory_ = DiscreteTrajectory<Barycentric>::ReadFromMessage(
        message.history(),
        /*tracked=*/{&vessel->psychohistory_, &vessel->prediction_});
    vessel->backstory_ = vessel->trajectory_.segments().begin();
    CHECK(vessel->backstory_ == std::prev(vessel->psychohistory_));
    vessel->downsampling_parameters_ = DefaultDownsamplingParameters();
  } else if (is_pre_हरीश_चंद्र) {
    DiscreteTrajectorySegmentIterator<Barycentric> history;
    vessel->trajectory_ = DiscreteTrajectory<Barycentric>::ReadFromMessage(
        message.history(),
        /*tracked=*/{&history,
                     &vessel->psychohistory_,
                     &vessel->prediction_});
    vessel->backstory_ = std::prev(vessel->psychohistory_);
    vessel->downsampling_parameters_ = DefaultDownsamplingParameters();
  } else {
    vessel->trajectory_ = DiscreteTrajectory<Barycentric>::ReadFromMessage(
        message.history(),
        /*tracked=*/{&vessel->backstory_,
                     &vessel->psychohistory_,
                     &vessel->prediction_});
    vessel->is_collapsible_ = message.is_collapsible();

    vessel->checkpointer_ =
        Checkpointer<serialization::Vessel>::ReadFromMessage(
            vessel->MakeCheckpointerWriter(),
            vessel->MakeCheckpointerReader(),
            message.checkpoint());
    if (message.has_downsampling_parameters()) {
      vessel->downsampling_parameters_ =
          DiscreteTrajectorySegment<Barycentric>::DownsamplingParameters{
              .max_dense_intervals =
                  message.downsampling_parameters().max_dense_intervals(),
              .tolerance = Length::ReadFromMessage(
                  message.downsampling_parameters().tolerance())};
    } else {
      vessel->downsampling_parameters_ = std::nullopt;
    }
  }

  // Necessary after Εὔδοξος because the ephemeris has not been prolonged
  // during deserialization.
  ephemeris->Prolong(vessel->prediction_->back().time).IgnoreError();

  if (is_pre_陈景润) {
    vessel->backstory_->SetDownsamplingUnconditionally(
        DefaultDownsamplingParameters());
  }

  for (const auto& flight_plan : message.flight_plans()) {
    // Starting with हरीश चंद्र we deserialize the flight plan lazily.
    vessel->flight_plans_.emplace_back(flight_plan);
  }
  vessel->selected_flight_plan_index_ =
      is_pre_hilbert ? static_cast<int>(vessel->flight_plans_.size()) - 1
                     : message.selected_flight_plan_index();

  // Figure out which was the last checkpoint to be "reanimated" by reading the
  // end of the trajectory from the serialized form.  Interestingly enough, that
  // checkpoint (that is, the non-collapsible segment) may overlap the beginning
  // of the trajectory that we just deserialized, in which case we must rebuild
  // the front part of the non-collapsible segment to make sure that the
  // trajectory doesn't start in the middle of a non-collapsible segment (the
  // integration of the preceding collapsible segment would not end at the right
  // time if it did).
  Instant const checkpoint =
      vessel->checkpointer_->checkpoint_at_or_after(
          vessel->trajectory().t_min());
  if (checkpoint != InfiniteFuture) {
    CHECK_OK(vessel->checkpointer_->ReadFromCheckpointAt(
        checkpoint,
        [checkpoint, &vessel](
            serialization::Vessel::Checkpoint const& message) {
          // This code is similar to the one in ReanimateOneCheckpoint except
          // that (1) we never need to reconstruct a collapsible segment; (2) we
          // may actually have to truncate the non-collapsible segment obtained
          // from the checkpoint.
          LOG(INFO) << "Restoring " << vessel->ShortDebugString()
                    << " to initial checkpoint at " << checkpoint;

          DiscreteTrajectorySegmentIterator<Barycentric> unused;
          auto reanimated_trajectory =
              DiscreteTrajectory<Barycentric>::ReadFromMessage(
                  message.non_collapsible_segment(),
                  /*tracked=*/{&unused});
          CHECK(!reanimated_trajectory.empty());
          CHECK_EQ(checkpoint, reanimated_trajectory.back().time);
          reanimated_trajectory.ForgetAfter(vessel->trajectory().t_min());
          if (!reanimated_trajectory.empty()) {
            vessel->trajectory_.Merge(std::move(reanimated_trajectory));
          }
          return absl::OkStatus();
        }));
  }
  vessel->oldest_reanimated_checkpoint_ = checkpoint;

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
      checkpointer_(make_not_null_unique<Checkpointer<serialization::Vessel>>(
          /*reader=*/nullptr,
          /*writer=*/nullptr)),
      reanimator_(/*action=*/nullptr, 0ms),
      reanimator_clientele_(InfiniteFuture),
      backstory_(trajectory_.segments().begin()),
      psychohistory_(trajectory_.segments().end()),
      prediction_(trajectory_.segments().end()),
      prognosticator_(nullptr, 20ms) {}

Checkpointer<serialization::Vessel>::Writer Vessel::MakeCheckpointerWriter() {
  return [this](not_null<serialization::Vessel::Checkpoint*> const message) {
    // The extremities of the |backstory_| are implicitly exact.  Note that
    // |backstory_->end()| might cause serialization of a 1-point psychohistory
    // or prediction (at the last time of the backstory).  To figure things out
    // when reading we must track the |backstory_|.
    trajectory_.WriteToMessage(message->mutable_non_collapsible_segment(),
                               backstory_->begin(),
                               backstory_->end(),
                               /*tracked=*/{backstory_},
                               /*exact=*/{});

    // Here the containing pile-up is the one for the collapsible segment.
    ForSomePart([message](Part& first_part) {
      first_part.containing_pile_up()->fixed_step_parameters().WriteToMessage(
          message->mutable_collapsible_fixed_step_parameters());
    });
  };
}

Checkpointer<serialization::Vessel>::Reader Vessel::MakeCheckpointerReader() {
  return [](serialization::Vessel::Checkpoint const& message) {
    return absl::OkStatus();
  };
}

absl::Status Vessel::Reanimate(Instant const desired_t_min) {
  // This method is very similar to Ephemeris::Reanimate.  See the comments
  // there for some of the subtle points.
  static_assert(base::is_serializable_v<Barycentric>);
  absl::btree_set<Instant> checkpoints;
  LOG(INFO) << "Reanimating " << ShortDebugString() << " until "
            << desired_t_min;

  Instant t_final;
  {
    absl::ReaderMutexLock l(&lock_);
    if (reanimated_trajectories_.empty()) {
      t_final = trajectory_.begin()->time;
    } else {
      t_final = reanimated_trajectories_.back().front().time;
    }

    Instant const oldest_checkpoint_to_reanimate =
        checkpointer_->checkpoint_at_or_before(desired_t_min);
    checkpoints = checkpointer_->all_checkpoints_between(
        oldest_checkpoint_to_reanimate, oldest_reanimated_checkpoint_);

    // The |oldest_reanimated_checkpoint_| has already been reanimated before,
    // we don't need it below.
    checkpoints.erase(oldest_reanimated_checkpoint_);
  }

  for (auto it = checkpoints.crbegin(); it != checkpoints.crend(); ++it) {
    Instant const& checkpoint = *it;
    RETURN_IF_ERROR(checkpointer_->ReadFromCheckpointAt(
        checkpoint,
        [this, t_initial = checkpoint, &t_final](
            serialization::Vessel::Checkpoint const& message) -> absl::Status {
          auto const status_or_t_final =
              ReanimateOneCheckpoint(message, t_initial, t_final);
          RETURN_IF_ERROR(status_or_t_final);
          t_final = status_or_t_final.value();
          return absl::OkStatus();
        }));
  }
  return absl::OkStatus();
}

absl::StatusOr<Instant> Vessel::ReanimateOneCheckpoint(
    serialization::Vessel::Checkpoint const& message,
    Instant const& t_initial,
    Instant const& t_final) {
  CHECK_LE(t_initial, t_final);
  LOG(INFO) << "Restoring " << ShortDebugString() << " to checkpoint at "
            << t_initial << " until " << t_final;

  // Restore the non-collapsible segment that was fully saved.  It was the
  // backstory when the checkpoint was taken.
  DiscreteTrajectorySegmentIterator<Barycentric> reanimated_backstory;
  auto reanimated_trajectory =
      DiscreteTrajectory<Barycentric>::ReadFromMessage(
          message.non_collapsible_segment(),
          /*tracked=*/{&reanimated_backstory});
  auto const collapsible_fixed_step_parameters =
      Ephemeris<Barycentric>::FixedStepParameters::ReadFromMessage(
          message.collapsible_fixed_step_parameters());
  CHECK(!reanimated_trajectory.empty());
  CHECK_EQ(t_initial, reanimated_trajectory.back().time);
  Instant const reanimated_trajectory_t_initial =
      reanimated_trajectory.front().time;
  std::int64_t const reanimated_trajectory_size = reanimated_trajectory.size();

  // Construct a new collapsible segment at the end of the non-collapsible
  // backstory and integrate it until |t_final|.
  ++reanimated_backstory;
  reanimated_trajectory.DeleteSegments(reanimated_backstory);
  auto const collapsible_segment = reanimated_trajectory.NewSegment();

  // Make sure that the ephemeris covers the times that we are going to
  // reanimate.
  ephemeris_->AwaitReanimation(t_initial);
  auto fixed_instance =
      ephemeris_->NewInstance({&reanimated_trajectory},
                              Ephemeris<Barycentric>::NoIntrinsicAccelerations,
                              collapsible_fixed_step_parameters);

  auto const status = ephemeris_->FlowWithFixedStep(t_final, *fixed_instance);
  RETURN_IF_ERROR(status);

  LOG(INFO) << "Burn from " << reanimated_trajectory.front().time << " to "
            << t_initial << " (" << reanimated_trajectory_size
            << " points), coast to " << t_final << " ("
            << collapsible_segment->size() << " points)";


  // Push the reanimated trajectory into the queue where it will be consumed by
  // RequestReanimation.
  {
    absl::MutexLock l(&lock_);
    reanimated_trajectories_.push(std::move(reanimated_trajectory));
    oldest_reanimated_checkpoint_ = t_initial;
  }

  return reanimated_trajectory_t_initial;
}

bool Vessel::DesiredTMinReachedOrFullyReanimated(
    Instant const& desired_t_min) {
  lock_.AssertReaderHeld();

  // Consume the reanimated trajectories and merge them into this trajectory.
  // This is the only place where the reanimation becomes externally visible,
  // thereby ensuring that the trajectory doesn't change, say, while clients
  // iterate over it.
  while (!reanimated_trajectories_.empty()) {
    trajectory_.Merge(std::move(reanimated_trajectories_.front()));
    reanimated_trajectories_.pop();
  }
  return trajectory_.t_min() <= desired_t_min ||
         oldest_reanimated_checkpoint_ == checkpointer_->oldest_checkpoint();
}

absl::StatusOr<DiscreteTrajectory<Barycentric>> Vessel::FlowPrognostication(
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

bool Vessel::IsCollapsible() const {
  PileUp* containing_pile_up = nullptr;
  std::set<not_null<Part*>> parts;
  for (const auto& [_, part] : parts_) {
    // We expect parts to be piled up.
    CHECK(part->is_piled_up());
    // Not collapsible if any part has a force applied to it (but a torque is
    // fine).
    if (part->intrinsic_force() != Vector<Force, Barycentric>{}) {
      return false;
    }
    parts.insert(part.get());
    if (containing_pile_up == nullptr) {
      containing_pile_up = part->containing_pile_up();
    } else {
      // All the parts should be in the same pile-up.
      CHECK_EQ(containing_pile_up, part->containing_pile_up());
    }
  }
  CHECK_NE(nullptr, containing_pile_up);
  for (const auto part : containing_pile_up->parts()) {
    // Not collapsible if the pile-up contains a part not in this vessel.
    if (!parts.contains(part)) {
      return false;
    }
  }
  return true;
}

bool Vessel::has_deserialized_flight_plan() const {
  return !flight_plans_.empty() &&
         std::holds_alternative<not_null<std::unique_ptr<FlightPlan>>>(
             selected_flight_plan());
}

Vessel::LazilyDeserializedFlightPlan& Vessel::selected_flight_plan() {
  CHECK_GE(selected_flight_plan_index_, 0);
  CHECK_LT(selected_flight_plan_index_, flight_plans_.size());
  return flight_plans_[selected_flight_plan_index_];
}

Vessel::LazilyDeserializedFlightPlan const&
Vessel::selected_flight_plan() const {
  CHECK_GE(selected_flight_plan_index_, 0);
  CHECK_LT(selected_flight_plan_index_, flight_plans_.size());
  return flight_plans_[selected_flight_plan_index_];
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
