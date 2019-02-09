
#include "ksp_plugin/vessel.hpp"

#include <algorithm>
#include <limits>
#include <list>
#include <string>
#include <vector>

#include "astronomy/epoch.hpp"
#include "base/map_util.hpp"
#include "ksp_plugin/integrators.hpp"
#include "ksp_plugin/pile_up.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/make_not_null.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_vessel {

using astronomy::InfiniteFuture;
using base::check_not_null;
using base::Contains;
using base::FindOrDie;
using base::make_not_null_unique;
using geometry::BarycentreCalculator;
using geometry::Position;
using quantities::IsFinite;
using quantities::Length;
using quantities::Time;
using quantities::si::Metre;

constexpr std::int64_t max_dense_intervals = 10'000;
constexpr Length downsampling_tolerance = 10 * Metre;

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
      history_(make_not_null_unique<DiscreteTrajectory<Barycentric>>()) {
  // Can't create the |psychohistory_| here because |history_| is empty;
}

Vessel::~Vessel() {
  LOG(INFO) << "Destroying vessel " << ShortDebugString();
  // Ask the prognosticator to shut down.  This may take a while.  Make sure
  // that we handle the case where |PrepareHistory| was not called.
  if (prognosticator_.joinable()) {
    {
      absl::MutexLock l(&prognosticator_lock_);
      CHECK(prognosticator_parameters_);
      prognosticator_parameters_->shutdown = true;
    }
    prognosticator_.join();
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

void Vessel::ClearAllIntrinsicForces() {
  for (auto const& pair : parts_) {
    auto const& part = pair.second;
    part->clear_intrinsic_force();
  }
}

void Vessel::PrepareHistory(Instant const& t) {
  CHECK(!parts_.empty());
  if (history_->Empty()) {
    LOG(INFO) << "Preparing history of vessel " << ShortDebugString()
              << " at " << t;
    BarycentreCalculator<DegreesOfFreedom<Barycentric>, Mass> calculator;
    ForAllParts([&calculator](Part& part) {
      calculator.Add(part.degrees_of_freedom(), part.mass());
    });
    CHECK(psychohistory_ == nullptr);
    history_->SetDownsampling(max_dense_intervals, downsampling_tolerance);
    history_->Append(t, calculator.Get());
    psychohistory_ = history_->NewForkAtLast();
    prediction_ = psychohistory_->NewForkAtLast();

    // Prepare the parameters for the |prognosticator_| and start it.
    {
      absl::MutexLock l(&prognosticator_lock_);
      prognosticator_parameters_ =
          PrognosticatorParameters{psychohistory_->last().time(),
                                   psychohistory_->last().degrees_of_freedom(),
                                   /*last_time=*/std::nullopt,
                                   prediction_adaptive_step_parameters_,
                                   /*shutdown=*/false};
    }
    prognosticator_ =
        std::thread(std::bind(&Vessel::RepeatedlyFlowPrognostication, this));
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
  for (auto const& pair : parts_) {
    action(*pair.second);
  }
}

DiscreteTrajectory<Barycentric> const& Vessel::psychohistory() const {
  return *psychohistory_;
}

DiscreteTrajectory<Barycentric> const& Vessel::prediction() const {
  return *prediction_;
}

void Vessel::set_prediction_adaptive_step_parameters(
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        prediction_adaptive_step_parameters) {
  prediction_adaptive_step_parameters_ = prediction_adaptive_step_parameters;
  absl::MutexLock l(&prognosticator_lock_);
  prognosticator_parameters_->adaptive_step_parameters =
      prediction_adaptive_step_parameters;
}

Ephemeris<Barycentric>::AdaptiveStepParameters const&
Vessel::prediction_adaptive_step_parameters() const {
  absl::ReaderMutexLock l(&prognosticator_lock_);
  return prognosticator_parameters_->adaptive_step_parameters;
}

FlightPlan& Vessel::flight_plan() const {
  CHECK(has_flight_plan());
  return *flight_plan_;
}

bool Vessel::has_flight_plan() const {
  return flight_plan_ != nullptr;
}

void Vessel::AdvanceTime() {
  // Squirrel away the prediction so that we can reattach it if we don't have a
  // prognostication.
  auto prediction = prediction_->DetachFork();

  history_->DeleteFork(psychohistory_);
  AppendToVesselTrajectory(&Part::history_begin,
                           &Part::history_end,
                           *history_);
  psychohistory_ = history_->NewForkAtLast();
  AppendToVesselTrajectory(&Part::psychohistory_begin,
                           &Part::psychohistory_end,
                           *psychohistory_);
  {
    absl::MutexLock l(&prognosticator_lock_);
    if (prognostication_ == nullptr) {
      AttachPrediction(std::move(prediction));
    } else {
      AttachPrediction(std::move(prognostication_));
    }
  }

  for (auto const& pair : parts_) {
    Part& part = *pair.second;
    part.ClearHistory();
  }
}

void Vessel::ForgetBefore(Instant const& time) {
  // Make sure that the history keeps at least one point and don't change the
  // psychohistory or prediction.  We cannot use the parts because they may have
  // been moved to the future already.
  history_->ForgetBefore(std::min(time, history_->last().time()));
  if (flight_plan_ != nullptr) {
    flight_plan_->ForgetBefore(time, [this]() { flight_plan_.reset(); });
  }
}

void Vessel::CreateFlightPlan(
    Instant const& final_time,
    Mass const& initial_mass,
    Ephemeris<Barycentric>::AdaptiveStepParameters const&
        flight_plan_adaptive_step_parameters,
    Ephemeris<Barycentric>::GeneralizedAdaptiveStepParameters const&
        flight_plan_generalized_adaptive_step_parameters) {
  auto const history_last = history_->last();
  flight_plan_ = std::make_unique<FlightPlan>(
      initial_mass,
      /*initial_time=*/history_last.time(),
      /*initial_degrees_of_freedom=*/history_last.degrees_of_freedom(),
      final_time,
      ephemeris_,
      flight_plan_adaptive_step_parameters,
      flight_plan_generalized_adaptive_step_parameters);
}

void Vessel::DeleteFlightPlan() {
  flight_plan_.reset();
}

void Vessel::FlowPrediction() {
  absl::MutexLock l(&prognosticator_lock_);
  prognosticator_parameters_->last_time = std::nullopt;
  if (prognostication_ != nullptr) {
    AttachPrediction(std::move(prognostication_));
  }
}

void Vessel::FlowPrediction(Instant const& time) {
  absl::MutexLock l(&prognosticator_lock_);
  prognosticator_parameters_->last_time = time;
  if (prognostication_ != nullptr) {
    AttachPrediction(std::move(prognostication_));
  }
}

void Vessel::WriteToMessage(not_null<serialization::Vessel*> const message,
                            PileUp::SerializationIndexForPileUp const&
                                serialization_index_for_pile_up) const {
  message->set_guid(guid_);
  message->set_name(name_);
  body_.WriteToMessage(message->mutable_body());
  prediction_adaptive_step_parameters_.WriteToMessage(
      message->mutable_prediction_adaptive_step_parameters());
  for (auto const& pair : parts_) {
    auto const& part = pair.second;
    part->WriteToMessage(message->add_parts(), serialization_index_for_pile_up);
  }
  for (auto const& part_id : kept_parts_) {
    CHECK(Contains(parts_, part_id));
    message->add_kept_parts(part_id);
  }
  {
    absl::ReaderMutexLock l(&prognosticator_lock_);
    history_->WriteToMessage(message->mutable_history(),
                             /*forks=*/{psychohistory_, prediction_});
  }
  if (flight_plan_ != nullptr) {
    flight_plan_->WriteToMessage(message->mutable_flight_plan());
  }
}

not_null<std::unique_ptr<Vessel>> Vessel::ReadFromMessage(
    serialization::Vessel const& message,
    not_null<Celestial const*> const parent,
    not_null<Ephemeris<Barycentric>*> const ephemeris,
    std::function<void(PartId)> const& deletion_callback) {
  bool const is_pre_cesàro = message.has_psychohistory_is_authoritative();
  bool const is_pre_chasles = message.has_prediction();
  bool const is_pre_陈景润 = !message.history().has_downsampling();

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
    for (auto it = psychohistory->Begin(); it != psychohistory->End(); ++it) {
      if (it == psychohistory->last() &&
          !message.psychohistory_is_authoritative()) {
        vessel->psychohistory_ = vessel->history_->NewForkAtLast();
        vessel->psychohistory_->Append(it.time(), it.degrees_of_freedom());
      } else {
        vessel->history_->Append(it.time(), it.degrees_of_freedom());
      }
    }
    if (message.psychohistory_is_authoritative()) {
      vessel->psychohistory_ = vessel->history_->NewForkAtLast();
    }
    vessel->prediction_ = vessel->psychohistory_->NewForkAtLast();
    vessel->FlowPrediction();
  } else if (is_pre_chasles) {
    vessel->history_ = DiscreteTrajectory<Barycentric>::ReadFromMessage(
        message.history(),
        /*forks=*/{&vessel->psychohistory_});
    vessel->prediction_ = vessel->psychohistory_->NewForkAtLast();
    vessel->FlowPrediction();
  } else {
    vessel->history_ = DiscreteTrajectory<Barycentric>::ReadFromMessage(
        message.history(),
        /*forks=*/{&vessel->psychohistory_, &vessel->prediction_});
    // Necessary after Εὔδοξος because the ephemeris has not been prolonged
    // during deserialization.  Doesn't hurt prior to Εὔδοξος.
    ephemeris->Prolong(vessel->prediction_->last().time());
  }

  if (is_pre_陈景润) {
    vessel->history_->SetDownsampling(max_dense_intervals,
                                      downsampling_tolerance);
  }

  if (message.has_flight_plan()) {
    vessel->flight_plan_ = FlightPlan::ReadFromMessage(message.flight_plan(),
                                                       ephemeris);
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

std::string Vessel::ShortDebugString() const {
  return name_ + " (" + guid_ + ")";
}

Vessel::Vessel()
    : body_(),
      prediction_adaptive_step_parameters_(DefaultPredictionParameters()),
      parent_(testing_utilities::make_not_null<Celestial const*>()),
      ephemeris_(testing_utilities::make_not_null<Ephemeris<Barycentric>*>()),
      history_(make_not_null_unique<DiscreteTrajectory<Barycentric>>()) {}

void Vessel::RepeatedlyFlowPrognostication() {
  for (;;) {
    // No point in going faster than 50 Hz.
    std::chrono::steady_clock::time_point const wakeup_time =
        std::chrono::steady_clock::now() + std::chrono::milliseconds(20);

    // The thread is only started after the parameters have been set, so we
    // should always find parameters here.
    std::optional<PrognosticatorParameters> prognosticator_parameters;
    {
      absl::ReaderMutexLock l(&prognosticator_lock_);
      CHECK(prognosticator_parameters_);
      prognosticator_parameters = prognosticator_parameters_;
    }

    if (prognosticator_parameters->shutdown) {
      break;
    }

    auto prognostication = std::make_unique<DiscreteTrajectory<Barycentric>>();
    prognostication->Append(
        prognosticator_parameters->first_time,
        prognosticator_parameters->first_degrees_of_freedom);
    if (prognosticator_parameters->last_time) {
      if (*prognosticator_parameters->last_time >
          prognostication->last().time()) {
        ephemeris_->FlowWithAdaptiveStep(
            prognostication.get(),
            Ephemeris<Barycentric>::NoIntrinsicAcceleration,
            *prognosticator_parameters->last_time,
            prognosticator_parameters->adaptive_step_parameters,
            FlightPlan::max_ephemeris_steps_per_frame,
            /*last_point_only=*/false);
      }
    } else {
      bool const reached_t = ephemeris_->FlowWithAdaptiveStep(
          prognostication.get(),
          Ephemeris<Barycentric>::NoIntrinsicAcceleration,
          ephemeris_->t_max(),
          prognosticator_parameters->adaptive_step_parameters,
          FlightPlan::max_ephemeris_steps_per_frame,
          /*last_point_only=*/false).ok();
      if (reached_t) {
        // This will prolong the ephemeris by |max_ephemeris_steps_per_frame|.
        ephemeris_->FlowWithAdaptiveStep(
          prognostication.get(),
          Ephemeris<Barycentric>::NoIntrinsicAcceleration,
          InfiniteFuture,
          prognosticator_parameters->adaptive_step_parameters,
          FlightPlan::max_ephemeris_steps_per_frame,
          /*last_point_only=*/false);
      }
    }

    // Publish the prognostication.
    {
      absl::MutexLock l(&prognosticator_lock_);
      prognostication_.swap(prognostication);
    }

    std::this_thread::sleep_until(wakeup_time);
  }
}

void Vessel::AppendToVesselTrajectory(
    TrajectoryIterator const part_trajectory_begin,
    TrajectoryIterator const part_trajectory_end,
    DiscreteTrajectory<Barycentric>& trajectory) {
  CHECK(!parts_.empty());
  std::vector<DiscreteTrajectory<Barycentric>::Iterator> its;
  std::vector<DiscreteTrajectory<Barycentric>::Iterator> ends;
  its.reserve(parts_.size());
  ends.reserve(parts_.size());
  for (auto const& pair : parts_) {
    Part& part = *pair.second;
    its.push_back((part.*part_trajectory_begin)());
    ends.push_back((part.*part_trajectory_end)());
  }

  // Loop over the times of the trajectory.
  for (;;) {
    auto const& it0 = its[0];
    bool const at_end_of_part_trajectory = it0 == ends[0];
    Instant const first_time = at_end_of_part_trajectory ? Instant()
                                                         : it0.time();

    // Loop over the parts at a given time.
    BarycentreCalculator<DegreesOfFreedom<Barycentric>, Mass> calculator;
    int i = 0;
    for (auto const& pair : parts_) {
      Part& part = *pair.second;
      auto& it = its[i];
      CHECK_EQ(at_end_of_part_trajectory, it == ends[i]);
      if (!at_end_of_part_trajectory) {
        calculator.Add(it.degrees_of_freedom(), part.mass());
        CHECK_EQ(first_time, it.time());
        ++it;
      }
      ++i;
    }

    if (at_end_of_part_trajectory) {
      return;
    }

    // Append the parts' barycentre to the trajectory.
    DegreesOfFreedom<Barycentric> const vessel_degrees_of_freedom =
        calculator.Get();
    trajectory.Append(first_time, vessel_degrees_of_freedom);
  }
}

void Vessel::AttachPrediction(
    not_null<std::unique_ptr<DiscreteTrajectory<Barycentric>>> trajectory) {
  trajectory->ForgetBefore(psychohistory_->last().time());
  if (trajectory->Empty()) {
    prediction_ = psychohistory_->NewForkAtLast();
  } else {
    prediction_ = trajectory.get();
    psychohistory_->AttachFork(std::move(trajectory));
  }
}

}  // namespace internal_vessel
}  // namespace ksp_plugin
}  // namespace principia
