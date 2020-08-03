
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
using base::Contains;
using base::Error;
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
             right.adaptive_step_parameters.speed_integration_tolerance() ||
         left.shutdown != right.shutdown;
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
      history_(make_not_null_unique<DiscreteTrajectory<Barycentric>>()) {
  // Can't create the |psychohistory_| and |prediction_| here because |history_|
  // is empty;
}

Vessel::~Vessel() {
  LOG(INFO) << "Destroying vessel " << ShortDebugString();
  // Ask the prognosticator to shut down.  This may take a while.  Make sure
  // that we handle the case where |PrepareHistory| was not called.
  if (prognosticator_.joinable()) {
    {
      absl::MutexLock l(&prognosticator_lock_);
      prognosticator_parameters_ =
          PrognosticatorParameters{Ephemeris<Barycentric>::Guard(ephemeris_),
                                   psychohistory_->back().time,
                                   psychohistory_->back().degrees_of_freedom,
                                   prediction_adaptive_step_parameters_,
                                   /*shutdown=*/true};
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

void Vessel::ClearAllIntrinsicForcesAndTorques() {
  for (auto const& [_, part] : parts_) {
    part->clear_intrinsic_force();
    part->clear_intrinsic_torque();
  }
}

void Vessel::PrepareHistory(Instant const& t) {
  CHECK(!parts_.empty());
  if (history_->Empty()) {
    LOG(INFO) << "Preparing history of vessel " << ShortDebugString()
              << " at " << t;
    BarycentreCalculator<DegreesOfFreedom<Barycentric>, Mass> calculator;
    ForAllParts([&calculator](Part& part) {
      calculator.Add(
          part.rigid_motion()({RigidPart::origin, RigidPart::unmoving}),
          part.mass());
    });
    CHECK(psychohistory_ == nullptr);
    history_->SetDownsampling(max_dense_intervals, downsampling_tolerance);
    history_->Append(t, calculator.Get());
    psychohistory_ = history_->NewForkAtLast();
    prediction_ = psychohistory_->NewForkAtLast();
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
  if (prognosticator_parameters_) {
    prognosticator_parameters_->adaptive_step_parameters =
        prediction_adaptive_step_parameters;
  }
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
  // Squirrel away the prediction so that we can reattach it if we don't have a
  // prognostication.
  auto prediction = prediction_->DetachFork();

  history_->DeleteFork(psychohistory_);
  AppendToVesselTrajectory(&Part::history_begin,
                           &Part::history_end,
                           *history_);
  psychohistory_ = history_->NewForkAtLast();

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
  // point.  See #2507.
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

  for (auto const& [_, part] : parts_) {
    part->ClearHistory();
  }
}

void Vessel::ForgetBefore(Instant const& time) {
  // Make sure that the history keeps at least one point and don't change the
  // psychohistory or prediction.  We cannot use the parts because they may have
  // been moved to the future already.
  history_->ForgetBefore(std::min(time, history_->back().time));
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
  flight_plan_.reset();
}

Status Vessel::RebaseFlightPlan(Mass const& initial_mass) {
  CHECK(has_flight_plan());
  Instant const new_initial_time = history_->back().time;
  int first_manœuvre_kept = 0;
  for (int i = 0; i < flight_plan_->number_of_manœuvres(); ++i) {
    auto const& manœuvre = flight_plan_->GetManœuvre(i);
    if (manœuvre.initial_time() < new_initial_time) {
      first_manœuvre_kept = i + 1;
      if (new_initial_time < manœuvre.final_time()) {
        return Status(Error::UNAVAILABLE,
                      u8"Cannot rebase during planned manœuvre execution");
      }
    }
  }
  not_null<std::unique_ptr<FlightPlan>> const original_flight_plan =
      std::move(flight_plan_);
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
    flight_plan_->Insert(manœuvre.burn(), i - first_manœuvre_kept);
  }
  return Status::OK;
}

void Vessel::RefreshPrediction() {
  absl::MutexLock l(&prognosticator_lock_);
  // The guard below ensures that the ephemeris will not be "forgotten before"
  // the end of the psychohistory between now and the time when the
  // prognosticator finishes the integration.  This ensures that the ephemeris'
  // |t_min| is never after the time that |FlowWithAdaptiveStep| tries to
  // integrate.
  // The guard will be destroyed either when the next set of parameters is
  // created or when the prognostication has been computed.
  // Note that we know that both |EventuallyForgetBefore| and
  // |RefreshPrediction| are called on the main thread, therefore the ephemeris
  // currently covers the last time of the psychohistory.  Were this to change,
  // this code might have to change.
  prognosticator_parameters_ =
      PrognosticatorParameters{Ephemeris<Barycentric>::Guard(ephemeris_),
                               psychohistory_->back().time,
                               psychohistory_->back().degrees_of_freedom,
                               prediction_adaptive_step_parameters_,
                               /*shutdown=*/false};
  if (synchronous_) {
    std::unique_ptr<DiscreteTrajectory<Barycentric>> prognostication;
    std::optional<PrognosticatorParameters> prognosticator_parameters;
    std::swap(prognosticator_parameters, prognosticator_parameters_);
    Status const status =
        FlowPrognostication(std::move(*prognosticator_parameters),
                            prognostication);
    SwapPrognostication(prognostication, status);
  } else {
    StartPrognosticatorIfNeeded();
  }
  if (prognostication_ != nullptr) {
    AttachPrediction(std::move(prognostication_));
  }
}

void Vessel::RefreshPrediction(Instant const& time) {
  RefreshPrediction();
  prediction_->ForgetAfter(time);
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
  history_->WriteToMessage(message->mutable_history(),
                           /*forks=*/{psychohistory_, prediction_});
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
    for (auto it = psychohistory->begin(); it != psychohistory->end();) {
      auto const& [time, degrees_of_freedom] = *it;
      ++it;
      if (it == psychohistory->end() &&
          !message.psychohistory_is_authoritative()) {
        vessel->psychohistory_ = vessel->history_->NewForkAtLast();
        vessel->psychohistory_->Append(time, degrees_of_freedom);
      } else {
        vessel->history_->Append(time, degrees_of_freedom);
      }
    }
    if (message.psychohistory_is_authoritative()) {
      vessel->psychohistory_ = vessel->history_->NewForkAtLast();
    }
    vessel->prediction_ = vessel->psychohistory_->NewForkAtLast();
  } else if (is_pre_chasles) {
    vessel->history_ = DiscreteTrajectory<Barycentric>::ReadFromMessage(
        message.history(),
        /*forks=*/{&vessel->psychohistory_});
    vessel->prediction_ = vessel->psychohistory_->NewForkAtLast();
  } else {
    vessel->history_ = DiscreteTrajectory<Barycentric>::ReadFromMessage(
        message.history(),
        /*forks=*/{&vessel->psychohistory_, &vessel->prediction_});
    // Necessary after Εὔδοξος because the ephemeris has not been prolonged
    // during deserialization.  Doesn't hurt prior to Εὔδοξος.
    ephemeris->Prolong(vessel->prediction_->back().time);
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

void Vessel::RefreshOrbitAnalysis(Time const& mission_duration) {
  if (!orbit_analyser_.has_value()) {
    // TODO(egg): perhaps we should get the history parameters from the plugin;
    // on the other hand, these are probably overkill for high orbits anyway,
    // and given that we know many things about our trajectory in the analyser,
    // perhaps we should pick something appropriate automatically instead.  The
    // default will do in the meantime.
    orbit_analyser_.emplace(ephemeris_, DefaultHistoryParameters());
  }
  orbit_analyser_->RequestAnalysis(psychohistory_->back().time,
                                   psychohistory_->back().degrees_of_freedom,
                                   mission_duration);
  orbit_analyser_->RefreshAnalysis();
}

double Vessel::progress_of_orbit_analysis() const {
  if (!orbit_analyser_.has_value()) {
    return 0;
  }
  return orbit_analyser_->progress_of_next_analysis();
}

OrbitAnalyser::Analysis* Vessel::orbit_analysis() {
  return orbit_analyser_->analysis();
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
      history_(make_not_null_unique<DiscreteTrajectory<Barycentric>>()) {}

void Vessel::StartPrognosticatorIfNeeded() {
  prognosticator_lock_.AssertHeld();
  if (!prognosticator_.joinable()) {
    prognosticator_ =
        std::thread(std::bind(&Vessel::RepeatedlyFlowPrognostication, this));
  }
}

void Vessel::RepeatedlyFlowPrognostication() {
  for (;;) {
    // No point in going faster than 50 Hz.
    std::chrono::steady_clock::time_point const wakeup_time =
        std::chrono::steady_clock::now() + std::chrono::milliseconds(20);

    std::optional<PrognosticatorParameters> prognosticator_parameters;
    {
      absl::ReaderMutexLock l(&prognosticator_lock_);
      if (!prognosticator_parameters_) {
        // No parameters, let's wait for them to appear.
        continue;
      }
      std::swap(prognosticator_parameters, prognosticator_parameters_);
    }

    if (prognosticator_parameters->shutdown) {
      break;
    }

    std::unique_ptr<DiscreteTrajectory<Barycentric>> prognostication;
    Status const status =
        FlowPrognostication(std::move(*prognosticator_parameters),
                            prognostication);
    {
      absl::MutexLock l(&prognosticator_lock_);
      SwapPrognostication(prognostication, status);
    }

    std::this_thread::sleep_until(wakeup_time);
  }
}

Status Vessel::FlowPrognostication(
    PrognosticatorParameters prognosticator_parameters,
    std::unique_ptr<DiscreteTrajectory<Barycentric>>& prognostication) {
  // The guard contained in |prognosticator_parameters| ensures that the |t_min|
  // of the ephemeris doesn't move in this function.
  prognostication = std::make_unique<DiscreteTrajectory<Barycentric>>();
  prognostication->Append(
      prognosticator_parameters.first_time,
      prognosticator_parameters.first_degrees_of_freedom);
  Status status;
  status = ephemeris_->FlowWithAdaptiveStep(
      prognostication.get(),
      Ephemeris<Barycentric>::NoIntrinsicAcceleration,
      ephemeris_->t_max(),
      prognosticator_parameters.adaptive_step_parameters,
      FlightPlan::max_ephemeris_steps_per_frame);
  bool const reached_t_max = status.ok();
  if (reached_t_max) {
    // This will prolong the ephemeris by |max_ephemeris_steps_per_frame|.
    status = ephemeris_->FlowWithAdaptiveStep(
        prognostication.get(),
        Ephemeris<Barycentric>::NoIntrinsicAcceleration,
        InfiniteFuture,
        prognosticator_parameters.adaptive_step_parameters,
        FlightPlan::max_ephemeris_steps_per_frame);
  }
  LOG_IF_EVERY_N(INFO, !status.ok(), 50)
      << "Prognostication from " << prognosticator_parameters.first_time
      << " finished at " << prognostication->back().time << " with "
      << status.ToString() << " for " << ShortDebugString();
  return status;
}

void Vessel::SwapPrognostication(
    std::unique_ptr<DiscreteTrajectory<Barycentric>>& prognostication,
    Status const& status) {
  prognosticator_lock_.AssertHeld();
  if (status.error() != Error::CANCELLED) {
    prognostication_.swap(prognostication);
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
  for (auto const& [_, part] : parts_) {
    its.push_back((*part.*part_trajectory_begin)());
    ends.push_back((*part.*part_trajectory_end)());
  }

  Instant fork_time;
  if (!trajectory.is_root()) {
    fork_time = trajectory.Fork()->time;
  }

  // Loop over the times of the trajectory.
  for (;;) {
    auto const& it0 = its[0];
    bool const at_end_of_part_trajectory = it0 == ends[0];
    Instant const first_time = at_end_of_part_trajectory ? Instant()
                                                         : it0->time;
    bool const after_fork_time = !at_end_of_part_trajectory &&
                                 (trajectory.is_root() ||
                                  first_time > fork_time);

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
    if (after_fork_time) {
      DegreesOfFreedom<Barycentric> const vessel_degrees_of_freedom =
          calculator.Get();
      trajectory.Append(first_time, vessel_degrees_of_freedom);
    }
  }
}

void Vessel::AttachPrediction(
    not_null<std::unique_ptr<DiscreteTrajectory<Barycentric>>> trajectory) {
  trajectory->ForgetBefore(psychohistory_->back().time);
  if (trajectory->Empty()) {
    prediction_ = psychohistory_->NewForkAtLast();
  } else {
    prediction_ = trajectory.get();
    psychohistory_->AttachFork(std::move(trajectory));
  }
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
