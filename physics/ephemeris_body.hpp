#pragma once

#include "physics/ephemeris.hpp"

#include <algorithm>
#include <functional>
#include <limits>
#include <set>
#include <vector>

#include "base/macros.hpp"
#include "base/map_util.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "physics/continuous_trajectory.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {

using base::FindOrDie;
using geometry::InnerProduct;
using geometry::R3Element;
using integrators::AdaptiveStepSize;
using integrators::IntegrationProblem;
using quantities::Abs;
using quantities::Exponentiation;
using quantities::Quotient;
using quantities::Time;
using quantities::si::Day;
using ::std::placeholders::_1;
using ::std::placeholders::_2;
using ::std::placeholders::_3;

namespace physics {

namespace {  // TODO(egg): this should be a named namespace (internal)

Time const kMaxTimeBetweenIntermediateStates = 180 * Day;

// If j is a unit vector along the axis of rotation, and r is the separation
// between the bodies, the acceleration computed here is:
//
//   -(J2 / |r|^5) (3 j (r.j) + r (3 - 15 (r.j)^2 / |r|^2) / 2)
//
// Where |r| is the norm of r and r.j is the inner product.
template<typename Frame>
FORCE_INLINE Vector<Quotient<Acceleration,
                             GravitationalParameter>, Frame> Order2ZonalEffect(
    OblateBody<Frame> const& body,
    Displacement<Frame> const& r,
    Exponentiation<Length, -2> const& one_over_r_squared,
    Exponentiation<Length, -3> const& one_over_r_cubed) {
  Vector<double, Frame> const& axis = body.axis();
  Length const r_axis_projection = InnerProduct(axis, r);
  auto const j2_over_r_fifth =
      body.j2_over_μ() * one_over_r_cubed * one_over_r_squared;
  Vector<Quotient<Acceleration,
                  GravitationalParameter>, Frame> const axis_effect =
      (-3 * j2_over_r_fifth * r_axis_projection) * axis;
  Vector<Quotient<Acceleration,
                  GravitationalParameter>, Frame> const radial_effect =
      (j2_over_r_fifth *
           (-1.5 +
            7.5 * r_axis_projection *
                  r_axis_projection * one_over_r_squared)) * r;
  return axis_effect + radial_effect;
}

// For mocking purposes.
template<typename Frame>
class DummyIntegrator
    : public FixedStepSizeIntegrator<
                 typename Ephemeris<Frame>::NewtonianMotionEquation> {
  using ODE = typename Ephemeris<Frame>::NewtonianMotionEquation;
  DummyIntegrator() : FixedStepSizeIntegrator<ODE>(
      serialization::FixedStepSizeIntegrator::DUMMY) {}

 public:
  void Solve(IntegrationProblem<ODE> const& problem,
             Time const& step) const override {
    LOG(FATAL) << "dummy";
  }

  static DummyIntegrator const& Instance() {
    static DummyIntegrator const instance;
    return instance;
  }
};

}  // namespace

template<typename Frame>
Ephemeris<Frame>::Ephemeris(
    std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies,
    std::vector<DegreesOfFreedom<Frame>> const& initial_state,
    Instant const& initial_time,
    FixedStepSizeIntegrator<NewtonianMotionEquation> const&
        planetary_integrator,
    Time const& step,
    Length const& fitting_tolerance)
    : planetary_integrator_(planetary_integrator),
      step_(step),
      fitting_tolerance_(fitting_tolerance) {
  CHECK(!bodies.empty());
  CHECK_EQ(bodies.size(), initial_state.size());

  last_state_.time = initial_time;

  for (int i = 0; i < bodies.size(); ++i) {
    auto& body = bodies[i];
    DegreesOfFreedom<Frame> const& degrees_of_freedom = initial_state[i];

    unowned_bodies_.push_back(body.get());

    auto const inserted = bodies_to_trajectories_.emplace(
                              body.get(),
                              std::make_unique<ContinuousTrajectory<Frame>>(
                                  step_, fitting_tolerance_));
    CHECK(inserted.second);
    ContinuousTrajectory<Frame>* const trajectory =
        inserted.first->second.get();
    trajectory->Append(initial_time, degrees_of_freedom);

    VLOG(1) << "Constructed trajectory " << trajectory
            << " for body with mass " << body->mass();

    if (body->is_oblate()) {
      // Inserting at the beginning of the vectors is O(N).
      oblate_bodies_.insert(oblate_bodies_.begin(), body.get());
      bodies_.insert(bodies_.begin(), std::move(body));
      trajectories_.insert(trajectories_.begin(), trajectory);
      last_state_.positions.insert(last_state_.positions.begin(),
                                   degrees_of_freedom.position());
      last_state_.velocities.insert(last_state_.velocities.begin(),
                                    degrees_of_freedom.velocity());
      ++number_of_oblate_bodies_;
    } else {
      // Inserting at the end of the vectors is O(1).
      spherical_bodies_.push_back(body.get());
      bodies_.push_back(std::move(body));
      trajectories_.push_back(trajectory);
      last_state_.positions.push_back(degrees_of_freedom.position());
      last_state_.velocities.push_back(degrees_of_freedom.velocity());
      ++number_of_spherical_bodies_;
    }
  }

  massive_bodies_equation_.compute_acceleration =
      std::bind(&Ephemeris::ComputeMassiveBodiesGravitationalAccelerations,
                this, _1, _2, _3);
}

template<typename Frame>
std::vector<MassiveBody const*> const& Ephemeris<Frame>::bodies() const {
  return unowned_bodies_;
}

template<typename Frame>
not_null<ContinuousTrajectory<Frame> const*> Ephemeris<Frame>::trajectory(
    not_null<MassiveBody const*> body) const {
  return FindOrDie(bodies_to_trajectories_, body).get();
}

template<typename Frame>
bool Ephemeris<Frame>::empty() const {
  for (auto const& pair : bodies_to_trajectories_) {
    auto const& trajectory = pair.second;
    if (trajectory->empty()) {
      return true;
    }
  }
  return false;
}

template<typename Frame>
Instant Ephemeris<Frame>::t_min() const {
  Instant t_min;
  for (auto const& pair : bodies_to_trajectories_) {
    auto const& trajectory = pair.second;
    t_min = std::max(t_min, trajectory->t_min());
  }
  return t_min;
}

template<typename Frame>
Instant Ephemeris<Frame>::t_max() const {
  Instant t_max = bodies_to_trajectories_.begin()->second->t_max();
  for (auto const& pair : bodies_to_trajectories_) {
    auto const& trajectory = pair.second;
    t_max = std::min(t_max, trajectory->t_max());
  }
  return t_max;
}

template<typename Frame>
FixedStepSizeIntegrator<
    typename Ephemeris<Frame>::NewtonianMotionEquation> const&
Ephemeris<Frame>::planetary_integrator() const {
  return planetary_integrator_;
}

template<typename Frame>
void Ephemeris<Frame>::ForgetAfter(Instant const & t) {
  auto it = std::lower_bound(
                intermediate_states_.begin(), intermediate_states_.end(), t,
                [](typename NewtonianMotionEquation::SystemState const& left,
                   Instant const& right) {
                  return left.time.value < right;
                });
  if (it == intermediate_states_.end()) {
    return;
  }
  CHECK_LE(t, it->time.value);

  int index = 0;
  for (auto const& trajectory : trajectories_) {
    trajectory->ForgetAfter(
        it->time.value,
        DegreesOfFreedom<Frame>(it->positions[index].value,
                                it->velocities[index].value));
    ++index;
  }
  last_state_ = *it;
  intermediate_states_.erase(it, intermediate_states_.end());
}

template<typename Frame>
void Ephemeris<Frame>::ForgetBefore(Instant const& t) {
  for (auto& pair : bodies_to_trajectories_) {
    ContinuousTrajectory<Frame>& trajectory = *pair.second;
    trajectory.ForgetBefore(t);
  }
}

template<typename Frame>
void Ephemeris<Frame>::Prolong(Instant const& t) {
  IntegrationProblem<NewtonianMotionEquation> problem;
  problem.equation = massive_bodies_equation_;
  problem.append_state =
      std::bind(&Ephemeris::AppendMassiveBodiesState, this, _1);

  // Note that |t| may be before the last time that we integrated and still
  // after |t_max()|.  In this case we want to make sure that the integrator
  // makes progress.
  problem.initial_state = &last_state_;
  if (t <= last_state_.time.value) {
    problem.t_final = last_state_.time.value + step_;
  } else {
    problem.t_final = t;
  }

  // Perform the integration.  Note that we may have to iterate until |t_max()|
  // actually reaches |t| because the last series may not be fully determined
  // after the first integration.
  while (t_max() < t) {
    planetary_integrator_.Solve(problem, step_);
    // Here |problem.initial_state| still points at |last_state_|, which is the
    // state at the end of the previous call to |Solve|.  It is therefore the
    // right initial state for the next call to |Solve|, if any.
    problem.t_final += step_;
  }
}

template<typename Frame>
void Ephemeris<Frame>::FlowWithAdaptiveStep(
    not_null<DiscreteTrajectory<Frame>*> const trajectory,
    IntrinsicAcceleration intrinsic_acceleration,
    Length const& length_integration_tolerance,
    Speed const& speed_integration_tolerance,
    AdaptiveStepSizeIntegrator<NewtonianMotionEquation> const& integrator,
    Instant const& t) {
  std::vector<not_null<DiscreteTrajectory<Frame>*>> const trajectories =
      {trajectory};
  std::vector<IntrinsicAcceleration> const intrinsic_accelerations =
      {std::move(intrinsic_acceleration)};
  if (empty() || t > t_max()) {
    Prolong(t);
  }

  std::vector<typename ContinuousTrajectory<Frame>::Hint> hints(bodies_.size());
  NewtonianMotionEquation massless_body_equation;
  massless_body_equation.compute_acceleration =
      std::bind(&Ephemeris::ComputeMasslessBodiesTotalAccelerations,
                this,
                std::cref(trajectories), std::cref(intrinsic_accelerations),
                _1, _2, _3, &hints);

  typename NewtonianMotionEquation::SystemState initial_state;
  auto const trajectory_last = trajectory->last();
  auto const last_degrees_of_freedom = trajectory_last.degrees_of_freedom();
  initial_state.time = trajectory_last.time();
  initial_state.positions.push_back(last_degrees_of_freedom.position());
  initial_state.velocities.push_back(last_degrees_of_freedom.velocity());

  IntegrationProblem<NewtonianMotionEquation> problem;
  problem.equation = massless_body_equation;
  problem.append_state =
      std::bind(&Ephemeris::AppendMasslessBodiesState,
                _1, std::cref(trajectories));
  problem.t_final = t;
  problem.initial_state = &initial_state;

  AdaptiveStepSize<NewtonianMotionEquation> step_size;
  step_size.first_time_step = problem.t_final - initial_state.time.value;
  step_size.safety_factor = 0.9;
  step_size.tolerance_to_error_ratio =
      std::bind(&Ephemeris<Frame>::ToleranceToErrorRatio,
                std::cref(length_integration_tolerance),
                std::cref(speed_integration_tolerance),
                _1, _2);

  integrator.Solve(problem, step_size);
}

template<typename Frame>
void Ephemeris<Frame>::FlowWithFixedStep(
    std::vector<not_null<DiscreteTrajectory<Frame>*>> const& trajectories,
    std::vector<IntrinsicAcceleration> const& intrinsic_accelerations,
    Time const& step,
    Instant const& t) {
  VLOG(1) << __FUNCTION__ << " " << NAMED(step) << " " << NAMED(t);
  if (empty() || t > t_max()) {
    Prolong(t);
  }

  std::vector<typename ContinuousTrajectory<Frame>::Hint> hints(bodies_.size());
  NewtonianMotionEquation massless_body_equation;
  massless_body_equation.compute_acceleration =
      std::bind(&Ephemeris::ComputeMasslessBodiesTotalAccelerations,
                this,
                std::cref(trajectories), std::cref(intrinsic_accelerations),
                _1, _2, _3, &hints);

  typename NewtonianMotionEquation::SystemState initial_state;
  for (auto const& trajectory : trajectories) {
    auto const trajectory_last = trajectory->last();
    auto const last_degrees_of_freedom = trajectory_last.degrees_of_freedom();
    initial_state.time = trajectory_last.time();
    initial_state.positions.push_back(last_degrees_of_freedom.position());
    initial_state.velocities.push_back(last_degrees_of_freedom.velocity());
  }

  IntegrationProblem<NewtonianMotionEquation> problem;
  problem.equation = massless_body_equation;

#if defined(WE_LOVE_228)
  typename NewtonianMotionEquation::SystemState last_state;
  problem.append_state =
      [&last_state](
          typename NewtonianMotionEquation::SystemState const& state) {
        last_state = state;
      };
#else
  problem.append_state =
      std::bind(&Ephemeris::AppendMasslessBodiesState,
                _1, std::cref(trajectories));
#endif
  problem.t_final = t;
  problem.initial_state = &initial_state;

  planetary_integrator_.Solve(problem, step);

#if defined(WE_LOVE_228)
  AppendMasslessBodiesState(last_state, trajectories);
#endif
}

template<typename Frame>
Vector<Acceleration, Frame> Ephemeris<Frame>::ComputeGravitationalAcceleration(
    not_null<DiscreteTrajectory<Frame>*> const trajectory,
    Instant const & t) {
  auto const it = trajectory->Find(t);
  DegreesOfFreedom<Frame> const& degrees_of_freedom = it.current()->second;

  // To avoid intrinsic accelerations.
  DiscreteTrajectory<Frame> empty_trajectory;

  std::vector<Vector<Acceleration, Frame>> accelerations(1);
  std::vector<typename ContinuousTrajectory<Frame>::Hint> hints(bodies_.size());
  ComputeMasslessBodiesGravitationalAccelerations(
      {&empty_trajectory},
      t,
      {degrees_of_freedom.position()},
      &accelerations,
      &hints);

  return accelerations[0];
}

template<typename Frame>
Vector<Acceleration, Frame> Ephemeris<Frame>::ComputeGravitationalAcceleration(
    not_null<MassiveBody const*> const body,
    Instant const& t) {
  bool const body_is_oblate = body->is_oblate();

  // |other_xxx_bodies| is |xxx_bodies_| without |body|.  Index 0 in |positions|
  // and |accelerations| corresponds to |body|, the other indices to
  // |other_xxx_bodies|.
  std::vector<not_null<MassiveBody const*>> other_oblate_bodies;
  std::vector<not_null<MassiveBody const*>> other_spherical_bodies;
  std::vector<Position<Frame>> positions;
  std::vector<Vector<Acceleration, Frame>> accelerations(bodies_.size());

  // Make room for |body|.
  positions.resize(1);

  // Fill |other_xxx_bodies| and evaluate the |positions|.
  std::vector<typename ContinuousTrajectory<Frame>::Hint> hints(bodies_.size());
  for (int b = 0; b < bodies_.size(); ++b) {
    auto const& other_body = bodies_[b];
    auto const& other_body_trajectory = trajectories_[b];
    if (other_body.get() == body) {
      positions[0] = other_body_trajectory->EvaluatePosition(t, &hints[b]);
    } else if (b < number_of_oblate_bodies_) {
      CHECK(other_body->is_oblate());
      other_oblate_bodies.push_back(other_body.get());
      positions.push_back(
          other_body_trajectory->EvaluatePosition(t, &hints[b]));
    } else {
      CHECK(!other_body->is_oblate());
      other_spherical_bodies.push_back(other_body.get());
      positions.push_back(
          other_body_trajectory->EvaluatePosition(t, &hints[b]));
    }
  }

  if (body_is_oblate) {
    ComputeGravitationalAccelerationByMassiveBodyOnMassiveBodies<
        true /*body1_is_oblate*/,
        true /*body2_is_oblate*/>(
        *body /*body1*/, 0 /*b1*/,
        other_oblate_bodies /*bodies2*/,
        1 /*b2_begin*/, other_oblate_bodies.size() + 1 /*b2_end*/,
        positions, &accelerations);
    ComputeGravitationalAccelerationByMassiveBodyOnMassiveBodies<
        true /*body1_is_oblate*/,
        false /*body2_is_oblate*/>(
        *body /*body1*/, 0 /*b1*/,
        other_spherical_bodies /*bodies2*/,
        other_oblate_bodies.size() + 1/*b2_begin*/, bodies_.size() /*b2_end*/,
        positions, &accelerations);
  } else {
    ComputeGravitationalAccelerationByMassiveBodyOnMassiveBodies<
        false /*body1_is_oblate*/,
        true /*body2_is_oblate*/>(
        *body /*body1*/, 0 /*b1*/,
        other_oblate_bodies /*bodies2*/,
        1 /*b2_begin*/, other_oblate_bodies.size() + 1 /*b2_end*/,
        positions, &accelerations);
    ComputeGravitationalAccelerationByMassiveBodyOnMassiveBodies<
        false /*body1_is_oblate*/,
        false /*body2_is_oblate*/>(
        *body /*body1*/, 0 /*b1*/,
        other_spherical_bodies /*bodies2*/,
        other_oblate_bodies.size() + 1 /*b2_begin*/, bodies_.size() /*b2_end*/,
        positions, &accelerations);
  }

  return accelerations[0];
}

template<typename Frame>
void Ephemeris<Frame>::WriteToMessage(
    not_null<serialization::Ephemeris*> const message) const {
  LOG(INFO) << __FUNCTION__;
  // The bodies are serialized in the order in which they were given at
  // construction.
  for (auto const& unowned_body : unowned_bodies_) {
    unowned_body->WriteToMessage(message->add_body());
  }
  // The trajectories are serialized in the order resulting from the separation
  // between oblate and spherical bodies.
  for (auto const& trajectory : trajectories_) {
    trajectory->WriteToMessage(message->add_trajectory());
  }
  planetary_integrator_.WriteToMessage(message->mutable_planetary_integrator());
  step_.WriteToMessage(message->mutable_step());
  fitting_tolerance_.WriteToMessage(message->mutable_fitting_tolerance());
  last_state_.WriteToMessage(message->mutable_last_state());
  LOG(INFO) << NAMED(message->SpaceUsed());
  LOG(INFO) << NAMED(message->ByteSize());
}

template<typename Frame>
std::unique_ptr<Ephemeris<Frame>> Ephemeris<Frame>::ReadFromMessage(
    serialization::Ephemeris const& message) {
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  for (auto const& body : message.body()) {
    bodies.push_back(MassiveBody::ReadFromMessage(body));
  }
  auto const& planetary_integrator =
      FixedStepSizeIntegrator<NewtonianMotionEquation>::ReadFromMessage(
          message.planetary_integrator());
  auto const step = Time::ReadFromMessage(message.step());
  auto const fitting_tolerance =
      Length::ReadFromMessage(message.fitting_tolerance());

  // Dummy initial state and time.  We'll overwrite them later.
  std::vector<DegreesOfFreedom<Frame>> const initial_state(
      bodies.size(),
      DegreesOfFreedom<Frame>(Position<Frame>(), Velocity<Frame>()));
  Instant const initial_time;
  std::unique_ptr<Ephemeris<Frame>> ephemeris =
      std::make_unique<Ephemeris<Frame>>(std::move(bodies),
                                         initial_state,
                                         initial_time,
                                         planetary_integrator,
                                         step,
                                         fitting_tolerance);
  ephemeris->last_state_ =
      NewtonianMotionEquation::SystemState::ReadFromMessage(
          message.last_state());
  int index = 0;
  ephemeris->bodies_to_trajectories_.clear();
  ephemeris->trajectories_.clear();
  for (auto const& trajectory : message.trajectory()) {
    not_null<MassiveBody const*> const body = ephemeris->bodies_[index].get();
    not_null<std::unique_ptr<ContinuousTrajectory<Frame>>>
        deserialized_trajectory =
            ContinuousTrajectory<Frame>::ReadFromMessage(trajectory);
    ephemeris->trajectories_.push_back(deserialized_trajectory.get());
    ephemeris->bodies_to_trajectories_.emplace(
        body, std::move(deserialized_trajectory));
    ++index;
  }
  return ephemeris;
}

template<typename Frame>
std::unique_ptr<Ephemeris<Frame>> Ephemeris<Frame>::ReadFromPreBourbakiMessages(
    google::protobuf::RepeatedPtrField<
        serialization::Plugin::CelestialAndProperties> const& messages,
    FixedStepSizeIntegrator<NewtonianMotionEquation> const&
        planetary_integrator,
    Time const& step,
    Length const& fitting_tolerance) {
  LOG(INFO) << "Reading "<< messages.SpaceUsedExcludingSelf()
            << " bytes in pre-Bourbaki compatibility mode ";
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies;
  std::vector<DegreesOfFreedom<Frame>> initial_state;
  std::vector<std::unique_ptr<DiscreteTrajectory<Frame>>> histories;
  std::set<Instant> initial_time;
  std::set<Instant> final_time;
  for (auto const& message : messages) {
    serialization::Celestial const& celestial = message.celestial();
    bodies.emplace_back(MassiveBody::ReadFromMessage(celestial.body()));
    histories.emplace_back(DiscreteTrajectory<Frame>::ReadFromMessage(
        celestial.history_and_prolongation().history()));
    auto const prolongation =
        DiscreteTrajectory<Frame>::ReadPointerFromMessage(
            celestial.history_and_prolongation().prolongation(),
            histories.back().get());
    initial_state.push_back(histories.back()->first().degrees_of_freedom());
    initial_time.insert(histories.back()->first().time());
    final_time.insert(prolongation->last().time());
  }
  CHECK_EQ(1, initial_time.size());
  CHECK_EQ(1, final_time.size());
  LOG(INFO) << "Initial time is " << *initial_time.cbegin()
            << ", final time is " << *final_time.cbegin();

  // Construct a new ephemeris using the bodies and initial states and time
  // extracted from the serialized celestials.
  auto ephemeris = std::make_unique<Ephemeris<Frame>>(std::move(bodies),
                                                      initial_state,
                                                      *initial_time.cbegin(),
                                                      planetary_integrator,
                                                      step,
                                                      fitting_tolerance);

  // Extend the continuous trajectories using the data from the discrete
  // trajectories.
  std::set<Instant> last_state_time;
  for (int i = 0; i < histories.size(); ++i) {
    auto* const body = ephemeris->unowned_bodies_[i];
    auto const& history = histories[i];

    // Apologies about the linear search here but we need the index |j| to fill
    // |last_state_| and I don't feel like introducing another data structure
    // just to speed up the compatibility case.
    int j = 0;
    for (; j < ephemeris->bodies_.size(); ++j) {
      if (body == ephemeris->bodies_[j].get()) {
        break;
      }
    }
    auto continuous_trajectory = ephemeris->trajectories_[j];

    auto it = history->first();
    Instant last_time = it.time();
    DegreesOfFreedom<Frame> last_degrees_of_freedom = it.degrees_of_freedom();
    for (; !it.at_end(); ++it) {
      Time const duration_since_last_time = it.time() - last_time;
      if (duration_since_last_time == step) {
        // A time in the discrete trajectory that is aligned on the continuous
        // trajectory.
        last_time = it.time();
        last_degrees_of_freedom = it.degrees_of_freedom();
        continuous_trajectory->Append(last_time, last_degrees_of_freedom);
      } else if (duration_since_last_time > step) {
        // A time in the discrete trajectory that is not aligned on the
        // continuous trajectory.  Stop here, we'll use prolong to recompute the
        // rest.
        break;
      }
    }

    // Fill the |last_state_| for this body.  It will be the starting state for
    // Prolong.
    last_state_time.insert(last_time);
    ephemeris->last_state_.positions[j] = last_degrees_of_freedom.position();
    ephemeris->last_state_.velocities[j] = last_degrees_of_freedom.velocity();
  }
  CHECK_EQ(1, last_state_time.size());
  ephemeris->last_state_.time = *last_state_time.cbegin();
  LOG(INFO) << "Last time in discrete trajectories is "
            << *last_state_time.cbegin();

  // Prolong the ephemeris to the final time.  This might create discrepancies
  // from the discrete trajectories.
  ephemeris->Prolong(*final_time.cbegin());

  return ephemeris;
}

template <typename Frame>
Ephemeris<Frame>::Ephemeris()
    : planetary_integrator_(DummyIntegrator<Frame>::Instance()) {}

template<typename Frame>
void Ephemeris<Frame>::AppendMassiveBodiesState(
    typename NewtonianMotionEquation::SystemState const& state) {
  last_state_ = state;
  int index = 0;
  for (auto& trajectory : trajectories_) {
    trajectory->Append(
        state.time.value,
        DegreesOfFreedom<Frame>(state.positions[index].value,
                                state.velocities[index].value));
    ++index;
  }

  // Record an intermediate state if we haven't done so for too long and this
  // time is a |t_max|.
  CHECK(!trajectories_.empty());
  Instant const t_max = trajectories_.front()->t_max();
  if (t_max == state.time.value) {
    Instant const t_last_intermediate_state =
        intermediate_states_.empty()
            ? Instant(-std::numeric_limits<double>::infinity() * Second)
            : intermediate_states_.back().time.value;
    CHECK_LE(t_last_intermediate_state, t_max);
    if (t_max - t_last_intermediate_state > kMaxTimeBetweenIntermediateStates) {
      intermediate_states_.push_back(state);
    }
  }
}

template<typename Frame>
void Ephemeris<Frame>::AppendMasslessBodiesState(
    typename NewtonianMotionEquation::SystemState const& state,
    std::vector<not_null<DiscreteTrajectory<Frame>*>> const& trajectories) {
  int index = 0;
  for (auto& trajectory : trajectories) {
    trajectory->Append(
        state.time.value,
        DegreesOfFreedom<Frame>(state.positions[index].value,
                                state.velocities[index].value));
    ++index;
  }
}

template<typename Frame>
template<bool body1_is_oblate,
         bool body2_is_oblate>
void Ephemeris<Frame>::
ComputeGravitationalAccelerationByMassiveBodyOnMassiveBodies(
    MassiveBody const& body1,
    size_t const b1,
    std::vector<not_null<MassiveBody const*>> const& bodies2,
    size_t const b2_begin,
    size_t const b2_end,
    std::vector<Position<Frame>> const& positions,
    not_null<std::vector<Vector<Acceleration, Frame>>*> const accelerations) {
  Vector<Acceleration, Frame>& acceleration_on_b1 = (*accelerations)[b1];
  GravitationalParameter const& μ1 = body1.gravitational_parameter();
  for (std::size_t b2 = std::max(b1 + 1, b2_begin); b2 < b2_end; ++b2) {
    Vector<Acceleration, Frame>& acceleration_on_b2 = (*accelerations)[b2];
    MassiveBody const& body2 = *bodies2[b2 - b2_begin];
    GravitationalParameter const& μ2 = body2.gravitational_parameter();

    Displacement<Frame> const Δq = positions[b1] - positions[b2];

    Exponentiation<Length, 2> const Δq_squared = InnerProduct(Δq, Δq);
    // NOTE(phl): Don't try to compute one_over_Δq_squared here, it makes the
    // non-oblate path slower.
    Exponentiation<Length, -3> const one_over_Δq_cubed =
        Sqrt(Δq_squared) / (Δq_squared * Δq_squared);

    auto const μ1_over_Δq_cubed = μ1 * one_over_Δq_cubed;
    acceleration_on_b2 += Δq * μ1_over_Δq_cubed;

    // Lex. III. Actioni contrariam semper & æqualem esse reactionem:
    // sive corporum duorum actiones in se mutuo semper esse æquales &
    // in partes contrarias dirigi.
    auto const μ2_over_Δq_cubed = μ2 * one_over_Δq_cubed;
    acceleration_on_b1 -= Δq * μ2_over_Δq_cubed;

    if (body1_is_oblate || body2_is_oblate) {
      Exponentiation<Length, -2> const one_over_Δq_squared = 1 / Δq_squared;
      if (body1_is_oblate) {
        Vector<Quotient<Acceleration,
                        GravitationalParameter>, Frame> const
            order_2_zonal_effect1 =
                Order2ZonalEffect<Frame>(
                    static_cast<OblateBody<Frame> const&>(body1),
                    Δq,
                    one_over_Δq_squared,
                    one_over_Δq_cubed);
        acceleration_on_b1 -= μ2 * order_2_zonal_effect1;
        acceleration_on_b2 += μ1 * order_2_zonal_effect1;
      }
      if (body2_is_oblate) {
        Vector<Quotient<Acceleration,
                        GravitationalParameter>, Frame> const
            order_2_zonal_effect2 =
                Order2ZonalEffect<Frame>(
                    static_cast<OblateBody<Frame> const&>(body2),
                    Δq,
                    one_over_Δq_squared,
                    one_over_Δq_cubed);
        acceleration_on_b1 -= μ2 * order_2_zonal_effect2;
        acceleration_on_b2 += μ1 * order_2_zonal_effect2;
      }
    }
  }
}

template<typename Frame>
template<bool body1_is_oblate>
void Ephemeris<Frame>::
ComputeGravitationalAccelerationByMassiveBodyOnMasslessBodies(
    Instant const& t,
    MassiveBody const& body1,
    size_t const b1,
    std::vector<Position<Frame>> const& positions,
    not_null<std::vector<Vector<Acceleration, Frame>>*> const accelerations,
    not_null<std::vector<typename ContinuousTrajectory<Frame>::Hint>*>
        const hints) {
  GravitationalParameter const& μ1 = body1.gravitational_parameter();

  for (size_t b2 = 0; b2 < positions.size(); ++b2) {
    Displacement<Frame> const Δq =
        trajectories_[b1]->EvaluatePosition(t, &(*hints)[b1]) - positions[b2];

    Exponentiation<Length, 2> const Δq_squared = InnerProduct(Δq, Δq);
    // NOTE(phl): Don't try to compute one_over_Δq_squared here, it makes the
    // non-oblate path slower.
    Exponentiation<Length, -3> const one_over_Δq_cubed =
        Sqrt(Δq_squared) / (Δq_squared * Δq_squared);

    auto const μ1_over_Δq_cubed = μ1 * one_over_Δq_cubed;
    (*accelerations)[b2] += Δq * μ1_over_Δq_cubed;

    if (body1_is_oblate) {
      Exponentiation<Length, -2> const one_over_Δq_squared = 1 / Δq_squared;
      Vector<Quotient<Acceleration,
                      GravitationalParameter>, Frame> const
          order_2_zonal_effect1 =
              Order2ZonalEffect<Frame>(
                  static_cast<OblateBody<Frame> const &>(body1),
                  Δq,
                  one_over_Δq_squared,
                  one_over_Δq_cubed);
      (*accelerations)[b2] += μ1 * order_2_zonal_effect1;
    }
  }
}

template<typename Frame>
void Ephemeris<Frame>::ComputeMassiveBodiesGravitationalAccelerations(
    Instant const& t,
    std::vector<Position<Frame>> const& positions,
    not_null<std::vector<Vector<Acceleration, Frame>>*> const accelerations) {
  accelerations->assign(accelerations->size(), Vector<Acceleration, Frame>());

  for (std::size_t b1 = 0; b1 < number_of_oblate_bodies_; ++b1) {
    MassiveBody const& body1 = *oblate_bodies_[b1];
    ComputeGravitationalAccelerationByMassiveBodyOnMassiveBodies<
        true /*body1_is_oblate*/,
        true /*body2_is_oblate*/>(
        body1, b1,
        oblate_bodies_ /*bodies2*/,
        0 /*b2_begin*/,
        number_of_oblate_bodies_ /*b2_end*/,
        positions,
        accelerations);
    ComputeGravitationalAccelerationByMassiveBodyOnMassiveBodies<
        true /*body1_is_oblate*/,
        false /*body2_is_oblate*/>(
        body1, b1,
        spherical_bodies_ /*bodies2*/,
        number_of_oblate_bodies_ /*b2_begin*/,
        number_of_oblate_bodies_ +
            number_of_spherical_bodies_ /*b2_end*/,
        positions,
        accelerations);
  }
  for (std::size_t b1 = number_of_oblate_bodies_;
       b1 < number_of_oblate_bodies_ +
            number_of_spherical_bodies_;
       ++b1) {
    MassiveBody const& body1 =
        *spherical_bodies_[b1 - number_of_oblate_bodies_];
    ComputeGravitationalAccelerationByMassiveBodyOnMassiveBodies<
        false /*body1_is_oblate*/,
        false /*body2_is_oblate*/>(
        body1, b1,
        spherical_bodies_ /*bodies2*/,
        number_of_oblate_bodies_ /*b2_begin*/,
        number_of_oblate_bodies_ +
            number_of_spherical_bodies_ /*b2_end*/,
        positions,
        accelerations);
  }
}

template<typename Frame>
void Ephemeris<Frame>::ComputeMasslessBodiesGravitationalAccelerations(
      std::vector<not_null<DiscreteTrajectory<Frame>*>> const& trajectories,
      Instant const& t,
      std::vector<Position<Frame>> const& positions,
      not_null<std::vector<Vector<Acceleration, Frame>>*> const accelerations,
      not_null<std::vector<typename ContinuousTrajectory<Frame>::Hint>*>
          const hints) {
  CHECK_EQ(trajectories.size(), positions.size());
  CHECK_EQ(trajectories.size(), accelerations->size());
  accelerations->assign(accelerations->size(), Vector<Acceleration, Frame>());

  for (std::size_t b1 = 0; b1 < number_of_oblate_bodies_; ++b1) {
    MassiveBody const& body1 = *oblate_bodies_[b1];
    ComputeGravitationalAccelerationByMassiveBodyOnMasslessBodies<
        true /*body1_is_oblate*/>(
        t,
        body1, b1,
        positions,
        accelerations,
        hints);
  }
  for (std::size_t b1 = number_of_oblate_bodies_;
       b1 < number_of_oblate_bodies_ +
            number_of_spherical_bodies_;
       ++b1) {
    MassiveBody const& body1 =
        *spherical_bodies_[b1 - number_of_oblate_bodies_];
    ComputeGravitationalAccelerationByMassiveBodyOnMasslessBodies<
        false /*body1_is_oblate*/>(
        t,
        body1, b1,
        positions,
        accelerations,
        hints);
  }
}

template<typename Frame>
void Ephemeris<Frame>::ComputeMasslessBodiesTotalAccelerations(
    std::vector<not_null<DiscreteTrajectory<Frame>*>> const& trajectories,
    std::vector<IntrinsicAcceleration> const& intrinsic_accelerations,
    Instant const& t,
    std::vector<Position<Frame>> const& positions,
    not_null<std::vector<Vector<Acceleration, Frame>>*> const accelerations,
    not_null<std::vector<typename ContinuousTrajectory<Frame>::Hint>*>
        const hints) {
  // First, the acceleration due to the gravitational field of the
  // massive bodies.
  ComputeMasslessBodiesGravitationalAccelerations(
      trajectories, t, positions, accelerations, hints);

  // Then, the intrinsic accelerations, if any.
  if (!intrinsic_accelerations.empty()) {
    CHECK_EQ(trajectories.size(), intrinsic_accelerations.size());
    for (int i = 0; i < intrinsic_accelerations.size(); ++i) {
      auto const intrinsic_acceleration = intrinsic_accelerations[i];
      if (intrinsic_acceleration != nullptr) {
        (*accelerations)[i] += intrinsic_acceleration(t);
      }
    }
  }
}

template<typename Frame>
double Ephemeris<Frame>::ToleranceToErrorRatio(
    Length const& length_integration_tolerance,
    Speed const& speed_integration_tolerance,
    Time const& current_step_size,
    typename NewtonianMotionEquation::SystemStateError const& error) {
  Length max_length_error;
  Speed max_speed_error;
  for (auto const& position_error : error.position_error) {
    max_length_error = std::max(max_length_error,
                                position_error.Norm());
  }
  for (auto const& velocity_error : error.velocity_error) {
    max_speed_error = std::max(max_speed_error,
                               velocity_error.Norm());
  }
  return std::min(length_integration_tolerance / max_length_error,
                  speed_integration_tolerance / max_speed_error);
}

}  // namespace physics
}  // namespace principia
