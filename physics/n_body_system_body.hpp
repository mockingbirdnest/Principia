#include "physics/n_body_system.hpp"

#include <cmath>
#include <functional>
#include <vector>

#include "geometry/r3_element.hpp"
#ifndef _MANAGED
#include "glog/logging.h"
#endif
#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"
#include "quantities/quantities.hpp"

using principia::geometry::R3Element;
using principia::integrators::SPRKIntegrator;
using principia::integrators::SymplecticIntegrator;
using principia::quantities::Acceleration;
using principia::quantities::Exponentiation;
using principia::quantities::GravitationalParameter;
using principia::quantities::Length;
using principia::quantities::Speed;

namespace principia {
namespace physics {

template<typename InertialFrame>
NBodySystem<InertialFrame>::NBodySystem(
    std::unique_ptr<Bodies>&& massive_bodies,
    std::unique_ptr<Bodies>&& massless_bodies,
    std::unique_ptr<Trajectories>&& trajectories)
    : massive_bodies_(massive_bodies == nullptr ?
                          std::unique_ptr<Bodies>(new Bodies) :
                          std::move(massive_bodies)),
      massless_bodies_(massless_bodies == nullptr ?
                           std::unique_ptr<Bodies>(new Bodies) :
                           std::move(massless_bodies)),
      trajectories_(std::move(trajectories)) {
  // Parameter checking.
  for (auto const& body : *massive_bodies_) {
#ifndef _MANAGED
    CHECK(!body->is_massless());
#endif
    bodies_.push_back(body.get());
  }
  for (auto const& body : *massless_bodies_) {
#ifndef _MANAGED
    CHECK(body->is_massless());
#endif
    bodies_.push_back(body.get());
  }
#ifndef _MANAGED
  CHECK_EQ(trajectories_->size(), bodies_.size());
#endif
}

template<typename InertialFrame>
std::vector<Body const*> NBodySystem<InertialFrame>::massless_bodies() const {
  std::vector<Body const*> result;
  for (auto const& body : *massless_bodies_) {
    result.push_back(body.get());
  }
  return result;
}

template<typename InertialFrame>
std::vector<Body const*> NBodySystem<InertialFrame>::massive_bodies() const {
  std::vector<Body const*> result;
  for (auto const& body : *massive_bodies_) {
    result.push_back(body.get());
  }
  return result;
}

template<typename InertialFrame>
std::vector<Body const*> NBodySystem<InertialFrame>::bodies() const {
  std::vector<Body const*> result;
  for (auto const& body : bodies_) {
    result.push_back(body);
  }
  return result;
}

template<typename InertialFrame>
std::vector<Trajectory<InertialFrame> const*>
NBodySystem<InertialFrame>::trajectories() const {
  std::vector<Trajectory<InertialFrame> const*> result;
  for (auto const& trajectory : *trajectories_) {
    result.push_back(trajectory.get());
  }
  return result;
}

template<typename InertialFrame>
void NBodySystem<InertialFrame>::Integrate(
    SymplecticIntegrator<Length, Speed> const& integrator,
    Time const& tmax,
    Time const& Δt,
    int const sampling_period) {
  SymplecticIntegrator<Length, Speed>::Parameters parameters;
  std::vector<SymplecticIntegrator<Length, Speed>::SystemState> solution;

  // Prepare the input data.
  std::unique_ptr<Time> reference_time;
  for (auto const& trajectory : *trajectories_) {
    // TODO(phl): Relation with bodies_?
    R3Element<Length> const& position =
        trajectory->last_position().coordinates();
    R3Element<Speed> const& velocity =
        trajectory->last_velocity().coordinates();
    Time const& time = trajectory->last_time();
    for (int i = 0; i < 3; ++i) {
      parameters.initial.positions.emplace_back(position[i]);
    }
    for (int i = 0; i < 3; ++i) {
      parameters.initial.momenta.emplace_back(velocity[i]);
    }
    parameters.initial.time = time;
    // All the positions/velocities must be for the same time.
    if (reference_time == nullptr) {
      reference_time.reset(new Time(time));
    } else {
#ifndef _MANAGED
      CHECK_EQ(*reference_time, time);
#endif
    }
  }

  parameters.tmax = tmax;
  parameters.Δt = Δt;
  parameters.sampling_period = sampling_period;
  dynamic_cast<const SPRKIntegrator<Length, Speed>*>(&integrator)->Solve(
      std::bind(&NBodySystem::ComputeGravitationalAccelerations, this,
                std::placeholders::_1,
                std::placeholders::_2,
                std::placeholders::_3),
      &ComputeGravitationalVelocities,
      parameters, &solution);

  // TODO(phl): Ignoring errors for now.
  // Loop over the time steps.
  for (std::size_t i = 0; i < solution.size(); ++i) {
    SymplecticIntegrator<Length, Speed>::SystemState const& state = solution[i];
    Time const& time = state.time.value;
#ifndef _MANAGED
    CHECK_EQ(state.positions.size(), state.momenta.size());
#endif
    // Loop over the dimensions.
    for (std::size_t k = 0, b = 0; k < state.positions.size(); k += 3, ++b) {
      // TODO(phl): bodies_ vs. trajectory.
      Trajectory<InertialFrame>* trajectory = (*trajectories_)[b].get();
      Vector<Length, InertialFrame> const position(
          R3Element<Length>(state.positions[k].value,
                            state.positions[k + 1].value,
                            state.positions[k + 2].value));
      Vector<Speed, InertialFrame> const velocity(
          R3Element<Speed>(state.momenta[k].value,
                           state.momenta[k + 1].value,
                           state.momenta[k + 2].value));
      trajectory->Append(time,
                         DegreesOfFreedom<InertialFrame>(position, velocity));
    }
  }
}

template<typename InertialFrame>
void NBodySystem<InertialFrame>::ComputeGravitationalAccelerations(
    Time const& t,
    std::vector<Length> const& q,
    std::vector<Acceleration>* result) const {
  result->assign(result->size(), Acceleration());

  // TODO(phl): Used to deal with proper accelerations here.
  for (std::size_t b1 = 0, three_b1 = 0;
       b1 < massive_bodies_->size();
       ++b1, three_b1 += 3) {
    GravitationalParameter const& body1_gravitational_parameter =
        (*massive_bodies_)[b1]->gravitational_parameter();
    for (std::size_t b2 = b1 + 1; b2 < massive_bodies_->size(); ++b2) {
      std::size_t const three_b2 = 3 * b2;
      Length const Δq0 = q[three_b1] - q[three_b2];
      Length const Δq1 = q[three_b1 + 1] - q[three_b2 + 1];
      Length const Δq2 = q[three_b1 + 2] - q[three_b2 + 2];

      Exponentiation<Length, 2> const squared_distance =
          Δq0 * Δq0 + Δq1 * Δq1 + Δq2 * Δq2;
      Exponentiation<Length, -3> const multiplier =
          Sqrt(squared_distance) / (squared_distance * squared_distance);

      auto const μ2OverRSquared =
          (*massive_bodies_)[b2]->gravitational_parameter() * multiplier;
      (*result)[three_b1] -= Δq0 * μ2OverRSquared;
      (*result)[three_b1 + 1] -= Δq1 * μ2OverRSquared;
      (*result)[three_b1 + 2] -= Δq2 * μ2OverRSquared;
      // Lex. III. Actioni contrariam semper & æqualem esse reactionem:
      // sive corporum duorum actiones in se mutuo semper esse æquales &
      // in partes contrarias dirigi.
      auto const μ1OverRSquared = body1_gravitational_parameter * multiplier;
      (*result)[three_b2] += Δq0 * μ1OverRSquared;
      (*result)[three_b2 + 1] += Δq1 * μ1OverRSquared;
      (*result)[three_b2 + 2] += Δq2 * μ1OverRSquared;
    }
    for (std::size_t b2 = 0; b2 < massless_bodies_->size(); ++b2) {
      std::size_t const three_b2 = 3 * b2;
      Length const Δq0 = q[three_b1] - q[three_b2];
      Length const Δq1 = q[three_b1 + 1] - q[three_b2 + 1];
      Length const Δq2 = q[three_b1 + 2] - q[three_b2 + 2];

      Exponentiation<Length, 2> const squared_distance =
          Δq0 * Δq0 + Δq1 * Δq1 + Δq2 * Δq2;
      Exponentiation<Length, -3> const multiplier =
          Sqrt(squared_distance) / (squared_distance * squared_distance);

      auto const μ1OverRSquared = body1_gravitational_parameter * multiplier;
      (*result)[three_b2] += Δq0 * μ1OverRSquared;
      (*result)[three_b2 + 1] += Δq1 * μ1OverRSquared;
      (*result)[three_b2 + 2] += Δq2 * μ1OverRSquared;
    }
  }
}

template<typename InertialFrame>
void NBodySystem<InertialFrame>::ComputeGravitationalVelocities(
    std::vector<Speed> const& p,
    std::vector<Speed>* result) {
  *result = p;
}

}  // namespace physics
}  // namespace principia
