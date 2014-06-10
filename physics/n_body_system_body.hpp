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
using principia::quantities::SIUnit;

namespace principia {
namespace physics {

template<typename InertialFrame>
NBodySystem<InertialFrame>::NBodySystem(
    std::vector<Body<InertialFrame>*> const* bodies)
    : bodies_(bodies) {}

template<typename InertialFrame>
NBodySystem<InertialFrame>::~NBodySystem() {
  for (Body<InertialFrame>* body : *bodies_) {
    delete body;
  }
}

template<typename InertialFrame>
std::vector<Body<InertialFrame>*> const& NBodySystem<InertialFrame>::bodies() {
  return *bodies_;
}

template<typename InertialFrame>
void NBodySystem<InertialFrame>::Integrate(
    SymplecticIntegrator<Length, Speed> const& integrator,
    Time const& tmax,
    Time const& Δt,
    int const sampling_period) {
  SymplecticIntegrator<Length, Speed>::Parameters parameters;
  SymplecticIntegrator<Length, Speed>::Solution solution;

  // Prepare the input data.
  std::unique_ptr<Time> reference_time;
  for (Body<InertialFrame> const* body : *bodies_) {
    R3Element<Length> const& position = body->positions().back().coordinates();
    R3Element<Speed> const& velocity = body->velocities().back().coordinates();
    Time const& time = body->times().back();
    for (int i = 0; i < 3; ++i) {
      parameters.q0.push_back(position[i]);
    }
    for (int i = 0; i < 3; ++i) {
      parameters.p0.push_back(velocity[i]);
    }
    parameters.t0 = time;
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

#ifndef _MANAGED
  // TODO(phl): Ignoring errors for now.
  CHECK_EQ(solution.position.size(), solution.momentum.size());
#endif
  std::vector<Time> const& t = solution.time.quantities;
  // Loop over the bodies.
  // TODO(phl): It looks like we are transposing in the integrator and then
  // transposing here again.
  for (size_t i = 0; i < solution.position.size(); i += 3) {
    Body<InertialFrame>* body = (*bodies_)[i / 3];
    std::vector<Length> const& q0 = solution.position[i + 0].quantities;
    std::vector<Length> const& q1 = solution.position[i + 1].quantities;
    std::vector<Length> const& q2 = solution.position[i + 2].quantities;
    std::vector<Speed> const& p0 = solution.momentum[i + 0].quantities;
    std::vector<Speed> const& p1 = solution.momentum[i + 1].quantities;
    std::vector<Speed> const& p2 = solution.momentum[i + 2].quantities;
#ifndef _MANAGED
    CHECK_EQ(t.size(), q0.size());
    CHECK_EQ(t.size(), q1.size());
    CHECK_EQ(t.size(), q2.size());
    CHECK_EQ(t.size(), p0.size());
    CHECK_EQ(t.size(), p1.size());
    CHECK_EQ(t.size(), p2.size());
#endif
    for (size_t j = 0; j < t.size(); ++j) {
      Vector<Length, InertialFrame> const position(
          R3Element<Length>(q0[j], q1[j], q2[j]));
      Vector<Speed, InertialFrame> const velocity(
          R3Element<Speed>(p0[j], p1[j], p2[j]));
      Time const time = t[j];
      body->AppendToTrajectory(position, velocity, time);
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
  for (size_t b1 = 0; b1 < bodies_->size(); ++b1) {
    for (size_t b2 = b1 + 1; b2 < bodies_->size(); ++b2) {
      if (!(*bodies_)[b1]->is_massless() || !(*bodies_)[b2]->is_massless()) {
        Length const Δq0 = q[3 * b1] - q[3 * b2];
        Length const Δq1 = q[3 * b1 + 1] - q[3 * b2 + 1];
        Length const Δq2 = q[3 * b1 + 2] - q[3 * b2 + 2];

        Exponentiation<Length, 2> const squared_distance =
            Δq0 * Δq0 + Δq1 * Δq1 + Δq2 * Δq2;
        Exponentiation<Length, 3> const denominator =
            squared_distance * Sqrt(squared_distance);

        if (!(*bodies_)[b2]->is_massless()) {
          auto const μ2OverRSquared =
              (*bodies_)[b2]->gravitational_parameter() / denominator;
          (*result)[3 * b1] -= Δq0 * μ2OverRSquared;
          (*result)[3 * b1 + 1] -= Δq1 * μ2OverRSquared;
          (*result)[3 * b1 + 2] -= Δq2 * μ2OverRSquared;
        }
        // Lex. III. Actioni contrariam semper & æqualem esse reactionem:
        // sive corporum duorum actiones in se mutuo semper esse æquales &
        // in partes contrarias dirigi.
        if (!(*bodies_)[b1]->is_massless()) {
          auto const μ1OverRSquared =
              (*bodies_)[b1]->gravitational_parameter() / denominator;
          (*result)[3 * b2] += Δq0 * μ1OverRSquared;
          (*result)[3 * b2 + 1] += Δq1 * μ1OverRSquared;
          (*result)[3 * b2 + 2] += Δq2 * μ1OverRSquared;
        }
      }
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
