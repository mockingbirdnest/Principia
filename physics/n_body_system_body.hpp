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
    std::unique_ptr<Bodies>&& massless_bodies)
    : massless_bodies_(std::move(massless_bodies)),
      massive_bodies_(std::move(massive_bodies)) {
  // Parameter checking.
  if (massive_bodies_ != nullptr) {
    for (auto const& body : *massive_bodies_) {
#ifndef _MANAGED
      CHECK(!body->is_massless());
#endif
      bodies_.push_back(body.get());
    }
  }
  if (massless_bodies_ != nullptr) {
    for (auto const& body : *massless_bodies_) {
#ifndef _MANAGED
      CHECK(body->is_massless());
#endif
      bodies_.push_back(body.get());
    }
  }
}

template<typename InertialFrame>
std::vector<Body<InertialFrame> const*>
NBodySystem<InertialFrame>::massless_bodies() const {
  std::vector<Body<InertialFrame> const*> result;
  for (auto const& body : *massless_bodies_) {
    result.push_back(body.get());
  }
  return result;
}

template<typename InertialFrame>
std::vector<Body<InertialFrame> const*>
NBodySystem<InertialFrame>::massive_bodies() const {
  std::vector<Body<InertialFrame> const*> result;
  for (auto const& body : *massive_bodies_) {
    result.push_back(body.get());
  }
  return result;
}

template<typename InertialFrame>
std::vector<Body<InertialFrame> const*>
NBodySystem<InertialFrame>::bodies() const {
  std::vector<Body<InertialFrame> const*> result;
  for (auto const& body : bodies_) {
    result.push_back(body);
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
  SymplecticIntegrator<Length, Speed>::Solution solution;

  // Prepare the input data.
  std::unique_ptr<Time> reference_time;
  for (Body<InertialFrame> const* body : bodies_) {
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
    Body<InertialFrame>* body = bodies_[i / 3];
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
  for (size_t b1 = 0, three_b1 = 0; b1 < bodies_.size(); ++b1, three_b1 += 3) {
    Body<InertialFrame> const& body1 = *bodies_[b1];
    bool const body1_is_massless = body1.is_massless();
    for (size_t b2 = b1 + 1; b2 < bodies_.size(); ++b2) {
      Body<InertialFrame> const& body2 = *bodies_[b2];
      bool const body2_is_massless = body2.is_massless();
      size_t const three_b2 = 3 * b2;
      if (!body1_is_massless || !body2_is_massless) {
        Length const Δq0 = q[three_b1] - q[three_b2];
        Length const Δq1 = q[three_b1 + 1] - q[three_b2 + 1];
        Length const Δq2 = q[three_b1 + 2] - q[three_b2 + 2];

        Exponentiation<Length, 2> const squared_distance =
            Δq0 * Δq0 + Δq1 * Δq1 + Δq2 * Δq2;
        Exponentiation<Length, -3> const multiplier =
            Sqrt(squared_distance) / (squared_distance * squared_distance);

        if (!body2_is_massless) {
          auto const μ2OverRSquared =
              body2.gravitational_parameter() * multiplier;
          (*result)[three_b1] -= Δq0 * μ2OverRSquared;
          (*result)[three_b1 + 1] -= Δq1 * μ2OverRSquared;
          (*result)[three_b1 + 2] -= Δq2 * μ2OverRSquared;
        }
        // Lex. III. Actioni contrariam semper & æqualem esse reactionem:
        // sive corporum duorum actiones in se mutuo semper esse æquales &
        // in partes contrarias dirigi.
        if (!body1_is_massless) {
          auto const μ1OverRSquared =
              body1.gravitational_parameter() * multiplier;
          (*result)[three_b2] += Δq0 * μ1OverRSquared;
          (*result)[three_b2 + 1] += Δq1 * μ1OverRSquared;
          (*result)[three_b2 + 2] += Δq2 * μ1OverRSquared;
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
