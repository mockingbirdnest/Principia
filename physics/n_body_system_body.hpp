#include "physics/n_body_system.hpp"

#include <cmath>
#include <functional>
#include <set>
#include <vector>

#include "geometry/named_quantities.hpp"
#include "geometry/r3_element.hpp"
#include "glog/logging.h"
#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"
#include "quantities/quantities.hpp"

using principia::geometry::Instant;
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
void NBodySystem<InertialFrame>::Integrate(
    SymplecticIntegrator<Length, Speed> const& integrator,
    Instant const& tmax,
    Time const& Δt,
    int const sampling_period,
    bool const tmax_is_exact,
    Trajectories const& trajectories) {
  SymplecticIntegrator<Length, Speed>::Parameters parameters;
  std::vector<SymplecticIntegrator<Length, Speed>::SystemState> solution;

  // TODO(phl): Use a position based on the first mantissa bits of the
  // center-of-mass referential and a time in the middle of the integration
  // interval.  In the integrator itself, all quantities are "vectors" relative
  // to these references.
  Position<InertialFrame> const reference_position;
  Instant const reference_time;

  // These objects are for checking the consistency of the parameters.
  std::set<Instant> times_in_trajectories;
  std::set<Body const*> bodies_in_trajectories;

  // Prepare the initial state of the integrator.  For efficiently computing the
  // accelerations, we need to separate the trajectories of massive bodies from
  // those of massless bodies.  They are put in this order in
  // |reordered_trajectories|.  The trajectories of massive and massless bodies
  // are put in |massive_trajectories| and |massless_trajectories|,
  // respectively.
  std::vector<Trajectory<InertialFrame>*> reordered_trajectories;
  std::vector<Trajectory<InertialFrame> const*> massive_trajectories;
  std::vector<Trajectory<InertialFrame> const*> massless_trajectories;
  // This loop ensures that the massive bodies precede the massless bodies in
  // the vectors representing the initial data.
  for (bool is_massless : {false, true}) {
    for (auto const& trajectory : trajectories) {
      // See if this trajectory should be processed in this iteration and
      // update the appropriate vector.
      Body const* const body = &trajectory->body();
      if (body->is_massless() != is_massless) {
        continue;
      }
      if (is_massless) {
        massless_trajectories.push_back(trajectory);
      } else {
        massive_trajectories.push_back(trajectory);
      }
      reordered_trajectories.push_back(trajectory);

      // Fill the initial position/velocity/time.
      R3Element<Length> const& position =
          (trajectory->last_position() - reference_position).coordinates();
      R3Element<Speed> const& velocity =
          trajectory->last_velocity().coordinates();
      Instant const& time = trajectory->last_time();
      for (int i = 0; i < 3; ++i) {
        parameters.initial.positions.emplace_back(position[i]);
      }
      for (int i = 0; i < 3; ++i) {
        parameters.initial.momenta.emplace_back(velocity[i]);
      }

      // Check that all trajectories are for different bodies.
      auto const inserted = bodies_in_trajectories.insert(body);
      CHECK(inserted.second) << "Multiple trajectories for the same body";
      // The final points of all trajectories must all be for the same time.
      times_in_trajectories.insert(time);
      CHECK_GE(1U, times_in_trajectories.size())
          << "Inconsistent last time in trajectories";
    }
  }
  {
    // Beyond this point we must not use the |trajectories| parameter as it is
    // in the wrong order with respect to the data passed to the integrator.  We
    // use this block to hide it.
    Trajectories const& trajectories = reordered_trajectories;

    parameters.initial.time = *times_in_trajectories.cbegin() - reference_time;
    parameters.tmax = tmax - reference_time;
    parameters.Δt = Δt;
    parameters.sampling_period = sampling_period;
    parameters.tmax_is_exact = tmax_is_exact;
    dynamic_cast<const SPRKIntegrator<Length, Speed>*>(&integrator)->Solve(
        std::bind(&NBodySystem::ComputeGravitationalAccelerations,
                  massive_trajectories,
                  massless_trajectories,
                  reference_time,
                  std::placeholders::_1,
                  std::placeholders::_2,
                  std::placeholders::_3),
        &ComputeGravitationalVelocities,
        parameters, &solution);

    // TODO(phl): Ignoring errors for now.
    // Loop over the time steps.
    for (std::size_t i = 0; i < solution.size(); ++i) {
      SymplecticIntegrator<Length, Speed>::SystemState const& state =
          solution[i];
      Instant const time = state.time.value + reference_time;
      CHECK_EQ(state.positions.size(), state.momenta.size());
      // Loop over the dimensions.
      for (std::size_t k = 0, t = 0; k < state.positions.size(); k += 3, ++t) {
        Vector<Length, InertialFrame> const position(
            R3Element<Length>(state.positions[k].value,
                              state.positions[k + 1].value,
                              state.positions[k + 2].value));
        Velocity<InertialFrame> const velocity(
            R3Element<Speed>(state.momenta[k].value,
                             state.momenta[k + 1].value,
                             state.momenta[k + 2].value));
        trajectories[t]->Append(
            time,
            DegreesOfFreedom<InertialFrame>(position + reference_position,
                                            velocity));
      }
    }
  }
}

template<typename InertialFrame>
void NBodySystem<InertialFrame>::ComputeGravitationalAccelerations(
    std::vector<Trajectory<InertialFrame> const*> const& massive_trajectories,
    std::vector<Trajectory<InertialFrame> const*> const& massless_trajectories,
    Instant const& reference_time,
    Time const& t,
    std::vector<Length> const& q,
    std::vector<Acceleration>* result) {
  result->assign(result->size(), Acceleration());
  size_t const number_of_massive_trajectories = massive_trajectories.size();
  size_t const number_of_massless_trajectories = massless_trajectories.size();

  // Declaring variables for values like 3 * b1 + 1, 3 * b2 + 1, etc. in the
  // code below brings no performance advantage as it seems that the compiler is
  // smart enough to figure common subexpressions.
  for (std::size_t b1 = 0, three_b1 = 0;
       b1 < number_of_massive_trajectories;
       ++b1, three_b1 += 3) {
    GravitationalParameter const& body1_gravitational_parameter =
        massive_trajectories[b1]->body().gravitational_parameter();
    for (std::size_t b2 = b1 + 1; b2 < massive_trajectories.size(); ++b2) {
      std::size_t const three_b2 = 3 * b2;
      Length const Δq0 = q[three_b1] - q[three_b2];
      Length const Δq1 = q[three_b1 + 1] - q[three_b2 + 1];
      Length const Δq2 = q[three_b1 + 2] - q[three_b2 + 2];

      Exponentiation<Length, 2> const squared_distance =
          Δq0 * Δq0 + Δq1 * Δq1 + Δq2 * Δq2;
      Exponentiation<Length, -3> const multiplier =
          Sqrt(squared_distance) / (squared_distance * squared_distance);

      auto const μ2OverRSquared =
          massive_trajectories[b2]->body().gravitational_parameter() *
          multiplier;
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
    for (size_t b2 = number_of_massive_trajectories;
         b2 < number_of_massive_trajectories + number_of_massless_trajectories;
         ++b2) {
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

  // Finally, take into account the intrinsic accelerations.
  for (size_t b2 = number_of_massive_trajectories;
       b2 < number_of_massive_trajectories + number_of_massless_trajectories;
       ++b2) {
    std::size_t const three_b2 = 3 * b2;
    Trajectory<InertialFrame> const* trajectory =
        massless_trajectories[b2 - number_of_massive_trajectories];
    if (trajectory->has_intrinsic_acceleration()) {
      R3Element<Acceleration> const acceleration =
          trajectory->evaluate_intrinsic_acceleration(
              t + reference_time).coordinates();
      (*result)[three_b2] += acceleration.x;
      (*result)[three_b2 + 1] += acceleration.y;
      (*result)[three_b2 + 2] += acceleration.z;
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
