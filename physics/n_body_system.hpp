#pragma once

#include <memory>
#include <set>
#include <vector>

#include "geometry/named_quantities.hpp"
#include "integrators/symplectic_integrator.hpp"
#include "physics/body.hpp"
#include "physics/trajectory.hpp"
#include "quantities/quantities.hpp"

using principia::geometry::Instant;
using principia::integrators::SymplecticIntegrator;
using principia::quantities::Acceleration;
using principia::quantities::Length;
using principia::quantities::Speed;
using principia::quantities::Time;

namespace principia {
namespace physics {

template<typename InertialFrame>
class NBodySystem {
 public:
  using Trajectories = std::vector<Trajectory<InertialFrame>*>;  // Not owned.

  NBodySystem() = default;
  ~NBodySystem() = default;

  // The |integrator| must already have been initialized.  All the
  // |trajectories| must have the same |last_time()| and must be for distinct
  // bodies.
  void Integrate(SymplecticIntegrator<Length, Speed> const& integrator,
                 Instant const& tmax,
                 Time const& Δt,
                 int const sampling_period,
                 bool const tmax_is_exact,
                 Trajectories const& trajectories);

 private:
  using ReadonlyTrajectories = std::vector<Trajectory<InertialFrame> const*>;

  // Computes the acceleration due to one body, |body1| (with index |b1| in the
  // |q| and |result| arrays) on the bodies with indices [b2_begin, b2_end[ in
  // |body2_trajectories|.  The template parameters specify what we know about
  // the bodies, and therefore what forces apply.
  template<bool body1_is_oblate,
           bool body2_is_oblate,
           bool body2_is_massive>
  static void ComputeOneBodyGravitationalAcceleration(
      Body<InertialFrame> const& body1,
      ReadonlyTrajectories const& body2_trajectories,
      size_t const b1,
      size_t const b2_begin,
      size_t const b2_end,
      std::vector<Length> const& q,
      std::vector<Acceleration>* result);

  // No transfer of ownership.
  static void ComputeGravitationalAccelerations(
      ReadonlyTrajectories const& massive_oblate_trajectories,
      ReadonlyTrajectories const& massive_spherical_trajectories,
      ReadonlyTrajectories const& massless_trajectories,
      Instant const& reference_time,
      Time const& t,
      std::vector<Length> const& q,
      std::vector<Acceleration>* result);

  // No transfer of ownership.
  static void ComputeGravitationalVelocities(std::vector<Speed> const& p,
                                             std::vector<Speed>* result);
};

}  // namespace physics
}  // namespace principia

#include "physics/n_body_system_body.hpp"
