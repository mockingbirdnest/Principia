#pragma once

#include <memory>
#include <set>
#include <vector>

#include "integrators/symplectic_integrator.hpp"
#include "physics/body.hpp"
#include "physics/trajectory.hpp"
#include "quantities/quantities.hpp"

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
  typedef std::vector<std::unique_ptr<Body> const> Bodies;
  typedef std::vector<Trajectory<InertialFrame>*> Trajectories;  // Not owned.

  // TODO(phl): Unclear relation between trajectories and bodies_.
  NBodySystem(Bodies&& massive_bodies,
              Bodies&& massless_bodies);
  ~NBodySystem() = default;

  // No transfer of ownership.
  std::vector<Body const*> massive_bodies() const;
  std::vector<Body const*> massless_bodies() const;
  std::vector<Body const*> bodies() const;

  // The |integrator| must already have been initialized.  All the
  // |trajectories| must have the same last_time() and must be for bodies passed
  // at construction.
  void Integrate(SymplecticIntegrator<Length, Speed> const& integrator,
                 Time const& tmax,
                 Time const& Δt,
                 int const sampling_period,
                 Trajectories const& trajectories);

 private:
  static void ComputeGravitationalAccelerations(
      std::vector<Trajectory<InertialFrame> const*> const& massive_trajectories,
      int const number_of_massless_trajectories,
      Time const& t,
      std::vector<Length> const& q,
      std::vector<Acceleration>* result);
  static void ComputeGravitationalVelocities(std::vector<Speed> const& p,
                                             std::vector<Speed>* result);

  Bodies const massive_bodies_;
  Bodies const massless_bodies_;

  // The pointers are not owned.
  std::set<Body const*> bodies_;
};

}  // namespace physics
}  // namespace principia

#include "physics/n_body_system_body.hpp"
