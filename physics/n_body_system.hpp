#pragma once

#include <memory>
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
  typedef std::vector<std::unique_ptr<Body>> Bodies;
  typedef std::vector<std::unique_ptr<Trajectory<InertialFrame>>> Trajectories;

  // TODO(phl): Unclear relation between trajectories and bodies_.
  NBodySystem(Bodies&& massive_bodies,
              Bodies&& massless_bodies,
              Trajectories&& trajectories);
  ~NBodySystem() = default;

  // No transfer of ownership.
  std::vector<Body const*> massive_bodies() const;
  std::vector<Body const*> massless_bodies() const;
  std::vector<Body const*> bodies() const;
  std::vector<Trajectory<InertialFrame> const*> trajectories() const;

  // The |integrator| must already have been initialized.
  void Integrate(SymplecticIntegrator<Length, Speed> const& integrator,
                 Time const& tmax,
                 Time const& Δt,
                 int const sampling_period);

 private:
  void ComputeGravitationalAccelerations(
      Time const& t,
      std::vector<Length> const& q,
      std::vector<Acceleration>* result) const;
  static void ComputeGravitationalVelocities(std::vector<Speed> const& p,
                                             std::vector<Speed>* result);

  Bodies const massive_bodies_;
  Bodies const massless_bodies_;
  Trajectories const trajectories_;

  // The pointers are not owned.  The massive bodies come first.
  std::vector<Body*> bodies_;
};

}  // namespace physics
}  // namespace principia

#include "physics/n_body_system_body.hpp"
