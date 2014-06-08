#pragma once

#include <memory>
#include <vector>

#include "integrators/symplectic_integrator.hpp"
#include "physics/body.hpp"
#include "physics/frame.hpp"
#include "quantities/quantities.hpp"

using principia::integrators::SymplecticIntegrator;
using principia::quantities::Acceleration;
using principia::quantities::Length;
using principia::quantities::Speed;
using principia::quantities::Time;

namespace principia {
namespace physics {

class NBodySystem {
 public:
  // Takes ownership of the set and the bodies.
  explicit NBodySystem(std::vector<Body<InertialFrame>*> const* bodies);
  ~NBodySystem();

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

  std::unique_ptr<std::vector<Body<InertialFrame>*> const> const bodies_;
};

}  // namespace physics
}  // namespace principia

#include "physics/n_body_system_body.hpp"
