#pragma once

#include <memory>
#include <vector>

#include "integrators/symplectic_integrator.hpp"
#include "physics/body.hpp"
#include "quantities/quantities.hpp"

using principia::integrators::SymplecticIntegrator;
using principia::quantities::Time;

namespace principia {
namespace physics {

template<typename InertialFrame>
class NBodySystem {
 public:
  // Takes ownership of the set and the bodies.
  explicit NBodySystem(std::vector<Body<InertialFrame>*> const* bodies);
  ~NBodySystem();

  // The |integrator| must already have been initialized.
  void Integrate(SymplecticIntegrator const& integrator,
                 Time const& tmax,
                 Time const& Δt,
                 int const sampling_period);

 private:
  void ComputeGravitationalAccelerations(double const t,
                                         std::vector<double> const& q,
                                         std::vector<double>* result) const;
  static void ComputeGravitationalVelocities(std::vector<double> const& p,
                                             std::vector<double>* result);

  std::unique_ptr<std::vector<Body<InertialFrame>*> const> const bodies_;
};

}  // namespace physics
}  // namespace principia

#include "physics/n_body_system_body.hpp"
