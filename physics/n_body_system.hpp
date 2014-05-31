#pragma once

#include <set>

#include "integrators/symplectic_integrator.hpp"
#include "physics/body.hpp"
#include "quantities/quantities.hpp"

using principia::integrators::SymplecticIntegrator;
using principia::quantities::Time;

namespace principia {
namespace physics {

template<typename Frame>
class NBodySystem {
 public:
  // Takes ownership of the set and the bodies.
  explicit NBodySystem(std::set<Body<Frame>*> const* bodies);
  ~NBodySystem();

  //TODO(phl): What about errors?
  // The |integrator| must already have been initialized.
  void Integrate(SymplecticIntegrator const& integrator,
                 Time const& tmax;
                 Time const& Δt);

 private:
  std::unique_ptr<std::set<Body<Frame>> const bodies_;
};

}  // namespace physics
}  // namespace principia
