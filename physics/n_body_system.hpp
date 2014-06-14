#pragma once

#include <memory>
#include <vector>

#include "integrators/symplectic_integrator.hpp"
#include "physics/body.hpp"
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
  typedef std::vector<std::unique_ptr<Body<InertialFrame>>> Bodies;

  // Takes ownership of the vectors.  Either pointer may be null.
  // TODO(phl): We would prefer to pass the unique_ptr<> by value, but that
  // confuses the compiler.  So for now, we'll use r-value references.
  NBodySystem(std::unique_ptr<Bodies>&& massive_bodies,
              std::unique_ptr<Bodies>&& massless_bodies);
  ~NBodySystem() = default;

  // No transfer of ownership.
  std::vector<Body<InertialFrame> const*> massive_bodies() const;
  std::vector<Body<InertialFrame> const*> massless_bodies() const;
  std::vector<Body<InertialFrame> const*> bodies() const;

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

  std::unique_ptr<Bodies const> const massive_bodies_;
  std::unique_ptr<Bodies const> const massless_bodies_;

  // The pointers are not owned.  The massive bodies come first.
  std::vector<Body<InertialFrame>*> bodies_;
};

}  // namespace physics
}  // namespace principia

#include "physics/n_body_system_body.hpp"
