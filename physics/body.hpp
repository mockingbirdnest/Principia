#pragma once

#include "geometry/grassmann.hpp"
#include "quantities/named_quantities.hpp"

using principia::geometry::Vector;
using principia::quantities::GravitationalParameter;
using principia::quantities::Length;
using principia::quantities::Momentum;
using principia::quantities::Time;

namespace principia {
namespace physics {

// TODO(phl): The frame used for the positions/momenta of the body.  How do we
// reify frame change and virtual forces?
template<typename Frame>
class Body {
 public:
  // We use the gravitational parameter μ = G M in order not to accumulate
  // unit roundoffs from repeated multiplications by G. Note that in KSP, the
  // gravitational parameter is computed from the mass as G M, but the mass is
  // itself computed from the radius and acceleration due to gravity at sea
  // level as M = g0 r^2 / G. This is silly (and introduces an---admittedly
  // tiny---error), so the gravitational parameter should ideally be computed
  // by the user as μ = g0 r^2. The generally accepted value for g0 in KSP
  // seems to be 9.81 m/s^2.
  explicit Body(GravitationalParameter const& gravitational_parameter);
  ~Body();

  void SetInitial(Vector<Length, Frame> const& position,
                  Vector<Momentum, Frame> const& momentum,
                  Time const& when);

  void AppendToTrajectory(std::vector<Vector<Length, Frame>> const& positions,
                          std::vector<Vector<Momentum, Frame>> const& momenta,
                          std::vector<Time> const& times);
};

}  // namespace physics
}  // namespace principia
