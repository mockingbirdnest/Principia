#pragma once

#include <vector>

#include "geometry/grassmann.hpp"
#include "quantities/named_quantities.hpp"

using principia::geometry::Vector;
using principia::quantities::GravitationalParameter;
using principia::quantities::Length;
using principia::quantities::Mass;
using principia::quantities::Speed;
using principia::quantities::Time;

namespace principia {
namespace physics {

// TODO(phl): The frame used for the positions/momenta of the body.  Not sure if
// this API is right: it might be more logical to templatize the methods that
// deal with positions and speed (but then we would need to understand frame
// changes).
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
  explicit Body(Mass const& mass);
  ~Body();

  // Returns the construction parameter.
  GravitationalParameter const& gravitational_parameter() const;
  Mass const& mass() const;

  // Returns true iff |gravitational_parameter| (or |mass|) returns 0.
  bool is_massless() const;

  // Appends one point to the trajectory of the body.
  void AppendToTrajectory(Vector<Length, Frame> const& position,
                          Vector<Speed, Frame> const& velocity,
                          Time const& time);

  // These functions return the series of positions/velocities/times for the
  // trajectory of the body.  All three vectors are guaranteed to have the same
  // length.
  std::vector<Vector<Length, Frame>> const& positions() const;
  std::vector<Vector<Speed, Frame>> const& velocities() const;
  std::vector<Time> const& times() const;

 private:
  GravitationalParameter const gravitational_parameter_;
  Mass const mass_;

  std::vector<Vector<Length, Frame>> positions_;
  std::vector<Vector<Speed, Frame>> velocities_;
  std::vector<Time> times_;
};

}  // namespace physics
}  // namespace principia

#include "physics/body_body.hpp"
