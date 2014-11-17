#pragma once

#include <vector>

#include "geometry/grassmann.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

using principia::geometry::Vector;
using principia::quantities::GravitationalParameter;
using principia::quantities::Length;
using principia::quantities::Mass;
using principia::quantities::Order2ZonalCoefficient;

namespace principia {
namespace physics {

template<typename Frame>
class Body {
 public:
  // We use the gravitational parameter μ = G M in order not to accumulate
  // unit roundoffs from repeated multiplications by G.
  explicit Body(GravitationalParameter const& gravitational_parameter);
  explicit Body(Mass const& mass);

  // A body with oblateness.  The body must not be massless.  The frame must be
  // inertial to ensure that the axis is well defined.  These constructors are
  // templatized just to enable SFINAE, clients should let the template
  // parameter default.
  template<typename F = Frame>
  Body(GravitationalParameter const& gravitational_parameter,
       double const j2,
       Length const& radius,
       std::enable_if_t<F::is_inertial, Vector<double, F>> const& axis);
  template<typename F = Frame>
  Body(Mass const& mass,
       double const j2,
       Length const& radius,
       std::enable_if_t<F::is_inertial, Vector<double, F>> const& axis);
  template<typename F = Frame>
  Body(GravitationalParameter const& gravitational_parameter,
       Order2ZonalCoefficient const& j2,
       std::enable_if_t<F::is_inertial, Vector<double, F>> const& axis);
  template<typename F = Frame>
  Body(Mass const& mass,
       Order2ZonalCoefficient const& j2,
       std::enable_if_t<F::is_inertial, Vector<double, F>> const& axis);

  ~Body() = default;

  // Returns the construction parameter.
  GravitationalParameter const& gravitational_parameter() const;
  Mass const& mass() const;

  // Returns the j2 coefficient.  Returns 0 for a non-oblate body.
  Order2ZonalCoefficient const& j2() const;

  // Returns the axis passed at construction.
  Vector<double, Frame> const& axis() const;

  // Returns true iff |gravitational_parameter| (or |mass|) returns 0.
  bool is_massless() const;

  // Returns true iff |j2| returns a non-0 value.
  bool is_oblate() const;

 private:
  GravitationalParameter const gravitational_parameter_;
  Mass const mass_;
  Order2ZonalCoefficient const j2_;
  Vector<double, Frame> const axis_;
};

}  // namespace physics
}  // namespace principia

#include "physics/body_body.hpp"
