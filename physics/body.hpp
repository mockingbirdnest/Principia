#pragma once

#include <vector>

#include "geometry/grassmann.hpp"
#include "quantities/named_quantities.hpp"

using principia::quantities::GravitationalParameter;
using principia::quantities::Mass;

namespace principia {
namespace physics {

class Body {
 public:
  // We use the gravitational parameter μ = G M in order not to accumulate
  // unit roundoffs from repeated multiplications by G.
  explicit Body(GravitationalParameter const& gravitational_parameter);
  explicit Body(Mass const& mass);
  ~Body() = default;

  // Returns the construction parameter.
  GravitationalParameter const& gravitational_parameter() const;
  Mass const& mass() const;

  // Returns true iff |gravitational_parameter| (or |mass|) returns 0.
  bool is_massless() const;

 private:
  GravitationalParameter const gravitational_parameter_;
  Mass const mass_;
};

}  // namespace physics
}  // namespace principia

#include "physics/body_body.hpp"
