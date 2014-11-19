#pragma once

#include "physics/body.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

using principia::quantities::GravitationalParameter;
using principia::quantities::Mass;

namespace principia {
namespace physics {

class MassiveBody : public Body {
 public:
  // We use the gravitational parameter μ = G M in order not to accumulate
  // unit roundoffs from repeated multiplications by G.  The parameter must not
  // be zero.
  explicit MassiveBody(GravitationalParameter const& gravitational_parameter);
  explicit MassiveBody(Mass const& mass);
  ~MassiveBody() = default;

  // Returns the construction parameter.
  GravitationalParameter const& gravitational_parameter() const;
  Mass const& mass() const;

  // Returns false.
  bool is_massless() const override;

  // Returns false.
  bool is_oblate() const override;

 private:
  GravitationalParameter const gravitational_parameter_;
  Mass const mass_;
};

}  // namespace physics
}  // namespace principia

#include "physics/massive_body_body.hpp"
