// The files containing the tree of of child classes of |Body| must be included
// in the order of inheritance to avoid circular dependencies.  This class will
// end up being reincluded as part of the implementation of its parent.
#ifndef PRINCIPIA_PHYSICS_BODY_HPP_
#include "physics/body.hpp"
#else
#ifndef PRINCIPIA_PHYSICS_MASSIVE_BODY_HPP_
#define PRINCIPIA_PHYSICS_MASSIVE_BODY_HPP_

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

  void WriteToMessage(not_null<serialization::Body*> message) const override;

  // Fails if the |MassiveBody| extension is absent from the message.
  static not_null<std::unique_ptr<MassiveBody>> ReadFromMessage(
      serialization::Body const& message);

  // Dispatches to |OblateBody<UncheckedInertialFrame>| if the |OblateBody|
  // extension is present in the message.
  static not_null<std::unique_ptr<MassiveBody>> ReadFromMessage(
      serialization::MassiveBody const& message);

 private:
  GravitationalParameter const gravitational_parameter_;
  Mass const mass_;
};

}  // namespace physics
}  // namespace principia

#include "physics/massive_body_body.hpp"

#endif  // PRINCIPIA_PHYSICS_MASSIVE_BODY_HPP_
#endif  // PRINCIPIA_PHYSICS_BODY_HPP_
