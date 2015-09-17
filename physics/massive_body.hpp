
// The files containing the tree of child classes of |Body| must be included in
// the order of inheritance to avoid circular dependencies.  This class will end
// up being reincluded as part of the implementation of its parent.
#ifndef PRINCIPIA_PHYSICS_BODY_HPP_
#include "physics/body.hpp"
#else
#ifndef PRINCIPIA_PHYSICS_MASSIVE_BODY_HPP_
#define PRINCIPIA_PHYSICS_MASSIVE_BODY_HPP_

#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

namespace principia {

using quantities::GravitationalParameter;
using quantities::Mass;

namespace physics {

class MassiveBody : public Body {
 public:
  // We use the gravitational parameter μ = G M in order not to accumulate
  // unit roundoffs from repeated multiplications by G.  The parameter must not
  // be zero.
  class Parameters {
   public:
    // The constructors are implicit on purpose.
    Parameters(
        GravitationalParameter const& gravitational_parameter);  // NOLINT
    Parameters(Mass const& mass);  // NOLINT(runtime/explicit)

   private:
    GravitationalParameter const gravitational_parameter_;
    Mass const mass_;
    friend class MassiveBody;
  };

  explicit MassiveBody(Parameters const& parameters);
  ~MassiveBody() = default;

  // Returns the construction parameter.
  GravitationalParameter const& gravitational_parameter() const;
  Mass const& mass() const;

  // Returns false.
  bool is_massless() const override;

  // Returns false.
  bool is_oblate() const override;

  void WriteToMessage(not_null<serialization::Body*> message) const override;

  virtual void WriteToMessage(
      not_null<serialization::MassiveBody*> message) const;

  // Both methods below dispatch to |OblateBody<UnknownFrame>| if the
  // |OblateBody| extension is present in the message.  Use |reinterpret_cast|
  // afterwards as appropriate if the frame is known.

  // |message.has_massless_body()| must be true.
  static not_null<std::unique_ptr<MassiveBody>> ReadFromMessage(
      serialization::Body const& message);

  static not_null<std::unique_ptr<MassiveBody>> ReadFromMessage(
      serialization::MassiveBody const& message);

 private:
  Parameters const parameters_;
};

}  // namespace physics
}  // namespace principia

#include "physics/massive_body_body.hpp"

#endif  // PRINCIPIA_PHYSICS_MASSIVE_BODY_HPP_
#endif  // PRINCIPIA_PHYSICS_BODY_HPP_
