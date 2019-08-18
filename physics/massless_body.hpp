
// The files containing the tree of child classes of |Body| must be included in
// the order of inheritance to avoid circular dependencies.
#ifndef PRINCIPIA_PHYSICS_BODY_HPP_
#include "physics/body.hpp"
#endif  // PRINCIPIA_PHYSICS_BODY_HPP
#ifndef PRINCIPIA_PHYSICS_MASSLESS_BODY_HPP_
#define PRINCIPIA_PHYSICS_MASSLESS_BODY_HPP_

#include <vector>

#include "serialization/physics.pb.h"

namespace principia {
namespace physics {
namespace internal_massless_body {

using base::not_null;

class MasslessBody : public Body {
 public:
  MasslessBody() = default;

  // Returns true.
  bool is_massless() const override;

  // Returns false.
  bool is_oblate() const override;

  void WriteToMessage(not_null<serialization::Body*> message) const override;

  virtual void WriteToMessage(
      not_null<serialization::MasslessBody*> message) const;

  // |message.has_massless_body()| must be true.
  static not_null<std::unique_ptr<MasslessBody>> ReadFromMessage(
      serialization::Body const& message);

  static not_null<std::unique_ptr<MasslessBody>> ReadFromMessage(
      serialization::MasslessBody const& message);
};

}  // namespace internal_massless_body

using internal_massless_body::MasslessBody;

}  // namespace physics
}  // namespace principia

#include "physics/massless_body_body.hpp"

#endif  // PRINCIPIA_PHYSICS_MASSLESS_BODY_HPP_
