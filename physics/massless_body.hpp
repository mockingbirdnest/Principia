// The files containing the tree of child classes of `Body` must be included in
// the order of inheritance to avoid circular dependencies.
#ifndef PRINCIPIA_PHYSICS_BODY_HPP_
#include "physics/body.hpp"
#endif  // PRINCIPIA_PHYSICS_BODY_HPP
#ifndef PRINCIPIA_PHYSICS_MASSLESS_BODY_HPP_
#define PRINCIPIA_PHYSICS_MASSLESS_BODY_HPP_

#include <memory>
#include <vector>

#include "base/not_null.hpp"
#include "serialization/physics.pb.h"

namespace principia {
namespace physics {
namespace _massless_body {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::physics::_body;

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

  // `message.has_massless_body()` must be true.
  static not_null<std::unique_ptr<MasslessBody>> ReadFromMessage(
      serialization::Body const& message);

  static not_null<std::unique_ptr<MasslessBody>> ReadFromMessage(
      serialization::MasslessBody const& message);
};

}  // namespace internal

using internal::MasslessBody;

}  // namespace _massless_body
}  // namespace physics
}  // namespace principia

#include "physics/massless_body_body.hpp"

#endif  // PRINCIPIA_PHYSICS_MASSLESS_BODY_HPP_
