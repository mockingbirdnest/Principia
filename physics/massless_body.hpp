// The files containing the tree of of child classes of |Body| must be included
// in the order of inheritance to avoid circular dependencies.  This class will
// end up being reincluded as part of the implementation of its parent.
#ifndef PRINCIPIA_PHYSICS_BODY_HPP_
#include "physics/body.hpp"
#else
#ifndef PRINCIPIA_PHYSICS_MASSLESS_BODY_HPP_
#define PRINCIPIA_PHYSICS_MASSLESS_BODY_HPP_

#include <vector>

#include "physics/body.hpp"

namespace principia {
namespace physics {

class MasslessBody : public Body {
 public:
  MasslessBody() = default;
  ~MasslessBody() = default;

  // Returns true.
  bool is_massless() const override;

  // Returns false.
  bool is_oblate() const override;

  void WriteToMessage(not_null<serialization::Body*> message) const override;
};

}  // namespace physics
}  // namespace principia

#include "physics/massless_body_body.hpp"

#endif  // PRINCIPIA_PHYSICS_MASSLESS_BODY_HPP_
#endif  // PRINCIPIA_PHYSICS_BODY_HPP
