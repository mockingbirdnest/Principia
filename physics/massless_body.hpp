#pragma once

#include <vector>

#include "physics/body.hpp

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
};

}  // namespace physics
}  // namespace principia

#include "physics/massless_body_body.hpp"
