#pragma once

namespace principia {
namespace physics {

class Body {
 public:
  ~Body() = default;

  // Returns true iff this body is massless.
  virtual bool is_massless() const = 0;

  // Returns true iff this body is oblate (which implies massive).
  virtual bool is_oblate() const = 0;

 protected:
  Body() = default;
};

}  // namespace physics
}  // namespace principia
