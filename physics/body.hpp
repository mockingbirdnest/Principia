#ifndef PRINCIPIA_PHYSICS_BODY_HPP_
#define PRINCIPIA_PHYSICS_BODY_HPP_

#include "base/not_null.hpp"
#include "serialization/physics.pb.h"

using principia::base::not_null;

namespace principia {
namespace physics {

class Body {
 public:
  ~Body() = default;

  // Returns true iff this body is massless.
  virtual bool is_massless() const = 0;

  // Returns true iff this body is oblate (which implies massive).
  virtual bool is_oblate() const = 0;

  // Returns true iff this body is compatible with the given frame (either
  // because it is spherical or because its axis is expressed in the same
  // frame).
  template<typename Frame>
  bool is_compatible_with() const;

  virtual void WriteToMessage(not_null<serialization::Body*> message) const = 0;

  // Dispatches to |MassiveBody| or |MasslessBody| depending on the contents of
  // the message.
  static not_null<std::unique_ptr<Body>> ReadFromMessage(
      serialization::Body const& message);

 protected:
  Body() = default;

 private:
  // A helper class which is here just so that we can specialize it on
  // |is_inertial|.  This is necessary because we cannot dynamic cast to
  // OblateBody<Frame> if |Frame| is not inertial.
  template<typename Frame, bool is_inertial>
  class CompatibilityHelper {
   public:
     static bool is_compatible_with(not_null<Body const*> const body);
   private:
     CompatibilityHelper() = delete;
  };
};

}  // namespace physics
}  // namespace principia

#include "physics/body_body.hpp"

#endif  // PRINCIPIA_PHYSICS_BODY_HPP_
