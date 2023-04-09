#pragma once

#include <functional>

#include "geometry/rotation.hpp"
#include "geometry/space.hpp"

namespace principia {
namespace physics {
namespace _frame_field {
namespace internal {

using namespace principia::geometry::_rotation;
using namespace principia::geometry::_space;

// A section of the frame bundle of the manifold |Position|, i.e., a smooth
// assignment of an orthonormal basis to the tangent space of positions at every
// |Position| q.
template<typename Frame, typename ThisFrame>
class FrameField {
 public:
  virtual ~FrameField() = default;

  virtual Rotation<ThisFrame, Frame> FromThisFrame(
      Position<Frame> const& q) const;

  virtual Rotation<Frame, ThisFrame> ToThisFrame(
      Position<Frame> const& q) const;

 protected:
  FrameField() = default;
};

// The identity everywhere.
template<typename Frame, typename ThisFrame>
class CoordinateFrameField : public FrameField<Frame, ThisFrame> {
  Rotation<ThisFrame, Frame> FromThisFrame(
      Position<Frame> const& q) const override;
};

}  // namespace internal

using internal::FrameField;
using internal::CoordinateFrameField;

}  // namespace _frame_field
}  // namespace physics
}  // namespace principia

namespace principia::physics {
using namespace principia::physics::_frame_field;
}  // namespace principia::physics

#include "physics/frame_field_body.hpp"
