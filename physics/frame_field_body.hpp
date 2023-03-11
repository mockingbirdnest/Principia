#pragma once

#include <functional>

#include "geometry/named_quantities.hpp"

namespace principia {
namespace physics {
namespace _frame_field {
namespace internal {

template<typename Frame, typename ThisFrame>
Rotation<ThisFrame, Frame> FrameField<Frame, ThisFrame>::FromThisFrame(
    Position<Frame> const& q) const {
  return ToThisFrame(q).Inverse();
}

template<typename Frame, typename ThisFrame>
Rotation<Frame, ThisFrame> FrameField<Frame, ThisFrame>::ToThisFrame(
    Position<Frame> const& q) const {
  return FromThisFrame(q).Inverse();
}

template<typename Frame, typename ThisFrame>
Rotation<ThisFrame, Frame>
CoordinateFrameField<Frame, ThisFrame>::FromThisFrame(
    Position<Frame> const& q) const {
  return Rotation<ThisFrame, Frame>::Identity();
}

}  // namespace internal
}  // namespace _frame_field
}  // namespace physics
}  // namespace principia
