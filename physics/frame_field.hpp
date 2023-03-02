#pragma once

#include <functional>

#include "geometry/named_quantities.hpp"

namespace principia {
namespace physics {
namespace internal_frame_field {

using namespace principia::geometry::_named_quantities;
using namespace principia::geometry::_rotation;

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

}  // namespace internal_frame_field

using internal_frame_field::FrameField;
using internal_frame_field::CoordinateFrameField;

}  // namespace physics
}  // namespace principia

#include "physics/frame_field_body.hpp"
