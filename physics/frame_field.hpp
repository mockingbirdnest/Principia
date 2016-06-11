
#pragma once

#include <functional>

#include "geometry/named_quantities.hpp"

namespace principia {
namespace physics {
namespace internal_frame_field {

using geometry::Instant;
using geometry::Position;
using geometry::Rotation;

// A section of the frame bundle of the manifold |Position|, i.e., a smooth
// assignment of an orthonormal basis to the tangent space of positions at every
// |Position| q.  The orthonormal basis is described as a rotation of the
// standard basis of |Frame|.
template<typename Frame>
using FrameField =
    std::function<Rotation<Frame, Frame>(Position<Frame> const& q)>;

// The identity everywhere.
template<typename Frame>
FrameField<Frame> CoordinateFrame();

}  // namespace internal_frame_field

using internal_frame_field::FrameField;
using internal_frame_field::CoordinateFrame;

}  // namespace physics
}  // namespace principia

#include "physics/frame_field_body.hpp"
