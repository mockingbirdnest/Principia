#pragma once

#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/rotation.hpp"
#include "physics/ephemeris.hpp"
#include "physics/frame_field.hpp"
#include "physics/rotating_body.hpp"

namespace principia {
namespace physics {
namespace internal_body_surface_frame_field {

using base::not_null;
using geometry::Instant;
using geometry::Position;
using geometry::Rotation;
using geometry::Vector;

// The z-axis goes from the point |q| to the centre of |body| at |t|  The
// x-axis is orthogonal to the z-axis and in the plane defined by the z-axis and
// the polar axis of |body|; its angle with the polar axis is less than π/2 in
// absolute value.  The y-axis is chosen to form a direct basis.  In other
// words, the frame field is (north, east, nadir).
template<typename Frame, typename ThisFrame>
class BodySurfaceFrameField : public FrameField<Frame, ThisFrame> {
 public:
  BodySurfaceFrameField(Ephemeris<Frame> const& ephemeris,
                        Instant const& t,
                        not_null<RotatingBody<Frame> const*> body);

  Rotation<ThisFrame, Frame> FromThisFrame(
      Position<Frame> const& q) const override;

 private:
  Vector<double, Frame> const body_axis_;
  Position<Frame> const body_position_;
};

}  // namespace internal_body_surface_frame_field

using internal_body_surface_frame_field::BodySurfaceFrameField;

}  // namespace physics
}  // namespace principia

#include "physics/body_surface_frame_field_body.hpp"
