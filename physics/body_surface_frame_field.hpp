#pragma once

#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/rotation.hpp"
#include "geometry/space.hpp"
#include "physics/ephemeris.hpp"
#include "physics/frame_field.hpp"
#include "physics/rotating_body.hpp"

namespace principia {
namespace physics {
namespace _body_surface_frame_field {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_rotation;
using namespace principia::geometry::_space;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_frame_field;
using namespace principia::physics::_rotating_body;

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

}  // namespace internal

using internal::BodySurfaceFrameField;

}  // namespace _body_surface_frame_field
}  // namespace physics
}  // namespace principia

#include "physics/body_surface_frame_field_body.hpp"
