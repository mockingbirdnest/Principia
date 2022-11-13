#pragma once

#include "physics/body_surface_frame_field.hpp"

#include "physics/continuous_trajectory.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace physics {
namespace internal_body_surface_frame_field {

using geometry::Bivector;
using geometry::Displacement;
using geometry::InnerProduct;
using geometry::Normalize;
using geometry::Wedge;
using quantities::Length;
using quantities::Sqrt;

template<typename Frame, typename ThisFrame>
BodySurfaceFrameField<Frame, ThisFrame>::BodySurfaceFrameField(
    Ephemeris<Frame> const& ephemeris,
    Instant const& t,
    not_null<RotatingBody<Frame> const*> const body)
    : body_axis_(body->polar_axis()),
      body_position_(
          ephemeris.trajectory(body)->EvaluatePosition(t)) {}

template<typename Frame, typename ThisFrame>
Rotation<ThisFrame, Frame>
BodySurfaceFrameField<Frame, ThisFrame>::FromThisFrame(
    Position<Frame> const& q) const {
  Displacement<Frame> const from_body_centre = q - body_position_;

  Vector<double, Frame> const zenith = Normalize(from_body_centre);
  double const axis_projection =
      InnerProduct(zenith, body_axis_);
  double const axis_projection² = axis_projection * axis_projection;

  // The unit vector |north| is directed along the polar axis.  Note that λ is
  // positive.
  double const λ = 1 / Sqrt(1 - axis_projection²);
  auto const μ = -λ * axis_projection;
  Vector<double, Frame> const north = λ * body_axis_ + μ * zenith;
  Vector<double, Frame> const nadir = -zenith;
  Bivector<double, Frame> const east = Wedge(nadir, north);

  return Rotation<ThisFrame, Frame>(north, east, nadir);
}

}  // namespace internal_body_surface_frame_field
}  // namespace physics
}  // namespace principia
