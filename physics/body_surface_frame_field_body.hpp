
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

template<typename Frame, typename ThisFrame>
BodySurfaceFrameField<Frame, ThisFrame>::BodySurfaceFrameField(
    Ephemeris<Frame> const& ephemeris,
    Instant const& t,
    not_null<RotatingBody<Frame> const*> const body)
    : body_axis_(body->polar_axis()),
      body_position_(
          ephemeris.trajectory(body)->EvaluatePosition(t, /*hint=*/nullptr)) {}

template<typename Frame, typename ThisFrame>
Rotation<ThisFrame, Frame>
BodySurfaceFrameField<Frame, ThisFrame>::FromThisFrame(
    Position<Frame> const& q) const {
  // A vector going from the centre of the body to the point in the field.
  Displacement<Frame> const displacement = q - body_position_;

  Vector<double, Frame> const normalized_displacement = Normalize(displacement);
  double const axis_projection =
      InnerProduct(normalized_displacement, body_axis_);
  double const axis_projection² = axis_projection * axis_projection;

  // The unit vector x directed along the polar axis is
  // |λ * body_axis + μ * normalized_displacement|.  Note that λ is positive.
  double const λ = 1 / Sqrt(1 + axis_projection²);
  auto const μ = -λ * axis_projection;

  Vector<double, Frame> const x = λ * body_axis_ + μ * normalized_displacement;
  Vector<double, Frame> const& z = normalized_displacement;
  Bivector<double, Frame> const y = Wedge(z, x);

  return Rotation<ThisFrame, Frame>(x, y, z);
}

}  // namespace internal_body_surface_frame_field
}  // namespace physics
}  // namespace principia
