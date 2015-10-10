#pragma once

#include "physics/dynamic_frame.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {

using geometry::Bivector;
using geometry::InnerProduct;
using geometry::Normalize;
using geometry::R3x3Matrix;
using geometry::Wedge;
using quantities::Sqrt;
using quantities::si::Metre;
using quantities::si::Second;

namespace physics {

template <typename InertialFrame, typename ThisFrame>
Rotation<Frenet<ThisFrame>, ThisFrame>
DynamicFrame<InertialFrame, ThisFrame>::FrenetFrame(
    Instant const& t,
    DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const {
  Velocity<ThisFrame> const& velocity = degrees_of_freedom.velocity();
  Vector<Acceleration, ThisFrame> const acceleration =
      GeometricAcceleration(t, degrees_of_freedom);
  Vector<Acceleration, ThisFrame> normal_acceleration = acceleration;
  velocity.Orthogonalize<Acceleration>(&normal_acceleration);
  Vector<double, ThisFrame> tangent = Normalize(velocity);
  Vector<double, ThisFrame> normal = Normalize(normal_acceleration);
  Bivector<double, ThisFrame> binormal = Wedge(tangent, normal);
  // Maps |tangent| to {1, 0, 0}, |normal| to {0, 1, 0}, and |binormal| to
  // {0, 0, 1}.
  return Rotation<Frenet<ThisFrame>, ThisFrame>(
      R3x3Matrix(tangent.coordinates(),
                 normal.coordinates(),
                 binormal.coordinates()).Transpose());
}

template<typename OtherFrame, typename ThisFrame>
InertialFrame<OtherFrame, ThisFrame>::InertialFrame(
    Velocity<OtherFrame> const& velocity,
    Position<OtherFrame> const& origin_at_epoch,
    Instant const& epoch,
    OrthogonalMap<OtherFrame, ThisFrame> const& orthogonal_map,
    std::function<Vector<Acceleration, OtherFrame>(
        Instant const& t,
        Position<OtherFrame> const& q)> gravity)
    : velocity_(velocity),
      origin_at_epoch_(origin_at_epoch),
      epoch_(epoch),
      orthogonal_map_(orthogonal_map),
      gravity_(gravity) {}

template<typename OtherFrame, typename ThisFrame>
RigidMotion<OtherFrame, ThisFrame>
InertialFrame<OtherFrame, ThisFrame>::ToThisFrameAtTime(
    Instant const& t) const {
  return RigidMotion<OtherFrame, ThisFrame>(
      RigidTransformation<OtherFrame, ThisFrame>(
          origin_at_epoch_ + (t - epoch_) * velocity_,
          ThisFrame::origin,
          orthogonal_map_),
      AngularVelocity<ThisFrame>(), -orthogonal_map_(velocity_));
}

template<typename OtherFrame, typename ThisFrame>
RigidMotion<ThisFrame, OtherFrame>
InertialFrame<OtherFrame, ThisFrame>::FromThisFrameAtTime(
    Instant const& t) const {
  return ToThisFrameAtTime(t).Inverse();
}

template<typename OtherFrame, typename ThisFrame>
Vector<Acceleration, ThisFrame>
InertialFrame<OtherFrame, ThisFrame>::GeometricAcceleration(
    Instant const& t,
    DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const {
  return orthogonal_map_(
      gravity_(t, FromThisFrameAtTime(t)(degrees_of_freedom).position()));
}


}  // namespace physics
}  // namespace principia
