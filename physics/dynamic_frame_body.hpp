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
std::unique_ptr<DynamicFrame<InertialFrame, Frenet<ThisFrame>>>
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
  Rotation<ThisFrame, Frenet<ThisFrame>> to_frenet_trihedron(
      R3x3Matrix(tangent.coordinates(),
                 normal.coordinates(),
                 binormal.coordinates()));

  AngularVelocity<ThisFrame> angular_velocity =
      binormal * Sqrt(InnerProduct(acceleration, acceleration) /
                      InnerProduct(velocity, velocity)) * Radian;

  Velocity<Frenet<ThisFrame>> this_frame_velocity(
      {-velocity.Norm(),
       0 * Metre / Second,
       0 * Metre / Second});
  AngularVelocity<Frenet<ThisFrame>> this_frame_angular_velocity(
      {0 * Radian / Second,
       0 * Radian / Second,
       -Sqrt(InnerProduct(acceleration, acceleration) /
             InnerProduct(velocity, velocity)) * Radian});
  RigidTransformation<ThisFrame, Frenet<ThisFrame>> rigid_transformation(
      degrees_of_freedom.position(),
      Frenet<ThisFrame>::origin,
      to_frenet_trihedron.Forget());
  RigidMotion<ThisFrame, Frenet<ThisFrame>> rigid_motion(
      rigid_transformation,
      this_frame_angular_velocity,
      this_frame_velocity);
  auto rigid_motion_inverse = rigid_motion.Inverse();
  auto geometric_acceleration = [this, to_frenet_trihedron, acceleration,
                                 rigid_motion_inverse](
      Instant const& t,
      DegreesOfFreedom<Frenet<ThisFrame>> const& dof) {
    return to_frenet_trihedron(
        GeometricAcceleration(t, rigid_motion_inverse(dof)) - acceleration);
  };
  std::unique_ptr<DynamicFrame<InertialFrame, Frenet<ThisFrame>>> result =
      std::make_unique<
          InstantaneouslyDefinedFrame<InertialFrame, Frenet<InertialFrame>>>(
          t /*definition_instant*/,
          rigid_motion * ToThisFrameAtTime(t) /*to_this_frame*/,
          geometric_acceleration);
  return std::move(result);
}

template<typename InertialFrame, typename ThisFrame>
InstantaneouslyDefinedFrame<InertialFrame, ThisFrame>::
    InstantaneouslyDefinedFrame(
        Instant const& definition_instant,
        RigidMotion<InertialFrame, ThisFrame> const& to_this_frame,
        std::function<Vector<Acceleration, ThisFrame>(
            Instant const& t,
            DegreesOfFreedom<ThisFrame> const& degrees_of_freedom)>
            geometric_acceleration)
    : definition_instant_(definition_instant),
      to_this_frame_(to_this_frame),
      geometric_acceleration_(geometric_acceleration) {}

template<typename InertialFrame, typename ThisFrame>
RigidMotion<InertialFrame, ThisFrame>
InstantaneouslyDefinedFrame<InertialFrame, ThisFrame>::ToThisFrameAtTime(
    Instant const& t) const {
  CHECK_EQ(definition_instant_, t);
  return to_this_frame_;
}

template<typename InertialFrame, typename ThisFrame>
RigidMotion<ThisFrame, InertialFrame>
InstantaneouslyDefinedFrame<InertialFrame, ThisFrame>::FromThisFrameAtTime(
    Instant const& t) const {
  CHECK_EQ(definition_instant_, t);
  return to_this_frame_.Inverse();
}

template<typename InertialFrame, typename ThisFrame>
Vector<Acceleration, ThisFrame>
InstantaneouslyDefinedFrame<InertialFrame, ThisFrame>::GeometricAcceleration(
    Instant const& t,
    DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const {
  CHECK_EQ(definition_instant_, t);
  return geometric_acceleration_(t, degrees_of_freedom);
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
          OtherFrame::origin,
          origin_at_epoch_ + (t - epoch_) * velocity_,
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
