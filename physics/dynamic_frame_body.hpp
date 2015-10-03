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
DynamicFrame<InertialFrame, Frenet<ThisFrame>>
DynamicFrame<InertialFrame, ThisFrame>::FrenetFrame(
    Instant const& t,
    DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const {
  Velocity<ThisFrame> const& velocity = degrees_of_freedom.velocity();
  Vector<Acceleration, ThisFrame> const acceleration =
      GeometricAcceleration(t, degrees_of_freedom);
  Vector<Acceleration, ThisFrame> normal_acceleration = acceleration;
  velocity.Orthogonalize(&normal_acceleration);
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
      to_frenet_trihedron);
  RigidMotion<ThisFrame, Frenet<ThisFrame>> rigid_motion(
      rigid_transformation,
      this_frame_angular_velocity,
      this_frame_velocity);
  return DynamicFrame<InertialFrame, Frenet<ThisFrame>>();
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
  return geometric_acceleration_(degrees_of_freedom);
}

}  // namespace physics
}  // namespace principia
