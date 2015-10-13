#pragma once

#include "physics/barycentric_rotating_dynamic_frame.hpp"

#include "geometry/named_quantities.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {

using geometry::Displacement;
using geometry::R3x3Matrix;
using geometry::Velocity;
using geometry::Wedge;
using quantities::Length;
using quantities::Pow;
using quantities::Product;
using quantities::Speed;

namespace physics {

template<typename InertialFrame, typename ThisFrame>
BarycentricRotatingDynamicFrame<InertialFrame, ThisFrame>::
BarycentricRotatingDynamicFrame(
    not_null<Ephemeris<InertialFrame> const*> const ephemeris,
    not_null<MassiveBody const*> const primary,
    not_null<MassiveBody const*> const secondary)
    : ephemeris_(ephemeris),
      primary_(primary),
      secondary_(secondary),
      primary_trajectory_(ephemeris_->trajectory(primary_)),
      secondary_trajectory_(ephemeris_->trajectory(secondary_)) {}

template<typename InertialFrame, typename ThisFrame>
RigidMotion<InertialFrame, ThisFrame>
BarycentricRotatingDynamicFrame<InertialFrame, ThisFrame>::ToThisFrameAtTime(
    Instant const& t) const {
  DegreesOfFreedom<InertialFrame> const primary_degrees_of_freedom =
      primary_trajectory_->EvaluateDegreesOfFreedom(t, &primary_hint_);
  DegreesOfFreedom<InertialFrame> const secondary_degrees_of_freedom =
      secondary_trajectory_->EvaluateDegreesOfFreedom(t, &secondary_hint_);
  DegreesOfFreedom<InertialFrame> const barycentre_degrees_of_freedom =
      Barycentre<InertialFrame, GravitationalParameter>(
          {primary_degrees_of_freedom,
           secondary_degrees_of_freedom},
          {primary_->gravitational_parameter(),
           secondary_->gravitational_parameter()});

  Rotation<InertialFrame, ThisFrame>
      from_basis_of_inertial_frame_to_basis_of_this_frame =
          Rotation<InertialFrame, ThisFrame>::Identity();
  AngularVelocity<InertialFrame> angular_velocity;
  FromBasisOfInertialFrameToBasisOfThisFrame(
      barycentre_degrees_of_freedom,
      primary_degrees_of_freedom,
      secondary_degrees_of_freedom,
      &from_basis_of_inertial_frame_to_basis_of_this_frame,
      &angular_velocity);

  RigidTransformation<InertialFrame, ThisFrame> const
      rigid_transformation(centre_degrees_of_freedom.position(),
                           ThisFrame::origin,
                           from_basis_of_inertial_frame_to_basis_of_this_frame);
  return RigidMotion<InertialFrame, ThisFrame>(
             rigid_transformation,
             angular_velocity,
             barycentre_degrees_of_freedom.velocity());
}

template<typename InertialFrame, typename ThisFrame>
RigidMotion<ThisFrame, InertialFrame>
BarycentricRotatingDynamicFrame<InertialFrame, ThisFrame>::
FromThisFrameAtTime(Instant const& t) const {
  Rotation<ThisFrame, InertialFrame>
      from_standard_basis_to_basis_of_last_barycentric_frame =
          Rotation<ThisFrame, InertialFrame>::Identity();
  DegreesOfFreedom<InertialFrame> const& primary_degrees_of_freedom =
      primary_trajectory_->EvaluateDegreesOfFreedom(t, &primary_hint_);
  DegreesOfFreedom<InertialFrame> const& secondary_degrees_of_freedom =
      secondary_trajectory_->EvaluateDegreesOfFreedom(t, &secondary_hint_);
  DegreesOfFreedom<InertialFrame> const barycentre_degrees_of_freedom =
      Barycentre<InertialFrame, GravitationalParameter>(
          {primary_degrees_of_freedom,
           secondary_degrees_of_freedom},
          {primary_->gravitational_parameter(),
           secondary_->gravitational_parameter()});

  Rotation<ThisFrame, InertialFrame>
      from_basis_of_this_frame_to_basic_of_inertial_frame =
          Rotation<ThisFrame, InertialFrame>::Identity();
  AngularVelocity<ThisFrame> angular_velocity;
  FromBasisOfThisFrameToBasisOfInertialFrame<ThisFrame, InertialFrame>(
      barycentre_degrees_of_freedom,
      primary_degrees_of_freedom,
      secondary_degrees_of_freedom,
      &from_basis_of_inertial_frame_to_basis_of_this_frame,
      &angular_velocity);

  RigidTransformation<ThisFrame, InertialFrame> const
      rigid_transformation(
          ThisFrame::origin,
          barycentre_degrees_of_freedom.position(),
          from_basis_of_this_frame_to_basic_of_inertial_frame);
  return RigidMotion<InertialFrame, ThisFrame>(
             rigid_transformation,
             angular_velocity,
             barycentre_degrees_of_freedom.velocity());
}

template<typename InertialFrame, typename ThisFrame>
Vector<Acceleration, ThisFrame>
BarycentricRotatingDynamicFrame<InertialFrame, ThisFrame>::
GeometricAcceleration(
    Instant const& t,
    DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const {
}

template<typename InertialFrame, typename ThisFrame>
void BarycentricRotatingDynamicFrame<InertialFrame, ThisFrame>::
FromBasisOfInertialFrameToBasisOfThisFrame(
    DegreesOfFreedom<InertialFrame> const& barycentre_degrees_of_freedom,
    DegreesOfFreedom<InertialFrame> const& primary_degrees_of_freedom,
    DegreesOfFreedom<InertialFrame> const& secondary_degrees_of_freedom,
    not_null<Rotation<InertialFrame, ThisFrame>*> const rotation,
    not_null<AngularVelocity<InertialFrame>*> const angular_velocity) {
  RelativeDegreesOfFreedom<InertialFrame> const reference =
      primary_degrees_of_freedom - barycentre_degrees_of_freedom;
  Displacement<InertialFrame> const& reference_direction =
      reference.displacement();
  Velocity<ThisFrame> reference_normal = reference.velocity();
  reference_direction.template Orthogonalize<Speed, ThisFrame>(
      &reference_normal);
  Bivector<Product<Length, Speed>, ThisFrame> const reference_binormal =
      Wedge(reference_direction, reference_normal);
  *rotation = Rotation<InertialFrame, ThisFrame>(
                  R3x3Matrix(Normalize(reference_direction).coordinates(),
                             Normalize(reference_normal).coordinates(),
                             Normalize(reference_binormal).coordinates()));
  *angular_velocity =
      (Radian / Pow<2>(reference_direction.Norm())) * reference_binormal;
}

template<typename InertialFrame, typename ThisFrame>
void BarycentricRotatingDynamicFrame<InertialFrame, ThisFrame>::
FromBasisOfThisFrameToBasisOfInertialFrame(
    DegreesOfFreedom<InertialFrame> const& barycentre_degrees_of_freedom,
    DegreesOfFreedom<InertialFrame> const& primary_degrees_of_freedom,
    DegreesOfFreedom<InertialFrame> const& secondary_degrees_of_freedom,
    not_null<Rotation<InertialFrame, ThisFrame>*> const rotation,
    not_null<AngularVelocity<InertialFrame>*> const angular_velocity) {
  Rotation<InertialFrame, ThisFrame>
      from_basis_of_inertial_frame_to_basis_of_this_frame =
          Rotation<InertialFrame, ThisFrame>::Identity();
  AngularVelocity<InertialFrame> angular_velocity;
  FromBasisOfInertialFrameToBasisOfThisFrame(
      *last_barycentre_degrees_of_freedom,
      primary_degrees_of_freedom,
      secondary_degrees_of_freedom,
      &from_basis_of_inertial_frame_to_basis_of_this_frame,
      &angular_velocity);
  *rotation = from_basis_of_inertial_frame_to_basis_of_this_frame.Inverse();
}

}  // namespace physics
}  // namespace principia
