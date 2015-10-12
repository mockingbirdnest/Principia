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
      from_basis_of_barycentric_frame_to_standard_basis =
          Rotation<InertialFrame, ThisFrame>::Identity();
  Bivector<AngularFrequency, InertialFrame> angular_frequency;
  FromBasisOfInertialFrameToBasisOfThisFrame(
      barycentre_degrees_of_freedom,
      primary_degrees_of_freedom,
      secondary_degrees_of_freedom,
      &from_basis_of_barycentric_frame_to_standard_basis,
      &angular_frequency);

}

template<typename InertialFrame, typename ThisFrame>
RigidMotion<ThisFrame, InertialFrame>
BarycentricRotatingDynamicFrame<InertialFrame, ThisFrame>::
FromThisFrameAtTime(Instant const& t) const {
  DegreesOfFreedom<InertialFrame> const centre_degrees_of_freedom =
      centre_trajectory_->EvaluateDegreesOfFreedom(t, &hint_);
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
    not_null<Bivector<AngularFrequency, InertialFrame>*> const
        angular_frequency) {
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
  *angular_frequency =
      (Radian / Pow<2>(reference_direction.Norm())) * reference_binormal;
}

}  // namespace physics
}  // namespace principia
