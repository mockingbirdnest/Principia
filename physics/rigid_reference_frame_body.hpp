#pragma once

#include "physics/rigid_reference_frame.hpp"

#include "physics/barycentric_rotating_reference_frame.hpp"
#include "physics/body_centred_body_direction_reference_frame.hpp"
#include "physics/body_centred_non_rotating_reference_frame.hpp"
#include "physics/body_surface_reference_frame.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace _rigid_reference_frame {
namespace internal {

using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_r3x3_matrix;
using namespace principia::physics::_barycentric_rotating_reference_frame;
using namespace principia::physics::_body_centred_body_direction_reference_frame;  // NOLINT
using namespace principia::physics::_body_centred_non_rotating_reference_frame;
using namespace principia::physics::_body_surface_reference_frame;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

template<typename InertialFrame, typename ThisFrame>
SimilarMotion<InertialFrame, ThisFrame>
RigidReferenceFrame<InertialFrame, ThisFrame>::ToThisFrameAtTimeSimilarly(
    Instant const& t) const {
  return ToThisFrameAtTime(t).template Forget<SimilarMotion>();
}

template<typename InertialFrame, typename ThisFrame>
SimilarMotion<ThisFrame, InertialFrame>
RigidReferenceFrame<InertialFrame, ThisFrame>::FromThisFrameAtTimeSimilarly(
    Instant const& t) const {
  return FromThisFrameAtTime(t).template Forget<SimilarMotion>();
}

template<typename InertialFrame, typename ThisFrame>
RigidMotion<InertialFrame, ThisFrame>
RigidReferenceFrame<InertialFrame, ThisFrame>::ToThisFrameAtTime(
    Instant const& t) const {
  return FromThisFrameAtTime(t).Inverse();
}

template<typename InertialFrame, typename ThisFrame>
RigidMotion<ThisFrame, InertialFrame>
RigidReferenceFrame<InertialFrame, ThisFrame>::FromThisFrameAtTime(
    Instant const& t) const {
  return ToThisFrameAtTime(t).Inverse();
}

template<typename InertialFrame, typename ThisFrame>
Vector<Acceleration, ThisFrame>
RigidReferenceFrame<InertialFrame, ThisFrame>::GeometricAcceleration(
    Instant const& t,
    DegreesOfFreedom<ThisFrame> const& degrees_of_freedom) const {
  Vector<Acceleration, ThisFrame> gravitational_acceleration;
  Vector<Acceleration, ThisFrame> linear_acceleration;
  Vector<Acceleration, ThisFrame> coriolis_acceleration;
  Vector<Acceleration, ThisFrame> centrifugal_acceleration;
  Vector<Acceleration, ThisFrame> euler_acceleration;
  ComputeGeometricAccelerations(t,
                                degrees_of_freedom,
                                gravitational_acceleration,
                                linear_acceleration,
                                coriolis_acceleration,
                                centrifugal_acceleration,
                                euler_acceleration);

  return gravitational_acceleration +
         (linear_acceleration + coriolis_acceleration +
          centrifugal_acceleration + euler_acceleration);
}

template<typename InertialFrame, typename ThisFrame>
Vector<Acceleration, ThisFrame>
RigidReferenceFrame<InertialFrame, ThisFrame>::
RotationFreeGeometricAccelerationAtRest(
    Instant const& t,
    Position<ThisFrame> const& position) const {
  Vector<Acceleration, ThisFrame> gravitational_acceleration;
  Vector<Acceleration, ThisFrame> linear_acceleration;
  Vector<Acceleration, ThisFrame> coriolis_acceleration;
  Vector<Acceleration, ThisFrame> centrifugal_acceleration;
  Vector<Acceleration, ThisFrame> euler_acceleration;
  ComputeGeometricAccelerations(t,
                                {position, ThisFrame::unmoving},
                                gravitational_acceleration,
                                linear_acceleration,
                                coriolis_acceleration,
                                centrifugal_acceleration,
                                euler_acceleration);

  DCHECK_EQ(coriolis_acceleration, (Vector<Acceleration, ThisFrame>{}));
  return gravitational_acceleration +
         (linear_acceleration + centrifugal_acceleration);
}

template<typename InertialFrame, typename ThisFrame>
SpecificEnergy
RigidReferenceFrame<InertialFrame, ThisFrame>::GeometricPotential(
    Instant const& t,
    Position<ThisFrame> const& position) const {
  AcceleratedRigidMotion<InertialFrame, ThisFrame> const motion =
      MotionOfThisFrame(t);
  RigidMotion<InertialFrame, ThisFrame> const& to_this_frame =
      motion.rigid_motion();
  RigidMotion<ThisFrame, InertialFrame> const from_this_frame =
      to_this_frame.Inverse();

  AngularVelocity<ThisFrame> const Ω = to_this_frame.orthogonal_map()(
      to_this_frame.template angular_velocity_of<ThisFrame>());
  Displacement<ThisFrame> const r = position - ThisFrame::origin;

  SpecificEnergy const gravitational_potential =
      GravitationalPotential(t,
                             from_this_frame.rigid_transformation()(position));
  SpecificEnergy const linear_potential = InnerProduct(
      r,
      to_this_frame.orthogonal_map()(
          motion.template acceleration_of_origin_of<ThisFrame>()));
  SpecificEnergy const centrifugal_potential = -0.5 * (Ω * r / Radian).Norm²();

  return gravitational_potential + (linear_potential + centrifugal_potential);
}

template<typename InertialFrame, typename ThisFrame>
not_null<std::unique_ptr<RigidReferenceFrame<InertialFrame, ThisFrame>>>
RigidReferenceFrame<InertialFrame, ThisFrame>::ReadFromMessage(
    serialization::ReferenceFrame const& message,
    not_null<Ephemeris<InertialFrame> const*> const ephemeris) {
  std::unique_ptr<RigidReferenceFrame> result;
  int extensions_found = 0;
  // NOTE(egg): the |static_cast|ing below is needed on MSVC, because the silly
  // compiler doesn't see the
  // |operator std::unique_ptr<RigidReferenceFrame>() &&|.
  if (message.HasExtension(
          serialization::BarycentricRotatingReferenceFrame::extension)) {
    ++extensions_found;
    result = static_cast<not_null<std::unique_ptr<RigidReferenceFrame>>>(
        BarycentricRotatingReferenceFrame<InertialFrame, ThisFrame>::
            ReadFromMessage(
                ephemeris,
                message.GetExtension(
                    serialization::BarycentricRotatingReferenceFrame::
                        extension)));
  }
  if (message.HasExtension(
          serialization::BodyCentredBodyDirectionReferenceFrame::extension)) {
    ++extensions_found;
    result = static_cast<not_null<std::unique_ptr<RigidReferenceFrame>>>(
        BodyCentredBodyDirectionReferenceFrame<InertialFrame, ThisFrame>::
            ReadFromMessage(
                ephemeris,
                message.GetExtension(
                    serialization::BodyCentredBodyDirectionReferenceFrame::
                        extension)));
  }
  if (message.HasExtension(
          serialization::BodyCentredNonRotatingReferenceFrame::extension)) {
    ++extensions_found;
    result = static_cast<not_null<std::unique_ptr<RigidReferenceFrame>>>(
        BodyCentredNonRotatingReferenceFrame<InertialFrame, ThisFrame>::
            ReadFromMessage(
                ephemeris,
                message.GetExtension(
                    serialization::BodyCentredNonRotatingReferenceFrame::
                        extension)));
  }
  if (message.HasExtension(
          serialization::BodySurfaceReferenceFrame::extension)) {
    ++extensions_found;
    result = static_cast<not_null<std::unique_ptr<RigidReferenceFrame>>>(
        BodySurfaceReferenceFrame<InertialFrame, ThisFrame>::
            ReadFromMessage(
                ephemeris,
                message.GetExtension(
                    serialization::BodySurfaceReferenceFrame::extension)));
  }
  CHECK_EQ(extensions_found, 1) << message.DebugString();
  return std::move(result);
}

template<typename InertialFrame, typename ThisFrame>
void RigidReferenceFrame<InertialFrame, ThisFrame>::
ComputeAngularDegreesOfFreedom(
    DegreesOfFreedom<InertialFrame> const& primary_degrees_of_freedom,
    DegreesOfFreedom<InertialFrame> const& secondary_degrees_of_freedom,
    Vector<Acceleration, InertialFrame> const& primary_acceleration,
    Vector<Acceleration, InertialFrame> const& secondary_acceleration,
    Rotation<InertialFrame, ThisFrame>& rotation,
    AngularVelocity<InertialFrame>& angular_velocity) {
  RelativeDegreesOfFreedom<InertialFrame> const reference =
       secondary_degrees_of_freedom - primary_degrees_of_freedom;

  Displacement<InertialFrame> const& r = reference.displacement();
  Velocity<InertialFrame> const ṙ = reference.velocity();
  Vector<Acceleration, InertialFrame> const r̈ =
      secondary_acceleration - primary_acceleration;

  Trihedron<Length, Speed> orthogonal;
  Trihedron<double, double> orthonormal;
  Trihedron<Length, Speed, 1> 𝛛orthogonal;
  Trihedron<double, double, 1> 𝛛orthonormal;

  ComputeTrihedra(r, ṙ, orthogonal, orthonormal);
  ComputeTrihedraDerivatives(
      r, ṙ, r̈, orthogonal, orthonormal, 𝛛orthogonal, 𝛛orthonormal);
  rotation = ComputeRotation(orthonormal);
  angular_velocity = ComputeAngularVelocity(orthonormal, 𝛛orthonormal);
}

template<typename InertialFrame, typename ThisFrame>
void RigidReferenceFrame<InertialFrame, ThisFrame>::ComputeTrihedra(
    Displacement<InertialFrame> const& r,
    Velocity<InertialFrame> const& ṙ,
    Trihedron<Length, Speed>& orthogonal,
    Trihedron<double, double>& orthonormal) {
  // Our orthogonal (but not orthonormal) trihedron for |ThisFrame|.
  Displacement<InertialFrame> const& T = r;
  Velocity<InertialFrame> const N = ṙ.OrthogonalizationAgainst(r);
  Bivector<Product<Length, Speed>, InertialFrame> const B = Wedge(T, N);

  // Our orthonormal trihedron.
  Vector<double, InertialFrame> const t = Normalize(T);
  Vector<double, InertialFrame> const n = Normalize(N);
  Bivector<double, InertialFrame> const b = Normalize(B);

  orthogonal = {.tangent = T, .normal = N, .binormal = B};
  orthonormal = {.tangent = t, .normal = n, .binormal = b};
}

template<typename InertialFrame, typename ThisFrame>
void RigidReferenceFrame<InertialFrame, ThisFrame>::ComputeTrihedraDerivatives(
    Displacement<InertialFrame> const& r,
    Velocity<InertialFrame> const& ṙ,
    Vector<Acceleration, InertialFrame> const& r̈,
    Trihedron<Length, Speed>& orthogonal,
    Trihedron<double, double>& orthonormal,
    Trihedron<Length, Speed, 1>& 𝛛orthogonal,
    Trihedron<double, double, 1>& 𝛛orthonormal) {
  auto const& T = orthogonal.tangent;
  auto const& N = orthogonal.normal;
  auto const& B = orthogonal.binormal;
  auto const& t = orthonormal.tangent;
  auto const& n = orthonormal.normal;
  auto const& b = orthonormal.binormal;

  // The derivatives of the |orthogonal| trihedron.
  Velocity<InertialFrame> const& Ṫ = ṙ;
  Vector<Acceleration, InertialFrame> const Ṅ =
      r̈ + 2 * Pow<2>(InnerProduct(r, ṙ) / r.Norm²()) * r -
      (r * (ṙ.Norm²() + InnerProduct(r, r̈)) + ṙ * InnerProduct(r, ṙ)) /
          r.Norm²();
  Bivector<Variation<Product<Length, Speed>>, InertialFrame> const Ḃ =
      Wedge(r, r̈);

  // For any multivector v this returns the derivative of the normalized v.
  auto derive_normalized = []<typename V>(V const& v, Variation<V> const& v̇) {
    return (v.Norm²() * v̇ - InnerProduct(v, v̇) * v) / Pow<3>(v.Norm());
  };

  // The derivatives of the |orthonormal| trihedron.
  Vector<Variation<double>, InertialFrame> const ṫ = derive_normalized(T, Ṫ);
  Vector<Variation<double>, InertialFrame> const ṅ = derive_normalized(N, Ṅ);
  Bivector<Variation<double>, InertialFrame> const ḃ = derive_normalized(B, Ḃ);

  𝛛orthogonal = {.tangent = Ṫ, .normal = Ṅ, .binormal = Ḃ};
  𝛛orthonormal = {.tangent = ṫ, .normal = ṅ, .binormal = ḃ};
}

template<typename InertialFrame, typename ThisFrame>
Rotation<InertialFrame, ThisFrame>
RigidReferenceFrame<InertialFrame, ThisFrame>::ComputeRotation(
    Trihedron<double, double> const& orthonormal) {
  return Rotation<InertialFrame, ThisFrame>(orthonormal.tangent,
                                            orthonormal.normal,
                                            orthonormal.binormal);
}

template<typename InertialFrame, typename ThisFrame>
AngularVelocity<InertialFrame>
RigidReferenceFrame<InertialFrame, ThisFrame>::ComputeAngularVelocity(
    Trihedron<double, double> const& orthonormal,
    Trihedron<double, double, 1>& 𝛛orthonormal) {
  auto const& t = orthonormal.tangent;
  auto const& n = orthonormal.normal;
  auto const& b = orthonormal.binormal;
  auto const& ṫ = 𝛛orthonormal.tangent;
  auto const& ṅ = 𝛛orthonormal.normal;
  auto const& ḃ = 𝛛orthonormal.binormal;

  return Radian * (Wedge(ṅ, b) * t + Wedge(ḃ, t) * n + InnerProduct(ṫ, n) * b);
}

template<typename InertialFrame, typename ThisFrame>
void RigidReferenceFrame<InertialFrame, ThisFrame>::
ComputeGeometricAccelerations(
    Instant const& t,
    DegreesOfFreedom<ThisFrame> const& degrees_of_freedom,
    Vector<Acceleration, ThisFrame>& gravitational_acceleration,
    Vector<Acceleration, ThisFrame>& linear_acceleration,
    Vector<Acceleration, ThisFrame>& coriolis_acceleration,
    Vector<Acceleration, ThisFrame>& centrifugal_acceleration,
    Vector<Acceleration, ThisFrame>& euler_acceleration) const {
  AcceleratedRigidMotion<InertialFrame, ThisFrame> const motion =
      MotionOfThisFrame(t);
  RigidMotion<InertialFrame, ThisFrame> const& to_this_frame =
      motion.rigid_motion();
  RigidMotion<ThisFrame, InertialFrame> const from_this_frame =
      to_this_frame.Inverse();

  // Beware, we want the angular velocity of ThisFrame as seen in the
  // InertialFrame, but pushed to ThisFrame.  Otherwise the sign is wrong.
  AngularVelocity<ThisFrame> const Ω = to_this_frame.orthogonal_map()(
      to_this_frame.template angular_velocity_of<ThisFrame>());
  Variation<AngularVelocity<ThisFrame>> const dΩ_over_dt =
      to_this_frame.orthogonal_map()(
          motion.template angular_acceleration_of<ThisFrame>());
  Displacement<ThisFrame> const r =
      degrees_of_freedom.position() - ThisFrame::origin;

  gravitational_acceleration = to_this_frame.orthogonal_map()(
      GravitationalAcceleration(t,
                                from_this_frame.rigid_transformation()(
                                    degrees_of_freedom.position())));
  linear_acceleration = -to_this_frame.orthogonal_map()(
      motion.template acceleration_of_origin_of<ThisFrame>());
  coriolis_acceleration = -2 * Ω * degrees_of_freedom.velocity() / Radian;
  centrifugal_acceleration = -Ω * (Ω * r) / Pow<2>(Radian);
  euler_acceleration = -dΩ_over_dt * r / Radian;
}

}  // namespace internal
}  // namespace _rigid_reference_frame
}  // namespace physics
}  // namespace principia
