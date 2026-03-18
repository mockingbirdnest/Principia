#pragma once

#include "physics/rigid_reference_frame.hpp"

#include <memory>
#include <utility>

#include "geometry/r3x3_matrix.hpp"
#include "numerics/elementary_functions.hpp"
#include "physics/barycentric_rotating_reference_frame.hpp"
#include "physics/body_centred_body_direction_reference_frame.hpp"
#include "physics/body_centred_non_rotating_reference_frame.hpp"
#include "physics/body_surface_reference_frame.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace _rigid_reference_frame {
namespace internal {

using namespace principia::geometry::_r3x3_matrix;
using namespace principia::numerics::_elementary_functions;
using namespace principia::physics::_barycentric_rotating_reference_frame;
using namespace principia::physics::_body_centred_body_direction_reference_frame;  // NOLINT
using namespace principia::physics::_body_centred_non_rotating_reference_frame;
using namespace principia::physics::_body_surface_reference_frame;
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
  // NOTE(egg): the `static_cast`ing below is needed on MSVC, because the silly
  // compiler doesn't see the
  // `operator std::unique_ptr<RigidReferenceFrame>() &&`.
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

  Trihedron<Length, ArealSpeed> orthogonal;
  Trihedron<double, double> orthonormal;
  Trihedron<Length, ArealSpeed, 1> 𝛛orthogonal;
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
    Trihedron<Length, ArealSpeed>& orthogonal,
    Trihedron<double, double>& orthonormal) {
  // Our orthogonal (but not orthonormal) trihedron for `ThisFrame`.
  Displacement<InertialFrame> const& F = r;
  Bivector<ArealSpeed, InertialFrame> const B = Wedge(r, ṙ);
  Vector<Product<Length, ArealSpeed>, InertialFrame> const N = B * F;

  // Our orthonormal trihedron.
  Vector<double, InertialFrame> const f = Normalize(F);
  Vector<double, InertialFrame> const n = Normalize(N);
  Bivector<double, InertialFrame> const b = Normalize(B);

  orthogonal = {.fore = F, .normal = N, .binormal = B};
  orthonormal = {.fore = f, .normal = n, .binormal = b};
}

template<typename InertialFrame, typename ThisFrame>
void RigidReferenceFrame<InertialFrame, ThisFrame>::ComputeTrihedraDerivatives(
    Displacement<InertialFrame> const& r,
    Velocity<InertialFrame> const& ṙ,
    Vector<Acceleration, InertialFrame> const& r̈,
    Trihedron<Length, ArealSpeed> const& orthogonal,
    Trihedron<double, double> const& /*orthonormal*/,
    Trihedron<Length, ArealSpeed, 1>& 𝛛orthogonal,
    Trihedron<double, double, 1>& 𝛛orthonormal) {
  auto const& F = orthogonal.fore;
  auto const& N = orthogonal.normal;
  auto const& B = orthogonal.binormal;

  // The derivatives of the `orthogonal` trihedron.
  Velocity<InertialFrame> const& Ḟ = ṙ;
  Bivector<Variation<ArealSpeed>, InertialFrame> const Ḃ = Wedge(r, r̈);
  Vector<Variation<Product<Length, ArealSpeed>>, InertialFrame> const Ṅ =
      Ḃ * F + B * Ḟ;

  // For any multivector v this returns the derivative of v / ‖v‖.
  auto 𝛛normalized = []<typename V>(V const& v, Variation<V> const& v̇) {
    return (v.Norm²() * v̇ - InnerProduct(v, v̇) * v) / Pow<3>(v.Norm());
  };

  // The derivatives of the `orthonormal` trihedron.
  Vector<Variation<double>, InertialFrame> const ḟ = 𝛛normalized(F, Ḟ);
  Vector<Variation<double>, InertialFrame> const ṅ = 𝛛normalized(N, Ṅ);
  Bivector<Variation<double>, InertialFrame> const ḃ = 𝛛normalized(B, Ḃ);

  𝛛orthogonal = {.fore = Ḟ, .normal = Ṅ, .binormal = Ḃ};
  𝛛orthonormal = {.fore = ḟ, .normal = ṅ, .binormal = ḃ};
}

template<typename InertialFrame, typename ThisFrame>
void RigidReferenceFrame<InertialFrame, ThisFrame>::ComputeTrihedraDerivatives2(
    Displacement<InertialFrame> const& r,
    Velocity<InertialFrame> const& ṙ,
    Vector<Acceleration, InertialFrame> const& r̈,
    Vector<Jerk, InertialFrame> const& r⁽³⁾,
    Trihedron<Length, ArealSpeed> const& orthogonal,
    Trihedron<double, double> const& /*orthonormal*/,
    Trihedron<Length, ArealSpeed, 1> const& 𝛛orthogonal,
    Trihedron<double, double, 1> const& /*𝛛orthonormal*/,
    Trihedron<Length, ArealSpeed, 2>& 𝛛²orthogonal,
    Trihedron<double, double, 2>& 𝛛²orthonormal) {
  auto const& F = orthogonal.fore;
  auto const& N = orthogonal.normal;
  auto const& B = orthogonal.binormal;
  auto const& Ḟ = 𝛛orthogonal.fore;
  auto const& Ṅ = 𝛛orthogonal.normal;
  auto const& Ḃ = 𝛛orthogonal.binormal;

  // The second derivatives of the `orthogonal` trihedron.
  Vector<Acceleration, InertialFrame> const& F̈ = r̈;
  Bivector<Variation<ArealSpeed, 2>, InertialFrame> const B̈ =
      Wedge(ṙ, r̈) + Wedge(r, r⁽³⁾);
  Vector<Variation<Product<Length, ArealSpeed>, 2>, InertialFrame> const N̈ =
      B̈ * F + 2 * Ḃ * Ḟ + B * F̈;

  // For any multivector v this returns the second derivative of v / ‖v‖.
  auto 𝛛²normalized = []<typename V>(V const& v,
                                     Variation<V> const& v̇,
                                     Variation<V, 2> const& v̈) {
    return v̈ / v.Norm() -
           (2 * InnerProduct(v, v̇) * v̇ + (v̇.Norm²() - InnerProduct(v, v̈)) * v) /
               Pow<3>(v.Norm()) +
           3 * v * Pow<2>(InnerProduct(v, v̇)) / Pow<5>(v.Norm());
  };

  // The second derivatives of the `orthonormal` trihedron.
  Vector<Variation<double, 2>, InertialFrame> const f̈ = 𝛛²normalized(F, Ḟ, F̈);
  Vector<Variation<double, 2>, InertialFrame> const n̈ = 𝛛²normalized(N, Ṅ, N̈);
  Bivector<Variation<double, 2>, InertialFrame> const b̈ = 𝛛²normalized(B, Ḃ, B̈);

  𝛛²orthogonal = {.fore = F̈, .normal = N̈, .binormal = B̈};
  𝛛²orthonormal = {.fore = f̈, .normal = n̈, .binormal = b̈};
}

template<typename InertialFrame, typename ThisFrame>
Rotation<InertialFrame, ThisFrame>
RigidReferenceFrame<InertialFrame, ThisFrame>::ComputeRotation(
    Trihedron<double, double> const& orthonormal) {
  return Rotation<InertialFrame, ThisFrame>(orthonormal.fore,
                                            orthonormal.normal,
                                            orthonormal.binormal);
}

template<typename InertialFrame, typename ThisFrame>
AngularVelocity<InertialFrame>
RigidReferenceFrame<InertialFrame, ThisFrame>::ComputeAngularVelocity(
    Trihedron<double, double> const& orthonormal,
    Trihedron<double, double, 1> const& 𝛛orthonormal) {
  auto const& f = orthonormal.fore;
  auto const& n = orthonormal.normal;
  auto const& b = orthonormal.binormal;
  auto const& ḟ = 𝛛orthonormal.fore;
  auto const& ṅ = 𝛛orthonormal.normal;
  auto const& ḃ = 𝛛orthonormal.binormal;

  return Radian * (Wedge(ṅ, b) * f + Wedge(ḃ, f) * n + InnerProduct(ḟ, n) * b);
}

template<typename InertialFrame, typename ThisFrame>
Bivector<AngularAcceleration, InertialFrame>
RigidReferenceFrame<InertialFrame, ThisFrame>::ComputeAngularAcceleration(
    Trihedron<double, double> const& orthonormal,
    Trihedron<double, double, 1> const& 𝛛orthonormal,
    Trihedron<double, double, 2> const& 𝛛²orthonormal) {
  auto const& f = orthonormal.fore;
  auto const& n = orthonormal.normal;
  auto const& b = orthonormal.binormal;
  auto const& ḟ = 𝛛orthonormal.fore;
  auto const& ṅ = 𝛛orthonormal.normal;
  auto const& ḃ = 𝛛orthonormal.binormal;
  auto const& f̈ = 𝛛²orthonormal.fore;
  auto const& n̈ = 𝛛²orthonormal.normal;
  auto const& b̈ = 𝛛²orthonormal.binormal;

  return Radian * (
      (Wedge(n̈, b) + Wedge(ṅ, ḃ)) * f + Wedge(ṅ, b) * ḟ +
      (Wedge(b̈, f) + Wedge(ḃ, ḟ)) * n + Wedge(ḃ, f) * ṅ +
      (InnerProduct(f̈, n) + InnerProduct(ḟ, ṅ)) * b + InnerProduct(ḟ, n) * ḃ);
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
