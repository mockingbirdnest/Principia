#pragma once

#include "physics/body_centred_body_direction_reference_frame.hpp"

#include <algorithm>
#include <memory>
#include <utility>

#include "geometry/orthogonal_map.hpp"
#include "geometry/r3x3_matrix.hpp"
#include "geometry/space_transformations.hpp"
#include "numerics/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace _body_centred_body_direction_reference_frame {
namespace internal {

using namespace principia::geometry::_orthogonal_map;
using namespace principia::geometry::_r3x3_matrix;
using namespace principia::geometry::_space_transformations;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

template<typename InertialFrame, typename ThisFrame>
BodyCentredBodyDirectionReferenceFrame<InertialFrame, ThisFrame>::
BodyCentredBodyDirectionReferenceFrame(
    not_null<Ephemeris<InertialFrame> const*> ephemeris,
    not_null<MassiveBody const*> const primary,
    not_null<MassiveBody const*> secondary)
    : ephemeris_(std::move(ephemeris)),
      primary_(primary),
      secondary_(std::move(secondary)),
      compute_gravitational_acceleration_on_primary_(
          [this](Position<InertialFrame> const& position, Instant const& t) {
            return ephemeris_->ComputeGravitationalAccelerationOnMassiveBody(
                primary_, t);
          }),
      compute_gravitational_jerk_on_primary_(
          [this](DegreesOfFreedom<InertialFrame> const& degrees_of_freedom,
                 Instant const& t) {
            return ephemeris_->ComputeGravitationalJerkOnMassiveBody(primary_,
                                                                     t);
          }),
      primary_trajectory_(
          [&t = *ephemeris_->trajectory(primary_)]() -> auto& { return t; }),
      secondary_trajectory_(ephemeris_->trajectory(secondary_)) {}

template<typename InertialFrame, typename ThisFrame>
BodyCentredBodyDirectionReferenceFrame<InertialFrame, ThisFrame>::
BodyCentredBodyDirectionReferenceFrame(
    not_null<Ephemeris<InertialFrame> const*> ephemeris,
    std::function<Trajectory<InertialFrame> const&()> primary_trajectory,
    not_null<MassiveBody const*> secondary)
    : ephemeris_(std::move(ephemeris)),
      primary_(nullptr),
      secondary_(std::move(secondary)),
      compute_gravitational_acceleration_on_primary_(
          [this](Position<InertialFrame> const& position, Instant const& t) {
            return ephemeris_->ComputeGravitationalAccelerationOnMasslessBody(
                position, t);
          }),
      compute_gravitational_jerk_on_primary_(
          [this](DegreesOfFreedom<InertialFrame> const& degrees_of_freedom,
                 Instant const& t) {
            return ephemeris_->ComputeGravitationalJerkOnMasslessBody(
                degrees_of_freedom, t);
          }),
      primary_trajectory_(std::move(primary_trajectory)),
      secondary_trajectory_(ephemeris_->trajectory(secondary_)) {}

template<typename InertialFrame, typename ThisFrame>
not_null<MassiveBody const*>
BodyCentredBodyDirectionReferenceFrame<InertialFrame, ThisFrame>::primary()
    const {
  return primary_;
}

template<typename InertialFrame, typename ThisFrame>
not_null<MassiveBody const*>
BodyCentredBodyDirectionReferenceFrame<InertialFrame, ThisFrame>::secondary()
    const {
  return secondary_;
}

template<typename InertialFrame, typename ThisFrame>
Instant BodyCentredBodyDirectionReferenceFrame<InertialFrame,
                                               ThisFrame>::t_min() const {
  // We depend on all bodies via the gravitational acceleration; in addition,
  // the primary may not be a massive body from the ephemeris.
  return std::max(primary_trajectory_().t_min(), ephemeris_->t_min());
}

template<typename InertialFrame, typename ThisFrame>
Instant BodyCentredBodyDirectionReferenceFrame<InertialFrame,
                                               ThisFrame>::t_max() const {
  // We depend on all bodies via the gravitational acceleration; in addition,
  // the primary may not be a massive body from the ephemeris.
  return std::min(primary_trajectory_().t_max(), ephemeris_->t_max());
}

template<typename InertialFrame, typename ThisFrame>
RigidMotion<InertialFrame, ThisFrame>
BodyCentredBodyDirectionReferenceFrame<InertialFrame, ThisFrame>::
ToThisFrameAtTime(Instant const& t) const {
  DegreesOfFreedom<InertialFrame> const primary_degrees_of_freedom =
      primary_trajectory_().EvaluateDegreesOfFreedom(t);
  DegreesOfFreedom<InertialFrame> const secondary_degrees_of_freedom =
      secondary_trajectory_->EvaluateDegreesOfFreedom(t);

  Vector<Acceleration, InertialFrame> const primary_acceleration =
      compute_gravitational_acceleration_on_primary_(
          primary_degrees_of_freedom.position(), t);
  Vector<Acceleration, InertialFrame> const secondary_acceleration =
      ephemeris_->ComputeGravitationalAccelerationOnMassiveBody(secondary_, t);

  return ToThisFrame(primary_degrees_of_freedom,
                     secondary_degrees_of_freedom,
                     primary_acceleration,
                     secondary_acceleration);
}

template<typename InertialFrame, typename ThisFrame>
void BodyCentredBodyDirectionReferenceFrame<InertialFrame, ThisFrame>::
WriteToMessage(not_null<serialization::ReferenceFrame*> const message) const {
  auto* const extension =
      message->MutableExtension(
          serialization::BodyCentredBodyDirectionReferenceFrame::extension);
  extension->set_primary(ephemeris_->serialization_index_for_body(primary_));
  extension->set_secondary(
      ephemeris_->serialization_index_for_body(secondary_));
}

template<typename InertialFrame, typename ThisFrame>
not_null<std::unique_ptr<
    BodyCentredBodyDirectionReferenceFrame<InertialFrame, ThisFrame>>>
BodyCentredBodyDirectionReferenceFrame<InertialFrame, ThisFrame>::
ReadFromMessage(
    not_null<Ephemeris<InertialFrame> const*> const ephemeris,
    serialization::BodyCentredBodyDirectionReferenceFrame const & message) {
  return std::make_unique<BodyCentredBodyDirectionReferenceFrame>(
      ephemeris,
      ephemeris->body_for_serialization_index(message.primary()),
      ephemeris->body_for_serialization_index(message.secondary()));
}

template<typename InertialFrame, typename ThisFrame>
Vector<Acceleration, InertialFrame>
BodyCentredBodyDirectionReferenceFrame<InertialFrame, ThisFrame>::
GravitationalAcceleration(Instant const& t,
                          Position<InertialFrame> const& q) const {
  return ephemeris_->ComputeGravitationalAccelerationOnMasslessBody(q, t);
}

template<typename InertialFrame, typename ThisFrame>
SpecificEnergy
BodyCentredBodyDirectionReferenceFrame<InertialFrame, ThisFrame>::
GravitationalPotential(Instant const& t,
                       Position<InertialFrame> const& q) const {
  return ephemeris_->ComputeGravitationalPotential(q, t);
}

template<typename InertialFrame, typename ThisFrame>
AcceleratedRigidMotion<InertialFrame, ThisFrame>
BodyCentredBodyDirectionReferenceFrame<InertialFrame, ThisFrame>::
MotionOfThisFrame(Instant const& t) const {
  DegreesOfFreedom<InertialFrame> const primary_degrees_of_freedom =
      primary_trajectory_().EvaluateDegreesOfFreedom(t);
  DegreesOfFreedom<InertialFrame> const secondary_degrees_of_freedom =
      secondary_trajectory_->EvaluateDegreesOfFreedom(t);

  Vector<Acceleration, InertialFrame> const primary_acceleration =
      compute_gravitational_acceleration_on_primary_(
          primary_degrees_of_freedom.position(), t);
  Vector<Acceleration, InertialFrame> const secondary_acceleration =
      ephemeris_->ComputeGravitationalAccelerationOnMassiveBody(secondary_, t);

  Vector<Jerk, InertialFrame> const primary_jerk =
      compute_gravitational_jerk_on_primary_(primary_degrees_of_freedom, t);
  Vector<Jerk, InertialFrame> const secondary_jerk =
      ephemeris_->ComputeGravitationalJerkOnMassiveBody(secondary_, t);

  auto const to_this_frame = ToThisFrame(primary_degrees_of_freedom,
                                         secondary_degrees_of_freedom,
                                         primary_acceleration,
                                         secondary_acceleration);

  Displacement<InertialFrame> const r =
      secondary_degrees_of_freedom.position() -
      primary_degrees_of_freedom.position();
  Velocity<InertialFrame> const ṙ =
      secondary_degrees_of_freedom.velocity() -
      primary_degrees_of_freedom.velocity();
  Vector<Acceleration, InertialFrame> const r̈ =
      secondary_acceleration - primary_acceleration;
  Vector<Jerk, InertialFrame> const r⁽³⁾ = secondary_jerk - primary_jerk;

  Trihedron<Length, ArealSpeed> orthogonal;
  Trihedron<double, double> orthonormal;
  Trihedron<Length, ArealSpeed, 1> 𝛛orthogonal;
  Trihedron<double, double, 1> 𝛛orthonormal;
  Trihedron<Length, ArealSpeed, 2> 𝛛²orthogonal;
  Trihedron<double, double, 2> 𝛛²orthonormal;

  Base::ComputeTrihedra(r, ṙ, orthogonal, orthonormal);
  Base::ComputeTrihedraDerivatives(r, ṙ, r̈,
                                   orthogonal, orthonormal,
                                   𝛛orthogonal, 𝛛orthonormal);
  Base::ComputeTrihedraDerivatives2(r, ṙ, r̈, r⁽³⁾,
                                    orthogonal, orthonormal,
                                    𝛛orthogonal, 𝛛orthonormal,
                                    𝛛²orthogonal, 𝛛²orthonormal);

  auto const angular_acceleration_of_to_frame =
      Base::ComputeAngularAcceleration(
          orthonormal, 𝛛orthonormal, 𝛛²orthonormal);

  Vector<Acceleration, InertialFrame> const& acceleration_of_to_frame_origin =
      primary_acceleration;
  return AcceleratedRigidMotion<InertialFrame, ThisFrame>(
             to_this_frame,
             angular_acceleration_of_to_frame,
             acceleration_of_to_frame_origin);}

template<typename InertialFrame, typename ThisFrame>
RigidMotion<InertialFrame, ThisFrame>
BodyCentredBodyDirectionReferenceFrame<InertialFrame, ThisFrame>::ToThisFrame(
    DegreesOfFreedom<InertialFrame> const& primary_degrees_of_freedom,
    DegreesOfFreedom<InertialFrame> const& secondary_degrees_of_freedom,
    Vector<Acceleration, InertialFrame> const& primary_acceleration,
    Vector<Acceleration, InertialFrame> const& secondary_acceleration) {
  Rotation<InertialFrame, ThisFrame> rotation =
      Rotation<InertialFrame, ThisFrame>::Identity();
  AngularVelocity<InertialFrame> angular_velocity;
  Base::ComputeAngularDegreesOfFreedom(primary_degrees_of_freedom,
                                       secondary_degrees_of_freedom,
                                       primary_acceleration,
                                       secondary_acceleration,
                                       rotation,
                                       angular_velocity);

  RigidTransformation<InertialFrame, ThisFrame> const
      rigid_transformation(primary_degrees_of_freedom.position(),
                           ThisFrame::origin,
                           rotation.template Forget<OrthogonalMap>());
  return RigidMotion<InertialFrame, ThisFrame>(
             rigid_transformation,
             angular_velocity,
             primary_degrees_of_freedom.velocity());
}

}  // namespace internal
}  // namespace _body_centred_body_direction_reference_frame
}  // namespace physics
}  // namespace principia
