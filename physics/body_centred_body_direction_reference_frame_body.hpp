#pragma once

#include "physics/body_centred_body_direction_reference_frame.hpp"

#include <algorithm>
#include <utility>

#include "geometry/r3x3_matrix.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace _body_centred_body_direction_reference_frame {
namespace internal {

using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_orthogonal_map;
using namespace principia::geometry::_r3x3_matrix;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
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
  return std::max(primary_trajectory_().t_min(),
                  secondary_trajectory_->t_min());
}

template<typename InertialFrame, typename ThisFrame>
Instant BodyCentredBodyDirectionReferenceFrame<InertialFrame,
                                               ThisFrame>::t_max() const {
  return std::min(primary_trajectory_().t_max(),
                  secondary_trajectory_->t_max());
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

  Rotation<InertialFrame, ThisFrame> rotation =
      Rotation<InertialFrame, ThisFrame>::Identity();
  AngularVelocity<InertialFrame> angular_velocity;
  ComputeAngularDegreesOfFreedom(primary_degrees_of_freedom,
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

  // TODO(egg): eventually we want to add the intrinsic acceleration here.
  Vector<Acceleration, InertialFrame> const primary_acceleration =
      compute_gravitational_acceleration_on_primary_(
          primary_degrees_of_freedom.position(), t);

  Vector<Acceleration, InertialFrame> const secondary_acceleration =
      ephemeris_->ComputeGravitationalAccelerationOnMassiveBody(secondary_, t);

  auto const to_this_frame = ToThisFrameAtTime(t);

  // TODO(egg): TeX and reference.
  RelativeDegreesOfFreedom<InertialFrame> const secondary_primary =
      secondary_degrees_of_freedom - primary_degrees_of_freedom;
  Displacement<InertialFrame> const& r = secondary_primary.displacement();
  Velocity<InertialFrame> const& ṙ = secondary_primary.velocity();
  Vector<Acceleration, InertialFrame> const r̈ =
      secondary_acceleration - primary_acceleration;
  AngularVelocity<InertialFrame> const& ω =
      to_this_frame.template angular_velocity_of<ThisFrame>();
  Variation<AngularVelocity<InertialFrame>> const
      angular_acceleration_of_to_frame =
          (Wedge(r, r̈) * Radian - 2 * ω * InnerProduct(r, ṙ)) / r.Norm²();

  Vector<Acceleration, InertialFrame> const& acceleration_of_to_frame_origin =
      primary_acceleration;
  return AcceleratedRigidMotion<InertialFrame, ThisFrame>(
             to_this_frame,
             angular_acceleration_of_to_frame,
             acceleration_of_to_frame_origin);
}

}  // namespace internal
}  // namespace _body_centred_body_direction_reference_frame
}  // namespace physics
}  // namespace principia
