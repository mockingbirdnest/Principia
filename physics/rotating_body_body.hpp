
#pragma once

#include "physics/rotating_body.hpp"

#include <algorithm>
#include <vector>

#include "physics/oblate_body.hpp"
#include "quantities/constants.hpp"

namespace principia {

using geometry::Exp;

namespace physics {

template<typename Frame>
RotatingBody<Frame>::Parameters::Parameters(
    Length const& mean_radius,
    Angle const& reference_angle,
    Instant const& reference_instant,
    AngularVelocity<Frame> const& angular_velocity)
    : mean_radius_(mean_radius),
      reference_angle_(reference_angle),
      reference_instant_(reference_instant),
      angular_velocity_(angular_velocity) {
  CHECK_NE(angular_velocity_.Norm(), 0.0 * SIUnit<AngularFrequency>())
      << "Rotating body cannot have zero angular velocity";
}

template<typename Frame>
RotatingBody<Frame>::RotatingBody(
    MassiveBody::Parameters const& massive_body_parameters,
    Parameters const& parameters)
    : MassiveBody(massive_body_parameters),
      parameters_(parameters) {}

template<typename Frame>
Length RotatingBody<Frame>::mean_radius() const {
  return parameters_.mean_radius_;
}

template<typename Frame>
AngularVelocity<Frame> const& RotatingBody<Frame>::angular_velocity() const {
  return parameters_.angular_velocity_;
}

template<typename Frame>
Angle RotatingBody<Frame>::AngleAt(Instant const& t) const {
  return parameters_.reference_angle_ +
         (t - parameters_.reference_instant_) *
             parameters_.angular_velocity_.Norm();
}

template<typename Frame>
Rotation<Frame, Frame> RotatingBody<Frame>::RotationAt(Instant const& t) const {
  return Exp((t - parameters_.reference_instant_) *
                 parameters_.angular_velocity_);
}

template<typename Frame>
bool RotatingBody<Frame>::is_massless() const {
  return false;
}

template<typename Frame>
bool RotatingBody<Frame>::is_oblate() const {
  return false;
}

template<typename Frame>
void RotatingBody<Frame>::WriteToMessage(
    not_null<serialization::Body*> const message) const {
  WriteToMessage(message->mutable_massive_body());
}

template<typename Frame>
void RotatingBody<Frame>::WriteToMessage(
    not_null<serialization::MassiveBody*> const message) const {
  MassiveBody::WriteToMessage(message);
  not_null<serialization::RotatingBody*> const rotating_body =
      message->MutableExtension(serialization::RotatingBody::extension);
  Frame::WriteToMessage(rotating_body->mutable_frame());
  parameters_.reference_angle_.WriteToMessage(
      rotating_body->mutable_reference_angle());
  parameters_.reference_instant_.WriteToMessage(
      rotating_body->mutable_reference_instant());
  parameters_.angular_velocity_.WriteToMessage(
      rotating_body->mutable_angular_velocity());
}

template<typename Frame>
not_null<std::unique_ptr<RotatingBody<Frame>>>
RotatingBody<Frame>::ReadFromMessage(
    serialization::RotatingBody const& message,
    MassiveBody::Parameters const& massive_body_parameters) {
  // For pre-Buffon compatibility, build a point mass.
  Parameters parameters(
                 message.has_mean_radius()
                     ? Length::ReadFromMessage(message.mean_radius())
                     : Length(),
                 Angle::ReadFromMessage(message.reference_angle()),
                 Instant::ReadFromMessage(message.reference_instant()),
                 AngularVelocity<Frame>::ReadFromMessage(
                     message.angular_velocity()));

  if (message.HasExtension(serialization::OblateBody::extension)) {
    serialization::OblateBody const& extension =
        message.GetExtension(serialization::OblateBody::extension);

    return OblateBody<Frame>::ReadFromMessage(extension,
                                              massive_body_parameters,
                                              parameters);
  } else {
    return std::make_unique<RotatingBody<Frame>>(massive_body_parameters,
                                                 parameters);
  }
}

}  // namespace physics
}  // namespace principia
