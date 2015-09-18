
#pragma once

#include "physics/rotating_body.hpp"

#include <algorithm>
#include <vector>

#include "physics/oblate_body.hpp"
#include "quantities/constants.hpp"

namespace principia {
namespace physics {

template<typename Frame>
RotatingBody<Frame>::Parameters::Parameters(
    Angle const& reference_angle,
    Instant const& reference_instant,
    AngularVelocity<Frame> const& angular_velocity)
    : reference_angle_(reference_angle),
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
      message->MutableExtension(serialization::RotatingBody::rotating_body);
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
  Parameters parameters(
                 Angle::ReadFromMessage(message.reference_angle()),
                 Instant::ReadFromMessage(message.reference_instant()),
                 AngularVelocity<Frame>::ReadFromMessage(
                     message.angular_velocity()));

  if (message.HasExtension(serialization::OblateBody::oblate_body)) {
    serialization::OblateBody const& extension =
        message.GetExtension(serialization::OblateBody::oblate_body);

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
