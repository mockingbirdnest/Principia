#pragma once

#include "physics/rotating_body.hpp"

#include <algorithm>
#include <vector>
#include "quantities/constants.hpp"

namespace principia {

using constants::GravitationalConstant;

namespace physics {

namespace {
double const kNormLow = 0.999;
double const kNormHigh = 1.001;
}  // namespace

template<typename Frame>
RotatingBody<Frame>::Parameters::Parameters(
    Angle const& reference_angle,
    Instant const& reference_time,
    AngularVelocity<Frame> const& angular_velocity)
    : reference_angle_(reference_angle),
      reference_time_(reference_time),
      angular_velocity_(angular_velocity) {
  CHECK_NE(angular_velocity_.Norm(), 0.0)
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
bool RotatingBody<Frame>::is_massless() const {
  return false;
}

template<typename Frame>
bool RotatingBody<Frame>::is_oblate() const {
  return false;
}

template<typename Frame>
inline void RotatingBody<Frame>::WriteToMessage(
    not_null<serialization::Body*> const message) const {
  WriteToMessage(message->mutable_massive_body());
}

template<typename Frame>
inline void RotatingBody<Frame>::WriteToMessage(
    not_null<serialization::MassiveBody*> const message) const {
  MassiveBody::WriteToMessage(message);
  not_null<serialization::RotatingBody*> const rotating_body =
      message->MutableExtension(serialization::RotatingBody::rotating_body);
  Frame::WriteToMessage(rotating_body->mutable_frame());
  parameters_.j2_->WriteToMessage(rotating_body->mutable_j2());
  parameters_.axis_.WriteToMessage(rotating_body->mutable_axis());
}


template<typename Frame>
not_null<std::unique_ptr<RotatingBody<Frame>>> RotatingBody<Frame>::ReadFromMessage(
    serialization::Body const& message) {
  CHECK(message.has_massive_body());
  return ReadFromMessage(message.massive_body());
}

template<typename Frame>
not_null<std::unique_ptr<RotatingBody<Frame>>> RotatingBody<Frame>::ReadFromMessage(
    serialization::MassiveBody const& message) {
  CHECK(message.HasExtension(serialization::RotatingBody::rotating_body));
  serialization::RotatingBody const& oblateness_information =
      message.GetExtension(serialization::RotatingBody::rotating_body);
  return std::make_unique<RotatingBody<Frame>>(
      GravitationalParameter::ReadFromMessage(
          message.gravitational_parameter()),
      RotatingBody<Frame>::Parameters(
          Order2ZonalCoefficient::ReadFromMessage(
              oblateness_information.j2()),
          Vector<double, Frame>::ReadFromMessage(
              oblateness_information.axis())));
}

}  // namespace physics
}  // namespace principia
