#pragma once

#include "physics/oblate_body.hpp"

#include <algorithm>
#include <vector>

#include "quantities/constants.hpp"

namespace principia {

using constants::GravitationalConstant;

namespace physics {

template<typename Frame>
OblateBody<Frame>::Parameters::Parameters(double const j2,
                                          Length const& radius)
    : j2_over_μ_(-j2 * radius * radius) {
  CHECK_NE(j2, 0.0) << "Oblate body cannot have zero j2";
}

template<typename Frame>
OblateBody<Frame>::Parameters::Parameters(Order2ZonalCoefficient const& j2)
    : j2_(j2) {
  CHECK_NE(j2, Order2ZonalCoefficient()) << "Oblate body cannot have zero j2";
}

template<typename Frame>
OblateBody<Frame>::OblateBody(
    MassiveBody::Parameters const& massive_body_parameters,
    typename RotatingBody<Frame>::Parameters const& rotating_body_parameters,
    Parameters const& parameters)
    : RotatingBody(massive_body_parameters, rotating_body_parameters),
      parameters_(parameters),
      axis_(
          Vector<double, Frame>(Normalize(angular_velocity()).coordinates())) {
  if (parameters_.j2_) {
    parameters_.j2_over_μ_ = *parameters_.j2_ / gravitational_parameter();
  }
  if (parameters_.j2_over_μ_) {
    parameters_.j2_ = *parameters_.j2_over_μ_ * gravitational_parameter();
  }
}

template<typename Frame>
Order2ZonalCoefficient const& OblateBody<Frame>::j2() const {
  return *parameters_.j2_;
}

template<typename Frame>
Quotient<Order2ZonalCoefficient,
         GravitationalParameter> const& OblateBody<Frame>::j2_over_μ() const {
  return *parameters_.j2_over_μ_;
}

template<typename Frame>
Vector<double, Frame> const& OblateBody<Frame>::axis() const {
  return axis_;
}

template<typename Frame>
bool OblateBody<Frame>::is_massless() const {
  return false;
}

template<typename Frame>
bool OblateBody<Frame>::is_oblate() const {
  return true;
}

template<typename Frame>
void OblateBody<Frame>::WriteToMessage(
    not_null<serialization::Body*> const message) const {
  WriteToMessage(message->mutable_massive_body());
}

template<typename Frame>
void OblateBody<Frame>::WriteToMessage(
    not_null<serialization::MassiveBody*> const message) const {
  RotatingBody<Frame>::WriteToMessage(message);
  not_null<serialization::OblateBody*> const oblate_body =
      message->MutableExtension(serialization::RotatingBody::rotating_body)->
               MutableExtension(serialization::OblateBody::oblate_body);
  parameters_.j2_->WriteToMessage(oblate_body->mutable_j2());
}

template<typename Frame>
not_null<std::unique_ptr<OblateBody<Frame>>> OblateBody<Frame>::ReadFromMessage(
      serialization::OblateBody const& message,
      MassiveBody::Parameters const& massive_body_parameters,
      typename RotatingBody<Frame>::Parameters const&
          rotating_body_parameters) {
  Parameters parameters(Order2ZonalCoefficient::ReadFromMessage(message.j2()));

  return std::make_unique<OblateBody<Frame>>(massive_body_parameters,
                                             rotating_body_parameters,
                                             parameters);
}

template<typename Frame>
not_null<std::unique_ptr<OblateBody<Frame>>> OblateBody<Frame>::ReadFromMessage(
    serialization::PreBrouwerOblateBody const& message,
    MassiveBody::Parameters const& massive_body_parameters) {
  auto const axis = Vector<double, Frame>::ReadFromMessage(message.axis());
  typename RotatingBody<Frame>::Parameters
      rotating_body_parameters(0 * Radian,
                               Instant(),
                               AngularVelocity<Frame>(
                                   axis.coordinates() * Radian / Second));
  Parameters parameters(Order2ZonalCoefficient::ReadFromMessage(message.j2()));

  return std::make_unique<OblateBody<Frame>>(massive_body_parameters,
                                             rotating_body_parameters,
                                             parameters);
}

}  // namespace physics
}  // namespace principia
