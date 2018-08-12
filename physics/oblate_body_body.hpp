
#pragma once

#include "physics/oblate_body.hpp"

#include <algorithm>
#include <vector>

#include "astronomy/epoch.hpp"
#include "quantities/constants.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace internal_oblate_body {

using astronomy::J2000;
using geometry::AngularVelocity;
using geometry::Instant;
using quantities::Angle;
using quantities::si::Radian;
using quantities::si::Second;

template<typename Frame>
OblateBody<Frame>::Parameters::Parameters(Order2ZonalCoefficient const& j2)
    : j2_(j2) {
  CHECK_LT(Order2ZonalCoefficient(), j2) << "Oblate body must have positive j2";
}

template<typename Frame>
OblateBody<Frame>::Parameters::Parameters(double const j2,
                                          Length const& reference_radius)
    : j2_over_μ_(j2 * reference_radius * reference_radius) {
  CHECK_LT(0.0, j2) << "Oblate body must have positive j2";
}

template<typename Frame>
OblateBody<Frame>::Parameters::Parameters(Order2ZonalCoefficient const& j2,
                                          Order3ZonalCoefficient const& j3)
    : j2_(j2), j3_(j3) {
  CHECK_LT(Order2ZonalCoefficient(), j2) << "Oblate body must have positive j2";
  CHECK_NE(Order3ZonalCoefficient(), j3) << "Oblate body cannot have zero j3";
}

template<typename Frame>
OblateBody<Frame>::Parameters::Parameters(double const j2,
                                          double const j3,
                                          Length const& reference_radius)
    : j2_over_μ_(j2 * reference_radius * reference_radius),
      j3_over_μ_(j3 * reference_radius * reference_radius) {
  CHECK_LT(0.0, j2) << "Oblate body must have positive j2";
  CHECK_NE(0.0, j3) << "Oblate body cannot have zero j3";
}

template<typename Frame>
OblateBody<Frame>::OblateBody(
    MassiveBody::Parameters const& massive_body_parameters,
    typename RotatingBody<Frame>::Parameters const& rotating_body_parameters,
    Parameters const& parameters)
    : RotatingBody<Frame>(massive_body_parameters, rotating_body_parameters),
      parameters_(parameters) {
  if (parameters_.j2_) {
    parameters_.j2_over_μ_ = *parameters_.j2_ / this->gravitational_parameter();
  }
  if (parameters_.j2_over_μ_) {
    parameters_.j2_ = *parameters_.j2_over_μ_ * this->gravitational_parameter();
  }
  if (parameters_.j3_) {
    parameters_.j3_over_μ_ = *parameters_.j3_ / this->gravitational_parameter();
  }
  if (parameters_.j3_over_μ_) {
    parameters_.j3_ = *parameters_.j3_over_μ_ * this->gravitational_parameter();
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
Order3ZonalCoefficient const & OblateBody<Frame>::j3() const {
  return *parameters_.j3_;
}

template<typename Frame>
Quotient<Order3ZonalCoefficient, GravitationalParameter> const&
    OblateBody<Frame>::j3_over_μ() const {
  return *parameters_.j3_over_μ;
}

template<typename Frame>
bool OblateBody<Frame>::has_j3() const {
  return parameters_.j3_.has_value();
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
      message->MutableExtension(serialization::RotatingBody::extension)->
               MutableExtension(serialization::OblateBody::extension);
  parameters_.j2_->WriteToMessage(oblate_body->mutable_j2());
  if (has_j3()) {
    parameters_.j3_->WriteToMessage(oblate_body->mutable_j3());
  }
}

template<typename Frame>
not_null<std::unique_ptr<OblateBody<Frame>>> OblateBody<Frame>::ReadFromMessage(
      serialization::OblateBody const& message,
      MassiveBody::Parameters const& massive_body_parameters,
      typename RotatingBody<Frame>::Parameters const&
          rotating_body_parameters) {
  std::unique_ptr<Parameters> parameters;
  if (message.has_j2()) {
    parameters = std::make_unique<Parameters>(
        Order2ZonalCoefficient::ReadFromMessage(message.j2()),
        Order3ZonalCoefficient::ReadFromMessage(message.j3()));
  } else {
    parameters = std::make_unique<Parameters>(
        Order2ZonalCoefficient::ReadFromMessage(message.j2()));
  }

  return std::make_unique<OblateBody<Frame>>(massive_body_parameters,
                                             rotating_body_parameters,
                                             *parameters);
}

}  // namespace internal_oblate_body
}  // namespace physics
}  // namespace principia
