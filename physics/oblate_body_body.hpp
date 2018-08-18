
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
OblateBody<Frame>::Parameters::Parameters(
    Degree2SphericalHarmonicCoefficient const& j2) : j2_(j2) {
  CHECK_LT(Degree2SphericalHarmonicCoefficient(), j2)
      << "Oblate body must have positive j2";
}

template<typename Frame>
OblateBody<Frame>::Parameters::Parameters(double const j2,
                                          Length const& reference_radius)
    : j2_over_μ_(j2 * reference_radius * reference_radius) {
  CHECK_LT(0.0, j2) << "Oblate body must have positive j2";
}

template<typename Frame>
OblateBody<Frame>::Parameters::Parameters(
    Degree2SphericalHarmonicCoefficient const& j2,
    Degree2SphericalHarmonicCoefficient const& c22,
    Degree2SphericalHarmonicCoefficient const& s22)
    : j2_(j2), c22_(c22), s22_(s22) {
  CHECK_LT(Degree2SphericalHarmonicCoefficient(), j2)
      << "Oblate body must have positive j2";
  CHECK_NE(Degree2SphericalHarmonicCoefficient(), c22)
      << "Oblate body cannot have zero c22";
  CHECK_NE(Degree2SphericalHarmonicCoefficient(), s22)
      << "Oblate body cannot have zero s22";
}

template<typename Frame>
OblateBody<Frame>::Parameters::Parameters(double const j2,
                                          double const c22,
                                          double const s22,
                                          Length const& reference_radius)
    : j2_over_μ_(j2 * reference_radius * reference_radius),
      c22_over_μ_(c22 * reference_radius * reference_radius),
      s22_over_μ_(s22 * reference_radius * reference_radius) {
  CHECK_LT(0.0, j2) << "Oblate body must have positive j2";
  CHECK_NE(0.0, c22) << "Oblate body cannot have zero c22";
  CHECK_NE(0.0, s22) << "Oblate body cannot have zero s22";
}

template<typename Frame>
OblateBody<Frame>::Parameters::Parameters(
    Degree2SphericalHarmonicCoefficient const& j2,
    Degree2SphericalHarmonicCoefficient const& c22,
    Degree2SphericalHarmonicCoefficient const& s22,
    Degree3SphericalHarmonicCoefficient const& j3)
    : j2_(j2), c22_(c22), s22_(s22), j3_(j3) {
  CHECK_LT(Degree2SphericalHarmonicCoefficient(), j2)
      << "Oblate body must have positive j2";
  CHECK_NE(Degree2SphericalHarmonicCoefficient(), c22)
      << "Oblate body cannot have zero c22";
  CHECK_NE(Degree2SphericalHarmonicCoefficient(), s22)
      << "Oblate body cannot have zero s22";
  CHECK_NE(Degree3SphericalHarmonicCoefficient(), j3)
      << "Oblate body cannot have zero j3";
}

template<typename Frame>
OblateBody<Frame>::Parameters::Parameters(double const j2,
                                          double const c22,
                                          double const s22,
                                          double const j3,
                                          Length const& reference_radius)
    : j2_over_μ_(j2 * reference_radius * reference_radius),
      c22_over_μ_(c22 * reference_radius * reference_radius),
      s22_over_μ_(s22 * reference_radius * reference_radius),
      j3_over_μ_(j3 * reference_radius * reference_radius * reference_radius) {
  CHECK_LT(0.0, j2) << "Oblate body must have positive j2";
  CHECK_NE(0.0, c22) << "Oblate body cannot have zero c22";
  CHECK_NE(0.0, s22) << "Oblate body cannot have zero s22";
  CHECK_NE(0.0, j3) << "Oblate body cannot have zero j3";
}

#define PRINCIPIA_FILL_OBLATE_BODY_PARAMETERS(name)                    \
  if (parameters_.name##_) {                                           \
    parameters_.name##_over_μ_ =                                       \
        *parameters_.name##_ / this->gravitational_parameter();        \
  }                                                                    \
  if (parameters_.name##_over_μ_) {                                    \
    parameters_.name##_ =                                              \
        *parameters_.name##_over_μ_ * this->gravitational_parameter(); \
  }

template<typename Frame>
OblateBody<Frame>::OblateBody(
    MassiveBody::Parameters const& massive_body_parameters,
    typename RotatingBody<Frame>::Parameters const& rotating_body_parameters,
    Parameters const& parameters)
    : RotatingBody<Frame>(massive_body_parameters, rotating_body_parameters),
      parameters_(parameters) {
  PRINCIPIA_FILL_OBLATE_BODY_PARAMETERS(j2);
  PRINCIPIA_FILL_OBLATE_BODY_PARAMETERS(c22);
  PRINCIPIA_FILL_OBLATE_BODY_PARAMETERS(s22);
  PRINCIPIA_FILL_OBLATE_BODY_PARAMETERS(j3);
}

#undef PRINCIPIA_FILL_OBLATE_BODY_PARAMETERS

template<typename Frame>
Degree2SphericalHarmonicCoefficient const& OblateBody<Frame>::j2() const {
  return *parameters_.j2_;
}

template<typename Frame>
Quotient<Degree2SphericalHarmonicCoefficient,
         GravitationalParameter> const& OblateBody<Frame>::j2_over_μ() const {
  return *parameters_.j2_over_μ_;
}

template<typename Frame>
Degree2SphericalHarmonicCoefficient const OblateBody<Frame>::c22() const {
  return parameters_.c22_.value_or(Degree2SphericalHarmonicCoefficient());
}

template<typename Frame>
Quotient<Degree2SphericalHarmonicCoefficient, GravitationalParameter> const
    OblateBody<Frame>::c22_over_μ() const {
  return parameters_.c22_over_μ_.value_or(
      Quotient<Degree2SphericalHarmonicCoefficient, GravitationalParameter>());
}

template<typename Frame>
Degree2SphericalHarmonicCoefficient const OblateBody<Frame>::s22() const {
  return parameters_.s22_.value_or(Degree2SphericalHarmonicCoefficient());
}

template<typename Frame>
Quotient<Degree2SphericalHarmonicCoefficient, GravitationalParameter> const
    OblateBody<Frame>::s22_over_μ() const {
  return parameters_.s22_over_μ_.value_or(
      Quotient<Degree2SphericalHarmonicCoefficient, GravitationalParameter>());
}

template<typename Frame>
Degree3SphericalHarmonicCoefficient const OblateBody<Frame>::j3() const {
  return parameters_.j3_.value_or(Degree3SphericalHarmonicCoefficient());
}

template<typename Frame>
Quotient<Degree3SphericalHarmonicCoefficient, GravitationalParameter> const
    OblateBody<Frame>::j3_over_μ() const {
  return parameters_.j3_over_μ_.value_or(
      Quotient<Degree3SphericalHarmonicCoefficient, GravitationalParameter>());
}

template<typename Frame>
bool OblateBody<Frame>::has_c22() const {
  return parameters_.c22_.has_value();
}

template<typename Frame>
bool OblateBody<Frame>::has_s22() const {
  return parameters_.s22_.has_value();
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
  if (has_c22()) {
    parameters_.c22_->WriteToMessage(oblate_body->mutable_c22());
  }
  if (has_s22()) {
    parameters_.s22_->WriteToMessage(oblate_body->mutable_s22());
  }
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
  if (message.has_c22() || message.has_s22()) {
    if (message.has_j3()) {
      parameters = std::make_unique<Parameters>(
          Degree2SphericalHarmonicCoefficient::ReadFromMessage(message.j2()),
          Degree2SphericalHarmonicCoefficient::ReadFromMessage(message.c22()),
          Degree2SphericalHarmonicCoefficient::ReadFromMessage(message.s22()),
          Degree3SphericalHarmonicCoefficient::ReadFromMessage(message.j3()));
    } else {
      parameters = std::make_unique<Parameters>(
          Degree2SphericalHarmonicCoefficient::ReadFromMessage(message.j2()),
          Degree2SphericalHarmonicCoefficient::ReadFromMessage(message.c22()),
          Degree2SphericalHarmonicCoefficient::ReadFromMessage(message.s22()));
    }
  } else {
    parameters = std::make_unique<Parameters>(
        Degree2SphericalHarmonicCoefficient::ReadFromMessage(message.j2()));
  }

  return std::make_unique<OblateBody<Frame>>(massive_body_parameters,
                                             rotating_body_parameters,
                                             *parameters);
}

}  // namespace internal_oblate_body
}  // namespace physics
}  // namespace principia
