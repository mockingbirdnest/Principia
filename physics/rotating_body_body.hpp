
#pragma once

#include "physics/rotating_body.hpp"

#include <algorithm>
#include <vector>

#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/rotation.hpp"
#include "physics/oblate_body.hpp"
#include "quantities/constants.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace physics {
namespace internal_rotating_body {

using geometry::Cross;
using geometry::DefinesFrame;
using geometry::EulerAngles;
using geometry::Exp;
using geometry::NormalizeOrZero;
using geometry::RadiusLatitudeLongitude;
using geometry::SphericalCoordinates;
using quantities::SIUnit;
using quantities::si::Radian;

template<typename Frame>
RotatingBody<Frame>::Parameters::Parameters(
    Length const& mean_radius,
    Angle const& reference_angle,
    Instant const& reference_instant,
    AngularFrequency const& angular_frequency,
    Angle const& right_ascension_of_pole,
    Angle const& declination_of_pole)
    : mean_radius_(mean_radius),
      reference_angle_(reference_angle),
      reference_instant_(reference_instant),
      angular_frequency_(angular_frequency),
      right_ascension_of_pole_(right_ascension_of_pole),
      declination_of_pole_(declination_of_pole) {
  CHECK_NE(angular_frequency_, 0.0 * SIUnit<AngularFrequency>())
      << "Rotating body cannot have zero angular velocity";
}

template<typename Frame>
RotatingBody<Frame>::RotatingBody(
    MassiveBody::Parameters const& massive_body_parameters,
    Parameters const& parameters)
    : MassiveBody(massive_body_parameters),
      parameters_(parameters),
      polar_axis_(RadiusLatitudeLongitude(
                      1.0,
                      parameters.declination_of_pole_,
                      parameters.right_ascension_of_pole_).ToCartesian()),
      angular_velocity_(polar_axis_.coordinates() *
                        parameters.angular_frequency_) {
  biequatorial_ = NormalizeOrZero(Vector<double, Frame>(
      {-polar_axis_.coordinates().y,
       polar_axis_.coordinates().x,
       0}));
  if (biequatorial_ == Vector<double, Frame>{}) {
    biequatorial_ = Vector<double, Frame>({1, 0, 0});
    equatorial_ = Vector<double, Frame>({0, 1, 0});
  } else {
    // TODO(phl): It is somewhat unpleasant that we have to make this a vector
    // when it would want to be a bivector.  The inner products in Geopotential
    // would then yield trivectors and we are not sure what that means.
    equatorial_ = Vector<double, Frame>(Cross(polar_axis_.coordinates(),
                                              biequatorial_.coordinates()));
  }
}

template<typename Frame>
Length RotatingBody<Frame>::mean_radius() const {
  return parameters_.mean_radius_;
}

template<typename Frame>
Vector<double, Frame> const& RotatingBody<Frame>::polar_axis() const {
  return polar_axis_;
}

template<typename Frame>
Vector<double, Frame> const& RotatingBody<Frame>::biequatorial() const {
  return biequatorial_;
}

template<typename Frame>
Vector<double, Frame> const& RotatingBody<Frame>::equatorial() const {
  return equatorial_;
}

template<typename Frame>
Angle const& RotatingBody<Frame>::right_ascension_of_pole() const {
  return parameters_.right_ascension_of_pole_;
}


template<typename Frame>
Angle const& RotatingBody<Frame>::declination_of_pole() const {
  return parameters_.declination_of_pole_;
}

template<typename Frame>
AngularFrequency const& RotatingBody<Frame>::angular_frequency() const {
  return parameters_.angular_frequency_;
}

template<typename Frame>
AngularVelocity<Frame> const& RotatingBody<Frame>::angular_velocity() const {
  return angular_velocity_;
}

template<typename Frame>
Angle RotatingBody<Frame>::AngleAt(Instant const& t) const {
  return parameters_.reference_angle_ +
         (t - parameters_.reference_instant_) * parameters_.angular_frequency_;
}

template<typename Frame>
Rotation<Frame, Frame> RotatingBody<Frame>::RotationAt(Instant const& t) const {
  return Exp((t - parameters_.reference_instant_) * angular_velocity_);
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
  parameters_.mean_radius_.WriteToMessage(
      rotating_body->mutable_mean_radius());
  parameters_.reference_angle_.WriteToMessage(
      rotating_body->mutable_reference_angle());
  parameters_.reference_instant_.WriteToMessage(
      rotating_body->mutable_reference_instant());
  parameters_.angular_frequency_.WriteToMessage(
      rotating_body->mutable_angular_frequency());
  parameters_.right_ascension_of_pole_.WriteToMessage(
      rotating_body->mutable_right_ascension_of_pole());
  parameters_.declination_of_pole_.WriteToMessage(
      rotating_body->mutable_declination_of_pole());
}

template<typename Frame>
not_null<std::unique_ptr<RotatingBody<Frame>>>
RotatingBody<Frame>::ReadFromMessage(
    serialization::RotatingBody const& message,
    MassiveBody::Parameters const& massive_body_parameters) {
  Parameters const parameters(
      Length::ReadFromMessage(message.mean_radius()),
      Angle::ReadFromMessage(message.reference_angle()),
      Instant::ReadFromMessage(message.reference_instant()),
      AngularFrequency::ReadFromMessage(message.angular_frequency()),
      Angle::ReadFromMessage(message.right_ascension_of_pole()),
      Angle::ReadFromMessage(message.declination_of_pole()));

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

}  // namespace internal_rotating_body
}  // namespace physics
}  // namespace principia
