#include "physics/oblate_body.hpp"

#include <algorithm>
#include <vector>

namespace principia {
namespace physics {

namespace {
double const kNormLow = 0.999;
double const kNormHigh = 1.001;
}  // namespace

template<typename Frame>
OblateBody<Frame>::OblateBody(
    GravitationalParameter const& gravitational_parameter,
    double const j2,
    Length const& radius,
    Vector<double, Frame> const& axis)
    : OblateBody(gravitational_parameter,
                 -j2 * gravitational_parameter * radius * radius,
                 axis) {}

template<typename Frame>
OblateBody<Frame>::OblateBody(
    Mass const& mass,
    double const j2,
    Length const& radius,
    Vector<double, Frame> const& axis)
    : OblateBody(mass,
                 -j2 * mass * GravitationalConstant * radius * radius,
                 axis) {}

template<typename Frame>
OblateBody<Frame>::OblateBody(
    GravitationalParameter const& gravitational_parameter,
    Order2ZonalCoefficient const& j2,
    Vector<double, Frame> const& axis)
    : MassiveBody(gravitational_parameter),
      j2_(j2),
      axis_(axis) {
  CHECK_NE(j2, Order2ZonalCoefficient()) << "Oblate cannot have zero j2";
  CHECK_GT(axis.Norm(), kNormLow) << "Axis must have norm one";
  CHECK_LT(axis.Norm(), kNormHigh) << "Axis must have norm one";
}

template<typename Frame>
OblateBody<Frame>::OblateBody(
    Mass const& mass,
    Order2ZonalCoefficient const& j2,
    Vector<double, Frame> const& axis)
    : MassiveBody(mass),
      j2_(j2),
      axis_(axis) {
  CHECK_NE(j2, Order2ZonalCoefficient()) << "Oblate cannot have zero j2";
  CHECK_GT(axis.Norm(), kNormLow) << "Axis must have norm one";
  CHECK_LT(axis.Norm(), kNormHigh) << "Axis must have norm one";
}

template<typename Frame>
Order2ZonalCoefficient const& OblateBody<Frame>::j2() const {
  return j2_;
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
inline void OblateBody<Frame>::WriteToMessage(
    not_null<serialization::Body*> const message) const {
  WriteToMessage(message->mutable_massive_body());
}

template<typename Frame>
inline void OblateBody<Frame>::WriteToMessage(
    not_null<serialization::MassiveBody*> const message) const {
  MassiveBody::WriteToMessage(message);
  not_null<serialization::OblateBody*> const oblate_body =
      message->MutableExtension(serialization::OblateBody::oblate_body);
  Frame::WriteToMessage(oblate_body->mutable_frame());
  j2_.WriteToMessage(oblate_body->mutable_j2());
  axis_.WriteToMessage(oblate_body->mutable_axis());
}


template<typename Frame>
not_null<std::unique_ptr<OblateBody<Frame>>> OblateBody<Frame>::ReadFromMessage(
    serialization::Body const& message) {
  CHECK(message.has_massive_body());
  return ReadFromMessage(message.massive_body());
}

template<typename Frame>
not_null<std::unique_ptr<OblateBody<Frame>>> OblateBody<Frame>::ReadFromMessage(
    serialization::MassiveBody const& message) {
  CHECK(message.HasExtension(serialization::OblateBody::oblate_body));
  serialization::OblateBody const& oblateness_information =
      message.GetExtension(serialization::OblateBody::oblate_body);
  return std::make_unique<OblateBody<Frame>>(
      GravitationalParameter::ReadFromMessage(
          message.gravitational_parameter()),
      Order2ZonalCoefficient::ReadFromMessage(
          oblateness_information.j2()),
      Vector<double, Frame>::ReadFromMessage(
          oblateness_information.axis()));
}

}  // namespace physics
}  // namespace principia
