#include "physics/oblate_body.hpp"

#include <algorithm>
#include <vector>

namespace principia {
namespace physics {

namespace {
double const kNormLow = kNormLow;
double const kNormHigh = kNormHigh;
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

}  // namespace physics
}  // namespace principia
