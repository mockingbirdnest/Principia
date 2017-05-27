
#pragma once

#include "geometry/rp2_element.hpp"

#include <cmath>
#include <string>

#include "quantities/quantities.hpp"

namespace principia {
namespace geometry {
namespace internal_rp2_element {

using quantities::DebugString;
using quantities::Infinity;

template<typename Scalar>
RP2Element<Scalar>::RP2Element(
    Scalar const& x, Scalar const& y, double z)
    : x_(x), y_(y), z_(z) {
  // [0:0:0] does not represent any point.
  CHECK(x != Scalar() || y != Scalar() || z != 0.0);

  // Normalize the sign of |z_| so that we return consistently-signed infinities
  // in the functions below.
  if (std::signbit(z_)) {
    x_ = -x_;
    y_ = -y_;
    z_ = -z_;
  }
}

template<typename Scalar>
bool RP2Element<Scalar>::is_at_infinity() const {
  return z_ == 0.0;
}

template<typename Scalar>
Scalar const RP2Element<Scalar>::x() const {
  if (x_ == Scalar() && z_ == 0.0) {
    return Infinity<Scalar>();
  } else {
    return x_ / z_;
  }
}

template<typename Scalar>
Scalar const RP2Element<Scalar>::y() const {
  if (y_ == Scalar() && z_ == 0.0) {
    return Infinity<Scalar>();
  } else {
    return y_ / z_;
  }
}

template<typename Scalar>
bool operator==(RP2Element<Scalar> const& left,
                RP2Element<Scalar> const& right) {
  if (left.z_ == 0.0 && right.z_ == 0.0) {
    return left.x_ * right.y_ == right.x_ * left.y_;
  } else {
    return left.x_ * right.z_ == right.x_ * left.z_ &&
           left.y_ * right.z_ == right.y_ * left.z_;
  }
}

template<typename Scalar>
bool operator!=(RP2Element<Scalar> const& left,
                RP2Element<Scalar> const& right) {
  return !(left == right);
}

template<typename Scalar>
std::string DebugString(RP2Element<Scalar> const & rp2_element) {
  return "[" + DebugString(rp2_element.x_) + ":" +
               DebugString(rp2_element.y_) + ":" +
               DebugString(rp2_element.z_) + "]";
}

template<typename Scalar>
std::ostream& operator<<(std::ostream& os,
                         RP2Element<Scalar> const& rp2_element) {
  os << DebugString(rp2_element);
  return os;
}

}  // namespace internal_rp2_element
}  // namespace geometry
}  // namespace principia
