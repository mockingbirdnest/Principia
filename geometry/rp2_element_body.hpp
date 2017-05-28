
#pragma once

#include "geometry/rp2_element.hpp"

#include <cmath>
#include <string>

#include "numerics/double_precision.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace geometry {
namespace internal_rp2_element {

using numerics::TwoProduct;
using quantities::DebugString;
using quantities::Infinity;
using quantities::Square;

template<typename Scalar>
RP2Element<Scalar>::RP2Element(
    Scalar const& x, Scalar const& y, double const z)
    : x_(x), y_(y), z_(z) {
  // [0:0:0] does not represent any point but we cannot reject it as it may
  // result from an underflow or a cancellation.

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
    // Returns an infinity of the right sign.
    return Infinity<Square<Scalar>>() / x_;
  } else {
    return x_ / z_;
  }
}

template<typename Scalar>
Scalar const RP2Element<Scalar>::y() const {
  if (y_ == Scalar() && z_ == 0.0) {
    // Returns an infinity of the right sign.
    return Infinity<Square<Scalar>>() / y_;
  } else {
    return y_ / z_;
  }
}

template<typename Scalar>
bool operator==(RP2Element<Scalar> const& left,
                RP2Element<Scalar> const& right) {
  bool const left_is_singular =
      left.x_ == Scalar() && left.y_ == Scalar() && left.z_ == 0.0;
  bool const right_is_singular =
      right.x_ == Scalar() && right.y_ == Scalar() && right.z_ == 0.0;
  if (left_is_singular || right_is_singular) {
    return left_is_singular && right_is_singular;
  } else {
    if (left.z_ == 0.0 && right.z_ == 0.0) {
      return TwoProduct(left.x_, right.y_) == TwoProduct(right.x_, left.y_);
    } else {
      return TwoProduct(left.x_, right.z_) == TwoProduct(right.x_, left.z_) &&
             TwoProduct(left.y_, right.z_) == TwoProduct(right.y_, left.z_);
    }
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
