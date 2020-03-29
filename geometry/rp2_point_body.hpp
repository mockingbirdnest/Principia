
#pragma once

#include "geometry/rp2_point.hpp"

#include <cmath>
#include <string>

#include "numerics/double_precision.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace geometry {
namespace internal_rp2_point {

using numerics::TwoProduct;
using quantities::DebugString;
using quantities::Infinity;
using quantities::Square;

template<typename Scalar, typename Frame>
RP2Point<Scalar, Frame>::RP2Point(
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

template<typename Scalar, typename Frame>
bool RP2Point<Scalar, Frame>::is_at_infinity() const {
  return z_ == 0.0;
}

template<typename Scalar, typename Frame>
Scalar RP2Point<Scalar, Frame>::x() const {
  if (x_ == Scalar() && z_ == 0.0) {
    // Returns an infinity of the right sign.
    return Infinity<Square<Scalar>>() / x_;
  } else {
    return x_ / z_;
  }
}

template<typename Scalar, typename Frame>
Scalar RP2Point<Scalar, Frame>::y() const {
  if (y_ == Scalar() && z_ == 0.0) {
    // Returns an infinity of the right sign.
    return Infinity<Square<Scalar>>() / y_;
  } else {
    return y_ / z_;
  }
}

template<typename Scalar, typename Frame>
bool operator==(RP2Point<Scalar, Frame> const& left,
                RP2Point<Scalar, Frame> const& right) {
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

template<typename Scalar, typename Frame>
bool operator!=(RP2Point<Scalar, Frame> const& left,
                RP2Point<Scalar, Frame> const& right) {
  return !(left == right);
}

template<typename Scalar, typename Frame>
std::string DebugString(RP2Point<Scalar, Frame> const & rp2_point) {
  return "[" + DebugString(rp2_point.x_) + ":" +
               DebugString(rp2_point.y_) + ":" +
               DebugString(rp2_point.z_) + "]";
}

template<typename Scalar, typename Frame>
std::ostream& operator<<(std::ostream& os,
                         RP2Point<Scalar, Frame> const& rp2_point) {
  os << DebugString(rp2_point);
  return os;
}

}  // namespace internal_rp2_point
}  // namespace geometry
}  // namespace principia
