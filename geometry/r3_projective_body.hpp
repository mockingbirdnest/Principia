
#pragma once

#include "geometry/r3_projective.hpp"

#include <cmath>
#include <string>

#include "quantities/quantities.hpp"

namespace principia {
namespace geometry {
namespace internal_r3_projective {

using quantities::Infinity;

template<typename Scalar>
R3Projective<Scalar>::R3Projective(
    R3Element<Scalar, Scalar, double> const& coordinates)
    : coordinates_(coordinates) {
  // (0, 0, 0) does not represent any point.
  static R3Element<Scalar, Scalar, double> const zero;
  CHECK(coordinates_ != zero);

  // Normalize the sign of |z| so that we return consistently-signed infinities
  // in the functions below.
  if (std::signbit(coordinates_.z)) {
    coordinates_.x = -coordinates_.x;
    coordinates_.y = -coordinates_.y;
    coordinates_.z = -coordinates_.z;
  }
}

template<typename Scalar>
bool R3Projective<Scalar>::is_at_infinity() const {
  return coordinates_.z == 0.0;
}

template<typename Scalar>
double R3Projective<Scalar>::point_at_infinity() const {
  CHECK(is_at_infinity());
  double const slope = coordinates_.y / coordinates_.x;
  if (std::isinf(slope)) {
    // Don't return distinct infinities in the vertical case, they would not
    // compare equal.
    return std::abs(slope);
  } else {
    return slope;
  }
}

template<typename Scalar>
Scalar const R3Projective<Scalar>::x() const {
  if (coordinates_.x == Scalar() && coordinates_.z == 0.0) {
    return Infinity<Scalar>();
  } else {
    return coordinates_.x / coordinates_.z;
  }
}

template<typename Scalar>
Scalar const R3Projective<Scalar>::y() const {
  if (coordinates_.y == Scalar() && coordinates_.z == 0.0) {
    return Infinity<Scalar>();
  } else {
    return coordinates_.y / coordinates_.z;
  }
}

template<typename Scalar>
bool operator==(R3Projective<Scalar> const& left,
                R3Projective<Scalar> const& right) {
  if (left.coordinates_.z == 0.0 && right.coordinates_.z == 0.0) {
    return left.coordinates_.x * right.coordinates_.y ==
           right.coordinates_.x * left.coordinates_.y;
  } else {
    return left.coordinates_.x * right.coordinates_.z ==
               right.coordinates_.x * left.coordinates_.z &&
           left.coordinates_.y * right.coordinates_.z ==
               right.coordinates_.y * left.coordinates_.z;
  }
}

template<typename Scalar>
bool operator!=(R3Projective<Scalar> const& left,
                R3Projective<Scalar> const& right) {
  return !(left == right);
}

template<typename Scalar>
std::string DebugString(R3Projective<Scalar> const & r3_projective) {
  return DebugString(r3_projective.coordinates_);
}

template<typename Scalar>
std::ostream& operator<<(std::ostream& os,
                         R3Projective<Scalar> const& r3_projective) {
  os << DebugString(r3_projective);
  return os;
}

}  // namespace internal_r3_projective
}  // namespace geometry
}  // namespace principia
