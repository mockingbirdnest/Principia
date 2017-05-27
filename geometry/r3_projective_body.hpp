
#pragma once

#include "geometry/projective.hpp"

namespace principia {
namespace geometry {

template<typename Scalar>
R3Projective<Scalar>::R3Projective(R3Element<Scalar> const& coordinates)
    : coordinates_(coordinates) {
  // (0, 0, 0) does not represent any point.
  CHECK(coordinates_ != R3Element<Scalar>());
}

template<typename Scalar>
bool R3Projective<Scalar>::at_infinity() const {
  return coordinates_.z == Scalar{};
}

template<typename Scalar>
inline Scalar R3Projective<Scalar>::slope() const {
  return coordinates_.y / coordinates_.x;
}

template<typename Scalar>
Scalar const R3Projective<Scalar>::x() const {
  return coordinates_.x / coordinates_.z;
}

template<typename Scalar>
Scalar const R3Projective<Scalar>::y() const {
  return coordinates_.x / coordinates_.z;
}

template<typename Scalar>
bool operator==(R3Projective<Scalar> const& left,
                R3Projective<Scalar> const& right) {
  return left.x * right.z == right.x * left.z &&
         left.y * right.z == right.y * left.z;
}

template<typename Scalar>
bool operator!=(R3Projective<Scalar> const& left,
                R3Projective<Scalar> const& right) {
  return !(left == right);
}

}  // namespace geometry
}  // namespace principia
