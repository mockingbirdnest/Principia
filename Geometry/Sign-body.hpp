#pragma once

namespace principia {
namespace geometry {

template<typename Scalar>
Sign::Sign(Scalar const& scalar) : negative_(scalar < 0) {}

bool Sign::Negative() const { 
  return negative_;
}

bool Sign::Positive() const { 
  return !negative_;
}

Sign operator*(Sign const& left, Sign const& right) {
  return Sign(left.negative_ == right.negative_ ? 1 : -1);
}

template<typename Scalar>
Scalar operator*(Sign const& left, Scalar const& right) {
  return left.negative_ ? -right : right;
}

}  // namespace geometry
}  // namespace principia
