#pragma once

namespace Principia {
namespace Geometry {

Sign::Sign(const bool positive) : negative_(!positive) {}

Sign::~Sign() {}

bool Sign::Negative() const { 
  return negative_;
}

bool Sign::Positive() const { 
  return !negative_;
}

Sign operator*(const Sign& left, const Sign& right) {
  return Sign(left.negative_ == right.negative_);
}

}  // namespace Geometry
}  // namespace Principia