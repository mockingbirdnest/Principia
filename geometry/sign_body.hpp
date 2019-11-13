
#pragma once

#include "geometry/sign.hpp"

#include <cmath>
#include <string>

#include "base/macros.hpp"

namespace principia {
namespace geometry {
namespace internal_sign {

using quantities::SIUnit;

// TODO(egg): Consider intrinsics.
inline Sign::Sign(double x) : negative_(std::signbit(x)) {}

template<typename Dimensions>
Sign::Sign(Quantity<Dimensions> const& x)
    : Sign(x / SIUnit<Quantity<Dimensions>>()) {}

template<typename T, typename>
constexpr Sign Sign::OfNonZero(T x) {
  if (x == 0) {
    LOG(FATAL) << "Sign::OfNonZero(" << x << ")";
  }
  return Sign(/*negative=*/x < 0);
}

constexpr bool Sign::Negative() const {
  return negative_;
}

constexpr bool Sign::Positive() const {
  return !negative_;
}

constexpr Sign::Sign(bool negative) : negative_(negative) {}

constexpr bool Sign::operator==(Sign const& other) const {
  return negative_ == other.negative_;
}

constexpr bool Sign::operator!=(Sign const& other) const {
  return negative_ != other.negative_;
}

inline void Sign::WriteToMessage(
    not_null<serialization::Sign*> const message) const {
  message->set_negative(negative_);
}

inline Sign Sign::ReadFromMessage(serialization::Sign const& message) {
  return Sign(/*negative=*/message.negative());
}

inline Sign operator*(Sign const& left, Sign const& right) {
  return Sign(/*negative=*/left.negative_ != right.negative_);
}

template<typename T>
FORCE_INLINE(inline) T operator*(Sign const& left, T const& right) {
  return left.negative_ ? -right : right;
}

inline std::string DebugString(Sign const& sign) {
  return sign.Negative() ? "-" : "+";
}

inline std::ostream& operator<<(std::ostream& out, Sign const& sign) {
  out << DebugString(sign);
  return out;
}

}  // namespace internal_sign
}  // namespace geometry
}  // namespace principia
