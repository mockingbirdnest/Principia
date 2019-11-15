
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
inline Sign::Sign(double const x) : negative_(std::signbit(x)) {}

template<typename Dimensions>
Sign::Sign(Quantity<Dimensions> const& x)
    : Sign(x / SIUnit<Quantity<Dimensions>>()) {}

template<typename T, typename>
constexpr Sign Sign::OfNonZero(T x) {
  CONSTEXPR_CHECK(x != 0) << x;
  return Sign(/*negative=*/x < 0);
}

constexpr Sign Sign::Positive() {
  return Sign(/*negative=*/false);
}

constexpr Sign Sign::Negative() {
  return Sign(/*negative=*/true);
}

constexpr bool Sign::is_positive() const {
  return !negative_;
}

constexpr bool Sign::is_negative() const {
  return negative_;
}

inline constexpr Sign Sign::operator+() const {
  return *this;
}

inline constexpr Sign Sign::operator-() const {
  return Sign(!negative_);
}

inline constexpr Sign::operator int() const {
  return *this * 1;
}

constexpr bool Sign::operator==(Sign const other) const {
  return negative_ == other.negative_;
}

constexpr bool Sign::operator!=(Sign const other) const {
  return negative_ != other.negative_;
}

inline void Sign::WriteToMessage(
    not_null<serialization::Sign*> const message) const {
  message->set_negative(negative_);
}

inline Sign Sign::ReadFromMessage(serialization::Sign const message) {
  return Sign(/*negative=*/message.negative());
}

constexpr Sign::Sign(bool const negative) : negative_(negative) {}

constexpr Sign operator*(Sign const left, Sign const right) {
  return Sign(/*negative=*/left.negative_ != right.negative_);
}

template<typename T>
constexpr T operator*(Sign const left, T const& right) {
  return left.negative_ ? -right : right;
}

inline std::string DebugString(Sign const sign) {
  return sign.is_negative() ? "-" : "+";
}

inline std::ostream& operator<<(std::ostream& out, Sign const sign) {
  out << DebugString(sign);
  return out;
}

}  // namespace internal_sign
}  // namespace geometry
}  // namespace principia
