
#pragma once

#include "quantities/quantities.hpp"

#include <cmath>
#include <cstdio>
#include <limits>
#include <string>

#include "base/macros.hpp"

namespace principia {
namespace quantities {
namespace internal_quantities {

using internal_dimensions::DimensionsAreSerializable;

template<typename D>
constexpr Quantity<D>::Quantity() : magnitude_(0) {}

template<typename D>
constexpr Quantity<D>::Quantity(uninitialized_t) {}

template<typename D>
constexpr Quantity<D>::Quantity(double const magnitude)
    : magnitude_(magnitude) {}

template<typename D>
Quantity<D>& Quantity<D>::operator+=(Quantity const& right) {
  magnitude_ += right.magnitude_;
  return *this;
}

template<typename D>
Quantity<D>& Quantity<D>::operator-=(Quantity const& right) {
  magnitude_ -= right.magnitude_;
  return *this;
}

template<typename D>
Quantity<D>& Quantity<D>::operator*=(double const right) {
  magnitude_ *= right;
  return *this;
}

template<typename D>
Quantity<D>& Quantity<D>::operator/=(double const right) {
  magnitude_ /= right;
  return *this;
}

// Additive group

template<typename D>
constexpr Quantity<D> Quantity<D>::operator+() const {
  return *this;
}

template<typename D>
constexpr Quantity<D> Quantity<D>::operator-() const {
  return Quantity(-magnitude_);
}

template<typename D>
FORCE_INLINE(constexpr) Quantity<D> Quantity<D>::operator+(
    Quantity const& right) const {
  return Quantity(magnitude_ + right.magnitude_);
}

template<typename D>
FORCE_INLINE(constexpr) Quantity<D> Quantity<D>::operator-(
    Quantity const& right) const {
  return Quantity(magnitude_ - right.magnitude_);
}

// Comparison operators

template<typename D>
constexpr bool Quantity<D>::operator>(Quantity const& right) const {
  return magnitude_ > right.magnitude_;
}

template<typename D>
constexpr bool Quantity<D>::operator<(Quantity const& right) const {
  return magnitude_ < right.magnitude_;
}

template<typename D>
constexpr bool Quantity<D>::operator>=(Quantity const& right) const {
  return magnitude_ >= right.magnitude_;
}

template<typename D>
constexpr bool Quantity<D>::operator<=(Quantity const& right) const {
  return magnitude_ <= right.magnitude_;
}

template<typename D>
constexpr bool Quantity<D>::operator==(Quantity const& right) const {
  return magnitude_ == right.magnitude_;
}

template<typename D>
constexpr bool Quantity<D>::operator!=(Quantity const& right) const {
  return magnitude_ != right.magnitude_;
}

template<typename D>
void Quantity<D>::WriteToMessage(
    not_null<serialization::Quantity*> const message) const {
  static_assert(internal_dimensions::DimensionsAreSerializable<D>::value,
                "Failed to check serializability");
  message->set_dimensions(D::representation);
  message->set_magnitude(magnitude_);
}

template<typename D>
Quantity<D> Quantity<D>::ReadFromMessage(
    serialization::Quantity const& message) {
  static_assert(internal_dimensions::DimensionsAreSerializable<D>::value,
                "Failed to check serializability");
  CHECK_EQ(D::representation, message.dimensions());
  return Quantity(message.magnitude());
}

// Multiplicative group

template<typename D>
constexpr Quantity<D> Quantity<D>::operator/(double const right) const {
  return Quantity(magnitude_ / right);
}

template<typename D>
constexpr Quantity<D> Quantity<D>::operator*(double const right) const {
  return Quantity(magnitude_ * right);
}

template<typename LDimensions, typename RDimensions>
FORCE_INLINE(constexpr) Product<Quantity<LDimensions>, Quantity<RDimensions>>
operator*(Quantity<LDimensions> const& left,
          Quantity<RDimensions> const& right) {
  return Product<Quantity<LDimensions>,
                 Quantity<RDimensions>>(left.magnitude_ * right.magnitude_);
}

template<typename LDimensions, typename RDimensions>
constexpr Quotient<Quantity<LDimensions>, Quantity<RDimensions>>
operator/(Quantity<LDimensions> const& left,
          Quantity<RDimensions> const& right) {
  return Quotient<Quantity<LDimensions>,
                  Quantity<RDimensions>>(left.magnitude_ / right.magnitude_);
}

template<typename RDimensions>
FORCE_INLINE(constexpr) Quantity<RDimensions> operator*(
    double const left,
    Quantity<RDimensions> const& right) {
  return Quantity<RDimensions>(left * right.magnitude_);
}

template<typename RDimensions>
constexpr Quotient<double, Quantity<RDimensions>> operator/(
    double const left,
    Quantity<RDimensions> const& right) {
  return Quotient<double, Quantity<RDimensions>>(left / right.magnitude_);
}

template<typename Q>
constexpr Q SIUnit() {
  static_assert(is_quantity<Q>::value, "Not a quantity");
  return Q(1);
}

template<>
constexpr double SIUnit<double>() {
  return 1;
}

template<typename Q>
constexpr Q Infinity() {
  return SIUnit<Q>() * std::numeric_limits<double>::infinity();
}

template<typename Q>
constexpr bool IsFinite(Q const& x) {
  return std::isfinite(x / SIUnit<Q>());
}

template<typename Q>
constexpr Q NaN() {
  return SIUnit<Q>() * std::numeric_limits<double>::quiet_NaN();
}

template<typename D>
std::string Format() {
  auto const format_unit = [](std::string const& name,
                              int const exponent) -> std::string {
    switch (exponent) {
      case 0:
        return "";
        break;
      case 1:
        return " " + name;
      default:
        return " " + name + "^" + std::to_string(exponent);
    }
  };

  // This string has a leading space if it's not empty.
  auto const format =
      format_unit("m", D::Length) + format_unit("kg", D::Mass) +
      format_unit("s", D::Time) + format_unit("A", D::Current) +
      format_unit("K", D::Temperature) + format_unit("mol", D::Amount) +
      format_unit("cd", D::LuminousIntensity) + format_unit("rad", D::Angle);

  if (format.empty()) {
    return format;
  } else {
    return format.substr(1, format.size() - 1);
  }
}

inline std::string DebugString(double const number, int const precision) {
  char result[50];
#if OS_WIN && PRINCIPIA_COMPILER_MSVC && (_MSC_VER < 1900)
  unsigned int old_exponent_format = _set_output_format(_TWO_DIGIT_EXPONENT);
  int const size = sprintf_s(result,
                             ("%+." + std::to_string(precision) + "e").c_str(),
                             number);
  _set_output_format(old_exponent_format);
#else
  int const size = snprintf(result, sizeof(result),
                            ("%+." + std::to_string(precision) + "e").c_str(),
                            number);
#endif
  CHECK_LE(0, size);
  return std::string(result, size);
}

template<typename D>
std::string DebugString(Quantity<D> const& quantity, int const precision) {
  return DebugString(quantity / SIUnit<Quantity<D>>(), precision) + " " +
         Format<D>();
}

template<typename D>
std::ostream& operator<<(std::ostream& out, Quantity<D> const& quantity) {
  return out << DebugString(quantity);
}

}  // namespace internal_quantities
}  // namespace quantities
}  // namespace principia
