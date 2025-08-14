#pragma once

#include "quantities/quantities.hpp"

#include <cmath>
#include <cstdio>
#include <limits>
#include <string>

namespace principia {
namespace quantities {
namespace _quantities {
namespace internal {

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

template<typename D>
void Quantity<D>::WriteToMessage(
    not_null<serialization::Quantity*> const message) const {
  static_assert(DimensionsAreSerializable<D>::value,
                "Failed to check serializability");
  message->set_dimensions(D::representation);
  message->set_magnitude(magnitude_);
}

template<typename D>
Quantity<D> Quantity<D>::ReadFromMessage(
    serialization::Quantity const& message) {
  static_assert(DimensionsAreSerializable<D>::value,
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
    double left,
    Quantity<RDimensions> const& right) {
  return Quotient<double, Quantity<RDimensions>>(left / right.magnitude_);
}

// Static specializations for frequently-used exponents, so that this gets
// turned into multiplications at compile time.

template<int exponent>
constexpr double Pow(double x) {
  // Use the Russian peasant algorithm for small exponents.
  if constexpr (exponent > 0 && exponent < 32) {
    // The end of the recursion is handled by the specializations below.
    auto const y = Pow<exponent / 2>(x);
    auto const y² = y * y;
    if constexpr (exponent % 2 == 1) {
      return y² * x;
    } else {
      return y²;
    }
  } else if constexpr (exponent < 0 && exponent > -32) {
    return 1 / Pow<-exponent>(x);
  } else {
    return std::pow(x, exponent);
  }
}

template<>
inline constexpr double Pow<0>(double) {
  return 1;
}

template<>
inline constexpr double Pow<1>(double x) {
  return x;
}

template<>
inline constexpr double Pow<2>(double x) {
  return x * x;
}

template<>
inline constexpr double Pow<3>(double x) {
  return x * x * x;
}

template<int exponent, typename Q>
constexpr Exponentiation<Q, exponent> Pow(Q const& x) {
  if constexpr (boost_cpp_rational<Q>) {
    // It seems that Boost does not define `pow` for `cpp_rational`.
    return cpp_rational(pow(numerator(x), exponent),
                        pow(denominator(x), exponent));
  } else if constexpr (boost_cpp_number<Q>) {
    return pow(x, exponent);
  } else {
    return SIUnit<Exponentiation<Q, exponent>>() *
           Pow<exponent>(x / SIUnit<Q>());
  }
}

inline __m128d ToM128D(double const x) {
  return _mm_set1_pd(x);
}

template<typename Dimensions>
__m128d ToM128D(Quantity<Dimensions> const x) {
  return _mm_set1_pd(x.magnitude_);
}

template<typename Q>
constexpr bool IsFinite(Q const& x) {
  return std::isfinite(x / SIUnit<Q>());
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
  std::string result;
  result.resize(50);
  int const size = snprintf(result.data(), result.size(),
                            ("%+." + std::to_string(precision) + "e").c_str(),
                            number);
  CHECK_LE(0, size);
  result.resize(size);
  return result;
}

inline std::string DebugString(M128D const number, int const precision) {
  return DebugString(static_cast<double>(number), precision);
}

template<typename D>
std::string DebugString(Quantity<D> const& quantity, int const precision) {
  return DebugString(quantity / SIUnit<Quantity<D>>(), precision) + " " +
         Format<D>();
}

template<boost_cpp_number N>
std::string DebugString(N const& number, int const precision) {
  return number.str();
}

template<typename D>
std::ostream& operator<<(std::ostream& out, Quantity<D> const& quantity) {
  return out << DebugString(quantity);
}

}  // namespace internal
}  // namespace _quantities
}  // namespace quantities
}  // namespace principia
