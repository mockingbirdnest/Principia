
#pragma once

#include <cmath>
#include <cstdio>
#include <string>

#include "base/macros.hpp"

namespace principia {
namespace quantities {
namespace internal_quantities {

template<int64_t LengthExponent, int64_t MassExponent, int64_t TimeExponent,
         int64_t CurrentExponent, int64_t TemperatureExponent,
         int64_t AmountExponent, int64_t LuminousIntensityExponent,
         int64_t AngleExponent>
struct Dimensions {
  enum {
    Length            = LengthExponent,
    Mass              = MassExponent,
    Time              = TimeExponent,
    Current           = CurrentExponent,
    Temperature       = TemperatureExponent,
    Amount            = AmountExponent,
    LuminousIntensity = LuminousIntensityExponent,
    Angle             = AngleExponent,
  };

  static int constexpr min_exponent = -16;
  static int constexpr max_exponent = 15;
  static int constexpr exponent_bits = 5;
  static int constexpr exponent_mask = 0x1F;

  static_assert(LengthExponent >= min_exponent &&
                LengthExponent <= max_exponent,
                "Invalid length exponent");
  static_assert(MassExponent >= min_exponent &&
                MassExponent <= max_exponent,
                "Invalid mass exponent");
  static_assert(TimeExponent >= min_exponent &&
                TimeExponent <= max_exponent,
                "Invalid time exponent");
  static_assert(CurrentExponent >= min_exponent &&
                CurrentExponent <= max_exponent,
                "Invalid current exponent");
  static_assert(TemperatureExponent >= min_exponent &&
                TemperatureExponent <= max_exponent,
                "Invalid temperature exponent");
  static_assert(AmountExponent >= min_exponent &&
                AmountExponent <= max_exponent,
                "Invalid amount exponent");
  static_assert(LuminousIntensityExponent >= min_exponent &&
                LuminousIntensityExponent <= max_exponent,
                "Invalid luminous intensity exponent");
  static_assert(AngleExponent >= min_exponent &&
                AngleExponent <= max_exponent,
                "Invalid angle exponent");

  // The NOLINT are because glint is confused by the binary and.  I kid you not.
  static int64_t constexpr representation =
      (LengthExponent & exponent_mask)                                 |  // NOLINT
      (MassExponent & exponent_mask)              << 1 * exponent_bits |  // NOLINT
      (TimeExponent & exponent_mask)              << 2 * exponent_bits |  // NOLINT
      (CurrentExponent & exponent_mask)           << 3 * exponent_bits |  // NOLINT
      (TemperatureExponent & exponent_mask)       << 4 * exponent_bits |  // NOLINT
      (AmountExponent & exponent_mask)            << 5 * exponent_bits |  // NOLINT
      (LuminousIntensityExponent & exponent_mask) << 6 * exponent_bits |  // NOLINT
      (AngleExponent & exponent_mask)             << 7 * exponent_bits;   // NOLINT
};

template<typename Q>
struct Collapse { using Type = Q; };
template<>
struct Collapse<Quantity<NoDimensions>> { using Type = double; };
template<typename Left, typename Right>
struct ProductGenerator {
  enum {
    Length            = Left::Dimensions::Length + Right::Dimensions::Length,
    Mass              = Left::Dimensions::Mass + Right::Dimensions::Mass,
    Time              = Left::Dimensions::Time + Right::Dimensions::Time,
    Current           = Left::Dimensions::Current + Right::Dimensions::Current,
    Temperature       = Left::Dimensions::Temperature +
                        Right::Dimensions::Temperature,
    Amount            = Left::Dimensions::Amount + Right::Dimensions::Amount,
    LuminousIntensity = Left::Dimensions::LuminousIntensity +
                        Right:: Dimensions::LuminousIntensity,
    Angle             = Left::Dimensions::Angle + Right::Dimensions::Angle,
  };
  using Type = typename Collapse<
      Quantity<Dimensions<Length, Mass, Time, Current, Temperature, Amount,
                          LuminousIntensity, Angle>>>::Type;
};
template<typename Left>
struct ProductGenerator<Left, double> { using Type = Left; };
template<typename Right>
struct ProductGenerator<double, Right> { using Type = Right; };
template<>
struct ProductGenerator<double, double> {
  using Type = double;
};
template<typename Left, typename Right>
struct QuotientGenerator {
  enum {
    Length            = Left::Dimensions::Length - Right::Dimensions::Length,
    Mass              = Left::Dimensions::Mass - Right::Dimensions::Mass,
    Time              = Left::Dimensions::Time - Right::Dimensions::Time,
    Current           = Left::Dimensions::Current - Right::Dimensions::Current,
    Temperature       = Left::Dimensions::Temperature -
                        Right::Dimensions::Temperature,
    Amount            = Left::Dimensions::Amount - Right::Dimensions::Amount,
    LuminousIntensity = Left::Dimensions::LuminousIntensity -
                        Right:: Dimensions::LuminousIntensity,
    Angle             = Left::Dimensions::Angle - Right::Dimensions::Angle,
  };
  using Type = typename Collapse<
      Quantity<Dimensions<Length, Mass, Time, Current, Temperature, Amount,
                          LuminousIntensity, Angle>>>::Type;
};
template<typename Left>
struct QuotientGenerator<Left, double> { using Type = Left; };
template<>
struct QuotientGenerator<double, double> {
  using Type = double;
};
template<typename Right>
struct QuotientGenerator<double, Right> {
  enum {
    Length            = -Right::Dimensions::Length,
    Mass              = -Right::Dimensions::Mass,
    Time              = -Right::Dimensions::Time,
    Current           = -Right::Dimensions::Current,
    Temperature       = -Right::Dimensions::Temperature,
    Amount            = -Right::Dimensions::Amount,
    LuminousIntensity = -Right::Dimensions::LuminousIntensity,
    Angle             = -Right::Dimensions::Angle,
  };
  using Type = Quantity<
      Dimensions<Length, Mass, Time, Current, Temperature, Amount,
                 LuminousIntensity, Angle>>;
};

template<typename T, int exponent>
struct ExponentiationGenerator<T, exponent, std::enable_if_t<(exponent > 1)>> {
  using Type = Product<typename ExponentiationGenerator<T, exponent - 1>::Type,
                       T>;
};

template<typename T, int exponent>
struct ExponentiationGenerator<T, exponent, std::enable_if_t<(exponent < 1)>>{
  using Type = Quotient<typename ExponentiationGenerator<T, exponent + 1>::Type,
                        T>;
};

template<typename T, int exponent>
struct ExponentiationGenerator<T, exponent, std::enable_if_t<(exponent == 1)>>{
  using Type = T;
};

template<typename D>
constexpr Quantity<D>::Quantity() : magnitude_(0) {}

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
FORCE_INLINE constexpr Quantity<D> Quantity<D>::operator+(
    Quantity const& right) const {
  return Quantity(magnitude_ + right.magnitude_);
}

template<typename D>
FORCE_INLINE constexpr Quantity<D> Quantity<D>::operator-(
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
  message->set_dimensions(D::representation);
  message->set_magnitude(magnitude_);
}

template<typename D>
Quantity<D> Quantity<D>::ReadFromMessage(
    serialization::Quantity const& message) {
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
FORCE_INLINE constexpr Product<Quantity<LDimensions>,
                                       Quantity<RDimensions>>
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
FORCE_INLINE constexpr Quantity<RDimensions> operator*(
    double const left,
    Quantity<RDimensions> const& right) {
  return Quantity<RDimensions>(left * right.magnitude_);
}

template<typename RDimensions>
constexpr typename Quantity<RDimensions>::Inverse operator/(
    double const left,
    Quantity<RDimensions> const& right) {
  return typename Quantity<RDimensions>::Inverse(left / right.magnitude_);
}

template<int exponent>
constexpr double Pow(double x) {
  return std::pow(x, exponent);
}

// Static specializations for frequently-used exponents, so that this gets
// turned into multiplications at compile time.

template<>
inline constexpr double Pow<-3>(double x) {
  return 1 / (x * x * x);
}

template<>
inline constexpr double Pow<-2>(double x) {
  return 1 / (x * x);
}

template<>
inline constexpr double Pow<-1>(double x) {
  return 1 / x;
}

template<>
inline constexpr double Pow<0>(double x) {
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

template<int exponent, typename D>
constexpr Exponentiation<Quantity<D>, exponent> Pow(
    Quantity<D> const& x) {
  return Exponentiation<Quantity<D>, exponent>(Pow<exponent>(x.magnitude_));
}

inline double Abs(double const x) {
  return std::abs(x);
}

template<typename D>
FORCE_INLINE Quantity<D> Abs(Quantity<D> const& quantity) {
  return Quantity<D>(std::abs(quantity.magnitude_));
}

template<typename D>
bool IsFinite(Quantity<D> const& x) {
  return std::isfinite(x.magnitude_);
}

template<typename Q>
constexpr Q SIUnit() {
  return Q(1);
}

template<>
constexpr double SIUnit<double>() {
  return 1;
}

inline std::string FormatUnit(std::string const& name, int const exponent) {
  switch (exponent) {
    case 0:
      return "";
      break;
    case 1:
      return " " + name;
    default:
      return " " + name + "^" + std::to_string(exponent);
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
  return DebugString(quantity.magnitude_, precision) +
      FormatUnit("m", D::Length) + FormatUnit("kg", D::Mass) +
      FormatUnit("s", D::Time) + FormatUnit("A", D::Current) +
      FormatUnit("K", D::Temperature) + FormatUnit("mol", D::Amount) +
      FormatUnit("cd", D::LuminousIntensity) + FormatUnit("rad", D::Angle);
}

template<typename D>
std::ostream& operator<<(std::ostream& out, Quantity<D> const& quantity) {
  return out << DebugString(quantity);
}

}  // namespace internal_quantities
}  // namespace quantities
}  // namespace principia
