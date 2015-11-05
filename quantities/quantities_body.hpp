#pragma once

#include <cmath>
#include <cstdio>
#include <string>

#include "base/macros.hpp"

namespace principia {
namespace quantities {

template<int64_t LengthExponent, int64_t MassExponent, int64_t TimeExponent,
         int64_t CurrentExponent, int64_t TemperatureExponent,
         int64_t AmountExponent, int64_t LuminousIntensityExponent,
         int64_t WindingExponent, int64_t AngleExponent,
         int64_t SolidAngleExponent>
struct Dimensions {
  enum {
    Length            = LengthExponent,
    Mass              = MassExponent,
    Time              = TimeExponent,
    Current           = CurrentExponent,
    Temperature       = TemperatureExponent,
    Amount            = AmountExponent,
    LuminousIntensity = LuminousIntensityExponent,
    Winding           = WindingExponent,
    Angle             = AngleExponent,
    SolidAngle        = SolidAngleExponent
  };

  static int constexpr kMinExponent = -16;
  static int constexpr kMaxExponent = 15;
  static int constexpr kExponentBits = 5;
  static int constexpr kExponentMask = 0x1F;

  static_assert(LengthExponent >= kMinExponent &&
                LengthExponent <= kMaxExponent,
                "Invalid length exponent");
  static_assert(MassExponent >= kMinExponent &&
                MassExponent <= kMaxExponent,
                "Invalid mass exponent");
  static_assert(TimeExponent >= kMinExponent &&
                TimeExponent <= kMaxExponent,
                "Invalid time exponent");
  static_assert(CurrentExponent >= kMinExponent &&
                CurrentExponent <= kMaxExponent,
                "Invalid current exponent");
  static_assert(TemperatureExponent >= kMinExponent &&
                TemperatureExponent <= kMaxExponent,
                "Invalid temperature exponent");
  static_assert(AmountExponent >= kMinExponent &&
                AmountExponent <= kMaxExponent,
                "Invalid amount exponent");
  static_assert(LuminousIntensityExponent >= kMinExponent &&
                LuminousIntensityExponent <= kMaxExponent,
                "Invalid luminous intensity exponent");
  static_assert(AngleExponent >= kMinExponent &&
                AngleExponent <= kMaxExponent,
                "Invalid angle exponent");
  static_assert(SolidAngleExponent >= kMinExponent &&
                SolidAngleExponent <= kMaxExponent,
                "Invalid solid angle exponent");
  static_assert(WindingExponent >= kMinExponent &&
                WindingExponent <= kMaxExponent,
                "Invalid winding exponent");

  // The NOLINT are because glint is confused by the binary and.  I kid you not.
  static int64_t constexpr representation =
      (LengthExponent & kExponentMask)                                 |  // NOLINT
      (MassExponent & kExponentMask)              << 1 * kExponentBits |  // NOLINT
      (TimeExponent & kExponentMask)              << 2 * kExponentBits |  // NOLINT
      (CurrentExponent & kExponentMask)           << 3 * kExponentBits |  // NOLINT
      (TemperatureExponent & kExponentMask)       << 4 * kExponentBits |  // NOLINT
      (AmountExponent & kExponentMask)            << 5 * kExponentBits |  // NOLINT
      (LuminousIntensityExponent & kExponentMask) << 6 * kExponentBits |  // NOLINT
      (AngleExponent & kExponentMask)             << 7 * kExponentBits |  // NOLINT
      (SolidAngleExponent & kExponentMask)        << 8 * kExponentBits |  // NOLINT
      (WindingExponent & kExponentMask)           << 9 * kExponentBits;   // NOLINT
};

namespace internal {
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
    Winding           = Left::Dimensions::Winding + Right::Dimensions::Winding,
    Angle             = Left::Dimensions::Angle + Right::Dimensions::Angle,
    SolidAngle        = Left::Dimensions::SolidAngle +
                        Right::Dimensions::SolidAngle
  };
  using Type = typename Collapse<
      Quantity<Dimensions<Length, Mass, Time, Current, Temperature, Amount,
                          LuminousIntensity, Winding, Angle,
                          SolidAngle>>>::Type;
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
    Winding           = Left::Dimensions::Winding - Right::Dimensions::Winding,
    Angle             = Left::Dimensions::Angle - Right::Dimensions::Angle,
    SolidAngle        = Left::Dimensions::SolidAngle -
                        Right::Dimensions::SolidAngle
  };
  using Type = typename Collapse<
      Quantity<Dimensions<Length, Mass, Time, Current, Temperature, Amount,
                          LuminousIntensity, Winding, Angle,
                          SolidAngle>>>::Type;
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
    Winding           = -Right::Dimensions::Winding,
    Angle             = -Right::Dimensions::Angle,
    SolidAngle        = -Right::Dimensions::SolidAngle
  };
  using Type = Quantity<
      Dimensions<Length, Mass, Time, Current, Temperature, Amount,
                 LuminousIntensity, Winding, Angle, SolidAngle>>;
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

}  // namespace internal

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
constexpr internal::Product<Quantity<LDimensions>, Quantity<RDimensions>>
operator*(Quantity<LDimensions> const& left,
          Quantity<RDimensions> const& right) {
  return Product<Quantity<LDimensions>,
                 Quantity<RDimensions>>(left.magnitude_ * right.magnitude_);
}

template<typename LDimensions, typename RDimensions>
constexpr internal::Quotient<Quantity<LDimensions>, Quantity<RDimensions>>
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
Quantity<D> Abs(Quantity<D> const& quantity) {
  return Quantity<D>(std::abs(quantity.magnitude_));
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
      FormatUnit("cd", D::LuminousIntensity) + FormatUnit("cycle", D::Winding) +
      FormatUnit("rad", D::Angle) + FormatUnit("sr", D::SolidAngle);
}

template<typename D>
std::ostream& operator<<(std::ostream& out, Quantity<D> const& quantity) {
  return out << DebugString(quantity);
}

}  // namespace quantities
}  // namespace principia
