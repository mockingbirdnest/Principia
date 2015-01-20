#pragma once

#include <cmath>
#include <string>

namespace principia {
namespace quantities {

template<int LengthExponent, int MassExponent, int TimeExponent,
         int CurrentExponent, int TemperatureExponent, int AmountExponent,
         int LuminousIntensityExponent, int WindingExponent,
         int AngleExponent, int SolidAngleExponent>
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
};

namespace type_generators {
template<typename Q>
struct Collapse { using ResultType = Q; };
template<>
struct Collapse<Quantity<NoDimensions>> { using ResultType = double; };
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
  using ResultType = typename Collapse<
      Quantity<Dimensions<Length, Mass, Time, Current, Temperature, Amount,
                          LuminousIntensity, Winding, Angle,
                          SolidAngle>>>::ResultType;
};
template<typename Left>
struct ProductGenerator<Left, double> { using ResultType = Left; };
template<typename Right>
struct ProductGenerator<double, Right> { using ResultType = Right; };
template<>
struct ProductGenerator<double, double> {
  using ResultType = double;
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
  using ResultType = typename Collapse<
      Quantity<Dimensions<Length, Mass, Time, Current, Temperature, Amount,
                          LuminousIntensity, Winding, Angle,
                          SolidAngle>>>::ResultType;
};
template<typename Left>
struct QuotientGenerator<Left, double> { using ResultType = Left; };
template<>
struct QuotientGenerator<double, double> {
  using ResultType = double;
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
  using ResultType = Quantity<
      Dimensions<Length, Mass, Time, Current, Temperature, Amount,
                 LuminousIntensity, Winding, Angle, SolidAngle>>;
};
template<typename Q, int Exponent, typename>
struct PowerGenerator {};
template<typename Q, int Exponent>
struct PowerGenerator<Q, Exponent, Range<(Exponent > 1)>> {
  using ResultType =
      Product<typename PowerGenerator<Q, Exponent - 1>::ResultType, Q>;
};
template<typename Q, int Exponent>
struct PowerGenerator<Q, Exponent, Range<(Exponent < 1)>>{
  using ResultType =
      Quotient<typename PowerGenerator<Q, Exponent + 1>::ResultType, Q>;
};
template<typename Q, int Exponent>
struct PowerGenerator<Q, Exponent, Range<(Exponent == 1)>>{
  using ResultType = Q;
};
}  // namespace type_generators

template<typename D>
inline Quantity<D>::Quantity() : magnitude_(0) {}

template<typename D>
inline Quantity<D>::Quantity(double const magnitude) : magnitude_(magnitude) {}

template<typename D>
inline Quantity<D>& Quantity<D>::operator+=(Quantity const& right) {
  magnitude_ += right.magnitude_;
  return *this;
}

template<typename D>
inline Quantity<D>& Quantity<D>::operator-=(Quantity const& right) {
  magnitude_ -= right.magnitude_;
  return *this;
}

template<typename D>
inline Quantity<D>& Quantity<D>::operator*=(double const right) {
  magnitude_ *= right;
  return *this;
}

template<typename D>
inline Quantity<D>& Quantity<D>::operator/=(double const right) {
  magnitude_ /= right;
  return *this;
}

// Additive group

template<typename D>
inline Quantity<D> Quantity<D>::operator+() const {
  return *this;
}

template<typename D>
inline Quantity<D> Quantity<D>::operator-() const {
  return Quantity(-magnitude_);
}

template<typename D>
inline Quantity<D> Quantity<D>::operator+(Quantity const& right) const {
  return Quantity(magnitude_ + right.magnitude_);
}

template<typename D>
inline Quantity<D> Quantity<D>::operator-(Quantity const& right) const {
  return Quantity(magnitude_ - right.magnitude_);
}

// Comparison operators

template<typename D>
inline bool Quantity<D>::operator>(Quantity const& right) const {
  return magnitude_ > right.magnitude_;
}

template<typename D>
inline bool Quantity<D>::operator<(Quantity const& right) const {
  return magnitude_ < right.magnitude_;
}

template<typename D>
inline bool Quantity<D>::operator>=(Quantity const& right) const {
  return magnitude_ >= right.magnitude_;
}

template<typename D>
inline bool Quantity<D>::operator<=(Quantity const& right) const {
  return magnitude_ <= right.magnitude_;
}

template<typename D>
inline bool Quantity<D>::operator==(Quantity const& right) const {
  return magnitude_ == right.magnitude_;
}

template<typename D>
inline bool Quantity<D>::operator!=(Quantity const& right) const {
  return magnitude_ != right.magnitude_;
}

template<typename D>
void Quantity<D>::SerializeTo(
    not_null<serialization::Quantity*> const quantity) const {
  quantity->set_dimensions(0);
  quantity->set_magnitude(magnitude_);
}

// Multiplicative group

template<typename D>
inline Quantity<D> Quantity<D>::operator/(double const right) const {
  return Quantity(magnitude_ / right);
}

template<typename D>
inline Quantity<D> Quantity<D>::operator*(double const right) const {
  return Quantity(magnitude_ * right);
}

template<typename LDimensions, typename RDimensions>
inline Product<Quantity<LDimensions>, Quantity<RDimensions>> operator*(
    Quantity<LDimensions> const& left,
    Quantity<RDimensions> const& right) {
  return Product<Quantity<LDimensions>,
                 Quantity<RDimensions>>(left.magnitude_ * right.magnitude_);
}

template<typename LDimensions, typename RDimensions>
inline Quotient<Quantity<LDimensions>, Quantity<RDimensions>> operator/(
    Quantity<LDimensions> const& left,
    Quantity<RDimensions> const& right) {
  return Quotient<Quantity<LDimensions>,
                  Quantity<RDimensions>>(left.magnitude_ / right.magnitude_);
}

template<typename RDimensions>
inline Quantity<RDimensions> operator*(double const left,
                                       Quantity<RDimensions> const& right) {
  return Quantity<RDimensions>(left * right.magnitude_);
}

template<typename RDimensions>
inline typename Quantity<RDimensions>::Inverse operator/(
    double const left,
    Quantity<RDimensions> const& right) {
  return typename Quantity<RDimensions>::Inverse(left / right.magnitude_);
}

template<int exponent>
inline double Pow(double x) {
  return std::pow(x, exponent);
}

// Static specializations for frequently-used exponents, so that this gets
// turned into multiplications at compile time.

template<>
inline double Pow<-3>(double x) {
  return 1 / (x * x * x);
}

template<>
inline double Pow<-2>(double x) {
  return 1 / (x * x);
}

template<>
inline double Pow<-1>(double x) {
  return 1 / x;
}

template<>
inline double Pow<0>(double x) {
  return 1;
}

template<>
inline double Pow<1>(double x) {
  return x;
}

template<>
inline double Pow<2>(double x) {
  return x * x;
}

template<>
inline double Pow<3>(double x) {
  return x * x * x;
}


template<int exponent, typename D>
Exponentiation<Quantity<D>, exponent> Pow(
    Quantity<D> const& x) {
  return Exponentiation<Quantity<D>, exponent>(
      Pow<exponent>(x.magnitude_));
}

inline double Abs(double const x) {
  return std::abs(x);
}

template<typename D>
Quantity<D> Abs(Quantity<D> const& quantity) {
  return Quantity<D>(std::abs(quantity.magnitude_));
}


template<typename Q>
inline Q SIUnit() {
  return Q(1);
}

template<>
inline double SIUnit<double>() {
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
#ifdef _MSC_VER
  unsigned int old_exponent_format = _set_output_format(_TWO_DIGIT_EXPONENT);
  sprintf_s(result, ("%+."+ std::to_string(precision) + "e").c_str(), number);
  _set_output_format(old_exponent_format);
#else
  snprintf(result, sizeof(result),
           ("%."+ std::to_string(precision) + "e").c_str(), number);
#endif
  return result;
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
