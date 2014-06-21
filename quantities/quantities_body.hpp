#pragma once

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
struct Collapse { typedef Q ResultType; };
template<>
struct Collapse<Quantity<NoDimensions>> { typedef double ResultType; };
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
  typedef typename Collapse<
      Quantity<Dimensions<Length, Mass, Time, Current, Temperature, Amount,
                          LuminousIntensity, Winding, Angle,
                          SolidAngle>>>::ResultType ResultType;
};
template<typename Left>
struct ProductGenerator<Left, double> { typedef Left ResultType; };
template<typename Right>
struct ProductGenerator<double, Right> { typedef Right ResultType; };
template<>
struct ProductGenerator<double, double> {
  typedef double ResultType;
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
  typedef typename Collapse<
      Quantity<Dimensions<Length, Mass, Time, Current, Temperature, Amount,
                          LuminousIntensity, Winding, Angle,
                          SolidAngle>>>::ResultType ResultType;
};
template<typename Left>
struct QuotientGenerator<Left, double> { typedef Left ResultType; };
template<>
struct QuotientGenerator<double, double> {
  typedef double ResultType;
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
  typedef Quantity<
      Dimensions<Length, Mass, Time, Current, Temperature, Amount,
                 LuminousIntensity, Winding, Angle, SolidAngle>> ResultType;
};
template<typename Q, int Exponent, typename>
struct PowerGenerator {};
template<typename Q, int Exponent>
struct PowerGenerator<Q, Exponent, Range<(Exponent > 1)>> {
  typedef Product<
      typename PowerGenerator<Q, Exponent - 1>::ResultType, Q> ResultType;
};
template<typename Q, int Exponent>
struct PowerGenerator<Q, Exponent, Range<(Exponent < 1)>>{
  typedef Quotient<
      typename PowerGenerator<Q, Exponent + 1>::ResultType, Q> ResultType;
};
template<typename Q, int Exponent>
struct PowerGenerator<Q, Exponent, Range<(Exponent == 1)>>{
  typedef Q ResultType;
};
}  // namespace type_generators

template<typename D>
inline Quantity<D>::Quantity() : magnitude_(0) {}

template<typename D>
inline Quantity<D>::Quantity(double const magnitude)
    : magnitude_(magnitude) {}

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

#pragma region Additive group

template<typename D>
inline Quantity<D> operator+(Quantity<D> const& right) {
  return Quantity<D>(+right.magnitude_);
}

template<typename D>
inline Quantity<D> operator-(Quantity<D> const& right) {
  return Quantity<D>(-right.magnitude_);
}

template<typename D>
inline Quantity<D> operator+(Quantity<D> const& left,
                             Quantity<D> const& right) {
  return Quantity<D>(left.magnitude_ + right.magnitude_);
}

template<typename D>
inline Quantity<D> operator-(Quantity<D> const& left,
                             Quantity<D> const& right) {
  return Quantity<D>(left.magnitude_ - right.magnitude_);
}

#pragma endregion

#pragma region Multiplicative group

template<typename DLeft, typename DRight>
inline Product <typename Quantity<DLeft>, typename Quantity <DRight>>
operator*(Quantity<DLeft> const& left,
          Quantity<DRight> const& right) {
  return Product<typename Quantity<DLeft>,
                 typename Quantity<DRight>>(left.magnitude_ *
                                            right.magnitude_);
}

template<typename DLeft, typename DRight>
inline Quotient<typename Quantity<DLeft>, typename Quantity <DRight>>
operator/(Quantity<DLeft> const& left,
          Quantity<DRight> const& right) {
  return Quotient<typename Quantity<DLeft>,
                  typename Quantity<DRight>>(left.magnitude_ /
                                             right.magnitude_);
}

template<typename D>
inline Quantity<D> operator*(Quantity<D> const& left,
                             double const right) {
  return Quantity<D>(left.magnitude_ * right);
}

template<typename D>
inline Quantity<D> operator*(double const left,
                             Quantity<D> const& right) {
  return Quantity<D>(left * right.magnitude_);
}

template<typename D>
inline Quantity<D> operator/(Quantity<D> const& left,
                             double const right) {
  return Quantity<D>(left.magnitude_ / right);
}

template<typename D>
inline typename Quantity<D>::Inverse operator/(double const left,
                                               Quantity<D> const& right) {
  return typename Quantity<D>::Inverse(left / right.magnitude_);
}

#pragma endregion
#pragma region Comparison operators

template<typename D>
inline bool operator>(Quantity<D> const& left, Quantity<D> const& right) {
  return left.magnitude_ > right.magnitude_;
}

template<typename D>
inline bool operator<(Quantity<D> const& left, Quantity<D> const& right) {
  return left.magnitude_ < right.magnitude_;
}

template<typename D>
inline bool operator>=(Quantity<D> const& left, Quantity<D> const& right) {
  return left.magnitude_ >= right.magnitude_;
}

template<typename D>
inline bool operator<=(Quantity<D> const& left, Quantity<D> const& right) {
  return left.magnitude_ <= right.magnitude_;
}

template<typename D>
inline bool operator==(Quantity<D> const& left, Quantity<D> const& right) {
  return left.magnitude_ == right.magnitude_;
}

template<typename D>
inline bool operator!=(Quantity<D> const& left, Quantity<D> const& right) {
  return left.magnitude_ != right.magnitude_;
}

#pragma endregion

template<int exponent>
inline double Pow(double x) {
  return std::pow(x, exponent);
}

// Static specialisations for frequently-used exponents, so that this gets
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
Exponentiation<Quantity<D>, exponent> Pow(Quantity<D> const& x) {
  return Exponentiation<Quantity<D>, exponent>(Pow<exponent>(x.magnitude_));
}

inline double Abs(double const x) {
  return std::abs(x);
}

template<typename D>
inline Quantity<D> Abs(Quantity<D> const& quantity) {
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

inline std::string DebugString(double const number,
                               unsigned char const precision) {
  char result[50];
  sprintf_s(result, ("%."+ std::to_string(precision) + "e").c_str(), number);
  return result;
}

template<typename D>
inline std::string DebugString(Quantity<D> const& quantity,
                               unsigned char const precision) {
  return DebugString(quantity.magnitude_, precision) +
      FormatUnit("m", D::Length) + FormatUnit("kg", D::Mass) +
      FormatUnit("s", D::Time) + FormatUnit("A", D::Current) +
      FormatUnit("K", D::Temperature) + FormatUnit("mol", D::Amount) +
      FormatUnit("cd", D::LuminousIntensity) +
      FormatUnit("cycle", D::Winding) + FormatUnit("rad", D::Angle) +
      FormatUnit("sr", D::SolidAngle);
}

template<typename D>
std::ostream& operator<<(std::ostream& out, Quantity<D> const& quantity) {
  return out << DebugString(quantity);
}

}  // namespace quantities
}  // namespace principia
