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
inline Quantity<D> Quantity<D>::SIUnit() {
  return Quantity<D>(1);
}

template<typename D>
template<int Exponent>
Exponentiation<Quantity<D>, Exponent> Quantity<D>::Pow() const {
  return Exponentiation<Quantity<D>,
                        Exponent>(pow(magnitude_, Exponent));
}

template<typename D>
inline Quantity<D>::Quantity(double const magnitude)
    : magnitude_(magnitude) {}

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
#pragma region Assigment operators
template<typename D>
inline void operator+=(Quantity<D>& left,  // NOLINT(runtime/references)
                       Quantity<D> const& right) {
  left.magnitude_ += right.magnitude_;
}

template<typename D>
inline void operator-=(Quantity<D>& left,  // NOLINT(runtime/references)
                       Quantity<D> const& right) {
  left.magnitude_ -= right.magnitude_;
}

template<typename D>
inline void operator*=(Quantity<D>& left,  // NOLINT(runtime/references)
                       double const right) {
  left.magnitude_ *= right;
}

template<typename D>
inline void operator/=(Quantity<D>& left,  // NOLINT(runtime/references)
                       double const right) {
  left.magnitude_ /= right;
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

template<typename D>
inline Quantity<D> Abs(Quantity<D> const& quantity) {
  return Quantity<D>(std::abs(quantity.magnitude_));
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

template<typename D>
inline std::string ToString(Quantity<D> const& quantity,
                            unsigned char const precision) {
  //FIX ME
  /*
  return ToString(quantity.magnitude_, precision) +
      FormatUnit("m", D::Length) + FormatUnit("kg", D::Mass) +
      FormatUnit("s", D::Time) + FormatUnit("A", D::Current) +
      FormatUnit("K", D::Temperature) + FormatUnit("mol", D::Amount) +
      FormatUnit("cd", D::LuminousIntensity) +
      FormatUnit("cycle", D::Winding) + FormatUnit("rad", D::Angle) +
      FormatUnit("sr", D::SolidAngle);*/
  return "";
}

template<typename D>
std::ostream& operator<<(std::ostream& out, Quantity<D> const& quantity) {
  return out << ToString(quantity);
}

}  // namespace quantities
}  // namespace principia
