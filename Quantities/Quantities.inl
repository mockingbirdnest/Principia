#pragma once

namespace Principia {
namespace Quantities {

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

namespace TypeGenerators {
template<typename Q>
struct Collapse { typedef Q ResultType; };
template<>
struct Collapse<Quantity<NoDimensions>> { typedef Dimensionless ResultType; };
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
struct ProductGenerator<Left, Dimensionless> { typedef Left ResultType; };
template<typename Right>
struct ProductGenerator<Dimensionless, Right> { typedef Right ResultType; };
template<>
struct ProductGenerator<Dimensionless, Dimensionless> {
 typedef Dimensionless ResultType;
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
struct QuotientGenerator<Left, Dimensionless> { typedef Left ResultType; };
template<>
struct QuotientGenerator<Dimensionless, Dimensionless> {
  typedef Dimensionless ResultType;
};
template<typename Right>
struct QuotientGenerator<Dimensionless, Right> {
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
}
namespace Factories {
inline Length Metres(Dimensionless const& number) { return Length(number); }
inline Mass Kilograms(Dimensionless const& number) { return Mass(number); }
inline Time Seconds(Dimensionless const& number) { return Time(number); }
inline Current Amperes(Dimensionless const& number) { return Current(number); }
inline Temperature Kelvins(Dimensionless const& number) { 
  return Temperature(number); 
}
inline Amount Moles(Dimensionless const& number) { return Amount(number); }
inline LuminousIntensity Candelas(Dimensionless const& number) {
 return LuminousIntensity(number);
}
inline Winding Cycles(Dimensionless const& number) { return Winding(number); }
inline Angle Radians(Dimensionless const& number) { return Angle(number); }
inline SolidAngle Steradians(Dimensionless const& number) { 
  return SolidAngle(number);
}
}
template<typename D>
template<int Exponent>
Exponentiation<Quantity<D>, Exponent> Quantity<D>::Pow() const {
  return Exponentiation<Quantity<D>, 
                        Exponent>(magnitude_.Pow(Exponent));
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
                             Dimensionless const& right) {
  return Quantity<D>(left.magnitude_ * right);
}
template<typename D>
inline Quantity<D> operator*(Dimensionless const& left, 
                             Quantity<D> const& right) {
  return Quantity<D>(left * right.magnitude_);
}
template<typename D>
inline Quantity<D> operator/(Quantity<D> const& left,
                             Dimensionless const& right) {
  return Quantity<D>(left.magnitude_ / right);
}
template<typename D>
inline Inverse<Quantity<D>> operator/(Dimensionless const& left,
                                      Quantity<D> const& right) {
  return Inverse<Quantity<D>>(left / right.magnitude_);
}
#pragma endregion
#pragma region Assigment operators
template<typename D>
inline void operator+=(Quantity<D>& left, Quantity<D> const& right) {
  left = left + right;
}
template<typename D>
inline void operator-=(Quantity<D>& left, Quantity<D> const& right) {
  left = left - right;
}
template<typename D>
inline void operator*=(Quantity<D>& left, Dimensionless const& right) {
  left = left * right;
}
template<typename D>
inline void operator/=(Quantity<D>& left, Dimensionless const& right) {
  left = left / right;
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
  return Quantity<D>(Abs(quantity.magnitude_));
}

template<typename D>
inline std::wstring ToString(Quantity<D> const& quantity, 
                             unsigned char const precision) {
  return ToString(quantity.magnitude_, precision) +
    (D::Length            != 0 ? L" m^" + std::to_wstring(D::Length) : L"") +
    (D::Mass              != 0 ? L" kg^" + std::to_wstring(D::Mass) : L"") +
    (D::Time              != 0 ? L" s^" + std::to_wstring(D::Time) : L"") +
    (D::Current           != 0 ? L" A^" + std::to_wstring(D::Current) : L"") +
    (D::Temperature       != 0 ? L" K^" + std::to_wstring(D::Temperature) 
                               : L"") +
    (D::Amount            != 0 ? L" mol^" + std::to_wstring(D::Amount) : L"") +
    (D::LuminousIntensity != 0 ? L" cd^" + std::to_wstring(D::LuminousIntensity) 
                               : L"") +
    (D::Winding           != 0 ? L" cycle^" + std::to_wstring(D::Winding) 
                               : L"") +
    (D::Angle             != 0 ? L" rad^" + std::to_wstring(D::Angle) : L"") +
    (D::SolidAngle        != 0 ? L" sr^" + std::to_wstring(D::SolidAngle) 
                               : L"");
}
}
}
