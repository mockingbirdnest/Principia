// PhysicalQuantities.ipp

#pragma once

namespace PhysicalQuantities {

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

#pragma region Type generators
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
#pragma endregion
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
  return Quantity<D>(left.magnitude_ + right.magnitude_);
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
#pragma region Base quantities
inline Length Factories::Metres(Dimensionless const& number) { return Length(number); }
inline Mass Factories::Kilograms(Dimensionless const& number) { return Mass(number); }
inline Time Factories::Seconds(Dimensionless const& number) { return Time(number); }
inline Current Factories::Amperes(Dimensionless const& number) { return Current(number); }
inline Temperature Factories::Kelvins(Dimensionless const& number) { return Temperature(number); }
inline Amount Factories::Moles(Dimensionless const& number) { return Amount(number); }
inline LuminousIntensity Factories::Candelas(Dimensionless const& number) {
 return LuminousIntensity(number);
}
inline Winding Factories::Cycles(Dimensionless const& number) { return Winding(number); }
inline Angle Factories::Radians(Dimensionless const& number) { return Angle(number); }
inline SolidAngle Factories::Steradians(Dimensionless const& number) { 
  return SolidAngle(number);
}
#pragma endregion

}