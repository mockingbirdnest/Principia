// PhysicalQuantities-body-inl.h

#pragma once

template<int LengthExponent, int MassExponent, int TimeExponent,
         int CurrentExponent, int TemperatureExponent, int AmountExponent,
         int LuminousIntensityExponent, int WindingExponent,
         int WrappingExponent>
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
    Wrapping          = WrappingExponent
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
    Wrapping          = Left::Dimensions::Wrapping + Right::Dimensions::Wrapping
  };
  typedef typename Collapse<Quantity<Dimensions<Length, Mass, Time, Current, Temperature,
                                       Amount, LuminousIntensity, Winding,
                                       Wrapping>>>::ResultType ResultType;
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
    Wrapping          = Left::Dimensions::Wrapping - Right::Dimensions::Wrapping
  };
  typedef typename Collapse<Quantity<Dimensions<Length, Mass, Time, Current, Temperature,
                                       Amount, LuminousIntensity, Winding,
                                       Wrapping >> >::ResultType ResultType;
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
    Wrapping          = -Right::Dimensions::Wrapping
  };
  typedef Quantity<Dimensions<Length, Mass, Time, Current, Temperature, Amount,
                              LuminousIntensity, Winding, Wrapping>> ResultType;
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
#pragma region Operators on units
template<typename Q>
inline Q operator*(Dimensionless const& left, Unit<Q> const& right) {
  return Q(left * right.value_);
}
template<typename Q_Left, typename Q_Right>
inline Unit<Product<Q_Left, Q_Right>> operator*(Unit<Q_Left> const& left,
                                                Unit<Q_Right> const& right) {
  return Unit<Product<Q_Left, Q_Right>>(left.value_ * right.value_);
}
template<typename Q_Left, typename Q_Right>
inline Unit<Quotient<Q_Left, Q_Right>> operator/(Unit<Q_Left> const& left, 
                                                 Unit<Q_Right> const& right) {
  return Unit<Quotient<Q_Left, Q_Right>>(left.value_ / right.value_);
}
#pragma endregion
#pragma region Base quantities
inline Length Metres(Dimensionless const& number) { return Length(number); }
inline Mass Kilograms(Dimensionless const& number) { return Mass(number); }
inline Time Seconds(Dimensionless const& number) { return Time(number); }
inline Current Amperes(Dimensionless const& number) { return Current(number); }
inline Temperature Kelvins(Dimensionless const& number) { return Temperature(number); }
inline Amount Moles(Dimensionless const& number) { return Amount(number); }
inline LuminousIntensity Candelas(Dimensionless const& number) {
 return LuminousIntensity(number);
}
inline Winding Cycles(Dimensionless const& number) { return Winding(number); }
inline Wrapping Globes(Dimensionless const& number) { return Wrapping(number); };
// The final semicolon is unneeded, but IntelliSense likes it.
#pragma endregion