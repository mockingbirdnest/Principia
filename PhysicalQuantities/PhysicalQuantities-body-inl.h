// PhysicalQuantities-body-inl.h

#pragma once

template<int LengthExponent, int MassExponent, int TimeExponent,
         int CurrentExponent, int TemperatureExponent, int AmountExponent,
         int LuminousIntensityExponent, int WindingExponent>
struct Dimensions {
  enum {
    Length            = LengthExponent,
    Mass              = MassExponent,
    Time              = TimeExponent,
    Current           = CurrentExponent,
    Temperature       = TemperatureExponent,
    Amount            = AmountExponent,
    LuminousIntensity = LuminousIntensityExponent,
    Winding           = WindingExponent
  };
};
#pragma region Type generators
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
    Winding           = Left::Dimensions::Winding + Right::Dimensions::Winding
  };
  typedef Quantity<Dimensions<Length, Mass, Time, Current, Temperature, Amount,
                              LuminousIntensity, Winding>> ResultType;
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
    Winding           = Left::Dimensions::Winding - Right::Dimensions::Winding
  };
  typedef Quantity<Dimensions<Length, Mass, Time, Current, Temperature, Amount,
                              LuminousIntensity, Winding>> ResultType;
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
template<typename D_Left, typename D_Right>
inline Product <typename Quantity<D_Left>, typename Quantity <D_Right>>
operator*(Quantity<D_Left> const& left,
          Quantity<D_Right> const& right) {
  return Product<typename Quantity<D_Left>, 
                 typename Quantity<D_Right>>(left.magnitude_ * 
                                             right.magnitude_);
}
template<typename D_Left, typename D_Right>
inline Quotient<typename Quantity<D_Left>, typename Quantity <D_Right>> 
operator/(Quantity<D_Left> const& left,
          Quantity<D_Right> const& right) {
  return Quotient<typename Quantity<D_Left>,
                  typename Quantity<D_Right>>(left.magnitude_ /
                                              right.magnitude_);
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
inline void operator*=(Quantity<D>& left, DimensionlessScalar const& right) {
  left = left * right;
}
template<typename D>
inline void operator/=(Quantity<D>& left, DimensionlessScalar const& right) {
  left = left / right;
}
#pragma endregion
#pragma region Operators on units
template<typename Q>
inline Q operator*(double const left, Unit<Q> const& right) {
  return Q(Dimensionless(left) * right.value_);
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
inline DimensionlessScalar Dimensionless(double const value) {
  return DimensionlessScalar(value);
}
inline double Value(DimensionlessScalar const &number) {
  return number.magnitude_;
}
inline Length Metres(double const number) { return Length(number); }
inline Mass Kilograms(double const number) { return Mass(number); }
inline Time Seconds(double const number) { return Time(number); }
inline Current Amperes(double const number) { return Current(number); }
inline Temperature Kelvins(double const number) { return Temperature(number); }
inline Amount Moles(double const number) { return Amount(number); }
inline LuminousIntensity Candelas(double const number) {
 return LuminousIntensity(number);
}
inline Winding Cycles(double const number) { return Winding(number); };
// The final semicolon is unneeded, but IntelliSense likes it.
#pragma endregion