// PhysicalQuantities-body-inl.h

#pragma once

template<int LengthExponent,
         int TimeExponent,
         int MassExponent,
         int TemperatureExponent>
struct Dimensions {
  enum {
    Length = LengthExponent,
    Time = TimeExponent,
    Mass = MassExponent,
    Temperature = TemperatureExponent
  };
};
#pragma region Type generators
template<typename Left, typename Right>
struct ProductGenerator {
  enum {
    Length = Left::Dimensions::Length + Right::Dimensions::Length,
    Time = Left::Dimensions::Time + Right::Dimensions::Time,
    Mass = Left::Dimensions::Mass + Right::Dimensions::Mass,
    Temperature = Left::Dimensions::Temperature + Right::Dimensions::Temperature
  };
  typedef Quantity<Dimensions<Length, Time, Mass, Temperature>> ResultType;
};
template<typename Left, typename Right>
struct QuotientGenerator {
  enum {
    Length = Left::Dimensions::Length - Right::Dimensions::Length,
    Time = Left::Dimensions::Time - Right::Dimensions::Time,
    Mass = Left::Dimensions::Mass - Right::Dimensions::Mass,
    Temperature = Left::Dimensions::Temperature - Right::Dimensions::Temperature
  };
  typedef Quantity<Dimensions<Length, Time, Mass, Temperature>> ResultType;
};
#pragma endregion
#pragma region Additive group
template<typename D>
inline Quantity<D> operator +(Quantity<D> right) {
  return Quantity<D>(+right.magnitude_);
}
template<typename D>
inline Quantity<D> operator -(Quantity<D> right) {
  return Quantity<D>(-right.magnitude_);
}
template<typename D>
inline Quantity<D> operator +(Quantity<D> left, Quantity<D> right) {
  return Quantity<D>(left.magnitude_ + right.magnitude_);
}
template<typename D>
inline Quantity<D> operator -(Quantity<D> left, Quantity<D> right) {
  return Quantity<D>(left.magnitude_ + right.magnitude_);
}
#pragma endregion
#pragma region Multiplicative group
template<typename D_Left, typename D_Right>
inline Product <typename Quantity<D_Left>, typename Quantity <D_Right>>
operator *(Quantity<D_Left> left, Quantity<D_Right> right) {
  return Product<typename Quantity<D_Left>, typename Quantity<D_Right>>(left.magnitude_ * right.magnitude_);
}
template<typename D_Left, typename D_Right>
inline Quotient<typename Quantity<D_Left>, typename Quantity <D_Right>> operator /(Quantity<D_Left> left,
                                                                                   Quantity<D_Right> right) {
  return Quotient<typename Quantity<D_Left>, typename Quantity<D_Right>>(left.magnitude_ / right.magnitude_);
}
#pragma endregion
#pragma region Assigment operators
template<typename D>
inline void operator += (Quantity<D> left, Quantity<D> right) {
  left = left + right;
}
template<typename D>
inline void operator -= (Quantity<D> left, Quantity<D> right) {
  left = left - right;
}
template<typename D>
inline void operator *= (Quantity<D> left, DimensionlessScalar right) {
  left = left * right;
}
template<typename D>
inline void operator /= (Quantity<D> left, DimensionlessScalar right) {
  left = left / right;
}
#pragma endregion
#pragma region Dimensionless scalars
inline DimensionlessScalar Dimensionless(double value) {
  return DimensionlessScalar(value);
}
inline double Value(DimensionlessScalar number) {
  return number.magnitude_;
};
#pragma endregion