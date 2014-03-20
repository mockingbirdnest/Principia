// PhysicalQuantities.h

#pragma once

using namespace System;

namespace PhysicalQuantities {
template<int LengthExponent,
          int TimeExponent,
          int MassExponent,
          int TemperatureExponent>
struct Dimensions;
template<typename Left, typename Right> struct ProductGenerator;
template<typename Left, typename Right> struct QuotientGenerator;
template<typename Left, typename Right> using Quotient
  = typename QuotientGenerator<Left, Right>::ResultType;
template<typename Left, typename Right> using Product
  = typename ProductGenerator<Left, Right>::ResultType;
template<typename D> struct Quantity;
template<typename Q> struct Unit;
#pragma region Base quantities
typedef Quantity<Dimensions<0, 0, 0, 0>> DimensionlessScalar;
typedef Quantity<Dimensions<1, 0, 0, 0>> Length;
typedef Quantity<Dimensions<0, 1, 0, 0>> Time;
typedef Quantity<Dimensions<0, 0, 1, 0>> Mass;
typedef Quantity<Dimensions<0, 0, 0, 1>> Temperature;
#pragma endregion
template<typename D>
struct Quantity {
public:
  Quantity() = default;
  typedef typename D Dimensions;
  friend double Value(DimensionlessScalar);
  friend DimensionlessScalar Dimensionless(double);
  friend Length Metres(double);
  friend Time Seconds(double);
  friend Mass Kilograms(double);
  friend Temperature Kelvins(double);
  template<typename D> friend Quantity<D> operator+ (Quantity<D>);
  template<typename D> friend Quantity<D> operator- (Quantity<D>);
  template<typename D> friend Quantity<D> operator+ (Quantity<D>, Quantity<D>);
  template<typename D> friend Quantity<D> operator- (Quantity<D>, Quantity<D>);
  template<typename D_Left, typename D_Right>
  friend Product<typename Quantity<D_Left>, typename Quantity <D_Right>> operator *(Quantity<D_Left>,
                                                                                    Quantity<D_Right>);
  template<typename D_Left, typename D_Right>
  friend Quotient<typename Quantity<D_Left>, typename Quantity <D_Right>> operator /(Quantity<D_Left>,
                                                                                      Quantity<D_Right>);
private:
  explicit Quantity(double magnitude) : magnitude_(magnitude) {};
  double   magnitude_;
};
template<typename Q>
struct Unit {
public:
  explicit Unit(Q value) : value_(value) {};
  template<typename Q> friend Q operator*(double, Unit<Q>);
  template<typename Q_Left, typename Q_Right>
  friend Unit<Product<Q_Left, Q_Right>> operator*(Unit<Q_Left>, Unit<Q_Right>);
  template<typename Q_Left, typename Q_Right>
  friend Unit<Quotient<Q_Left, Q_Right>> operator/(Unit<Q_Left>, Unit<Q_Right>);
private:
  Q value_;
};
template<typename Q>
inline Q operator*(double left, Unit<Q> right) {
  return Q(Dimensionless(left) * right.value_);
}
template<typename Q_Left, typename Q_Right>
inline Unit<Product<Q_Left, Q_Right>> operator*(Unit<Q_Left> left, Unit<Q_Right> right) {
  return Unit<Product<Q_Left, Q_Right>>(left.value_ * right.value_);
}
template<typename Q_Left, typename Q_Right>
inline Unit<Quotient<Q_Left, Q_Right>> operator/(Unit<Q_Left> left, Unit<Q_Right> right) {
  return Unit<Quotient<Q_Left, Q_Right>>(left.value_ / right.value_);
}
#include "PhysicalQuantities-body-inl.h"
#pragma region Common quantities
#pragma region General mechanics
typedef Quotient<Length, Time> Speed;
typedef Quotient<Speed, Time> Acceleration;
typedef Product<Mass, Speed> Momentum;
typedef Quotient<Momentum, Time> Force;
typedef Product<Force, Length> Energy;
#pragma endregion
#pragma region Thermodynamics
typedef Product<Length, Length> Surface;
typedef Product<Surface, Length> Volume;
typedef Quotient<Force, Surface> Pressure;
typedef Quotient<Energy, Temperature> Entropy;
#pragma endregion
#pragma endregion
}
