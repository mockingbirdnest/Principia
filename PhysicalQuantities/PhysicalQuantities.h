// PhysicalQuantities.h

#pragma once

namespace PhysicalQuantities {
template<int LengthExponent, int TimeExponent, int MassExponent,
         int TemperatureExponent>
struct Dimensions;
template<typename Left, typename Right> struct ProductGenerator;
template<typename Left, typename Right> struct QuotientGenerator;
template<typename Left, typename Right> 
using Quotient = typename QuotientGenerator<Left, Right>::ResultType;
template<typename Left, typename Right> 
using Product = typename ProductGenerator<Left, Right>::ResultType;
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
  friend double              Value(DimensionlessScalar const);
  friend DimensionlessScalar Dimensionless(double const);
  friend Length              Metres(double const);
  friend Time                Seconds(double const);
  friend Mass                Kilograms(double const);
  friend Temperature         Kelvins(double const);
  template<typename D> friend Quantity<D> operator+ (Quantity<D> const);
  template<typename D> friend Quantity<D> operator- (Quantity<D> const);
  template<typename D> friend Quantity<D> operator+ (Quantity<D> const, 
                                                     Quantity<D> const);
  template<typename D> friend Quantity<D> operator- (Quantity<D> const, 
                                                     Quantity<D> const);
  template<typename D_Left, typename D_Right>
  friend Product<typename Quantity<D_Left>, typename Quantity <D_Right>> 
    operator *(Quantity<D_Left> const, Quantity<D_Right> const);
  template<typename D_Left, typename D_Right>
  friend Quotient<typename Quantity<D_Left>, typename Quantity <D_Right>> 
    operator /(Quantity<D_Left> const, Quantity<D_Right> const);
private:
  explicit Quantity(double const magnitude) : magnitude_(magnitude) {};
  double   magnitude_;
};
template<typename D>
void operator += (Quantity<D>, Quantity<D> const);
template<typename D>
inline void operator -= (Quantity<D>, Quantity<D> const);
template<typename D>
inline void operator *= (Quantity<D>, DimensionlessScalar const);
template<typename D>
inline void operator /= (Quantity<D>, DimensionlessScalar const);

template<typename Q>
struct Unit {
public:
  explicit Unit(Q const value) : value_(value) {};
  template<typename Q> friend Q operator*(double const, Unit<Q> const);
  template<typename Q_Left, typename Q_Right>
  friend Unit<Product<Q_Left, Q_Right>> operator*(Unit<Q_Left> const,
                                                  Unit<Q_Right> const);
  template<typename Q_Left, typename Q_Right>
  friend Unit<Quotient<Q_Left, Q_Right>> operator/(Unit<Q_Left> const,
                                                   Unit<Q_Right> const);
private:
  Q const value_;
};
#include "PhysicalQuantities-body-inl.h"
}
