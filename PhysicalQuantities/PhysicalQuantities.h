// PhysicalQuantities.h
#include <functional>

#pragma once

using namespace System;

namespace PhysicalQuantities {
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
  template<typename Left, typename Right> struct ProductGenerator;
  template<typename Left, typename Right> struct QuotientGenerator;
  template<typename Left, typename Right>
  using Quotient = typename QuotientGenerator<Left, Right>::ResultType;
  template<typename Left, typename Right>
  using Product = typename ProductGenerator<Left, Right>::ResultType;
  template<typename D> class Quantity;
#pragma region Base quantities
  typedef Quantity<Dimensions<0, 0, 0, 0>> DimensionlessNumber;
  typedef Quantity<Dimensions<1, 0, 0, 0>> Length;
  typedef Quantity<Dimensions<0, 1, 0, 0>> Time;
  typedef Quantity<Dimensions<0, 0, 1, 0>> Mass;
  typedef Quantity<Dimensions<0, 0, 0, 1>> Temperature;
#pragma endregion
  template<typename D>
  class Quantity {
  public:
    Quantity() {};
    typedef typename D Dimensions;
    template <typename Dummy = std::enable_if<std::is_same<Quantity<D>, DimensionlessNumber>::value, double>::type>
    Quantity(double magnitude) : _magnitude(magnitude) {}
    friend double Value(DimensionlessNumber);
    friend Length Metres(double);
    friend Time Seconds(double);
    friend Mass Kilograms(double);
    friend Temperature Kelvins(double);
    template<typename D> friend Quantity<D> operator+ (Quantity<D>);
    template<typename D> friend Quantity<D> operator- (Quantity<D>);
    template<typename D>
    friend Quantity<D> operator+ (Quantity<D>, Quantity<D>);
    template<typename D>
    friend Quantity<D> operator- (Quantity<D>, Quantity<D>);
    template<typename D_Left, typename D_Right>
    friend Product<typename Quantity<D_Left>, typename Quantity <D_Right>>
      operator *(Quantity<D_Left>, Quantity<D_Right>);
    template<typename D_Right>
    friend Quantity<D_Right> operator *(double, Quantity<D_Right>);
    template<typename D_Left>
    friend Quantity<D_Left> operator *(Quantity<D_Left>, double);
    template<typename D_Left>
    friend Quantity<D_Left> operator /(Quantity<D_Left>, double);
    template<typename D_Left, typename D_Right>
    friend Quotient<typename Quantity<D_Left>, typename Quantity <D_Right>>
      operator /(Quantity<D_Left>, Quantity<D_Right>);
  private:
    explicit Quantity(double magnitude) : _magnitude(magnitude) {};
    double   _magnitude;
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
    return Quantity<D>(+right._magnitude);
  }
  template<typename D>
  inline Quantity<D> operator -(Quantity<D> right) {
    return Quantity<D>(-right._magnitude);
  }
  template<typename D>
  inline Quantity<D> operator +(Quantity<D> left, Quantity<D> right) {
    return Quantity<D>(left._magnitude + right._magnitude);
  }
  template<typename D>
  inline Quantity<D> operator -(Quantity<D> left, Quantity<D> right) {
    return Quantity<D>(left._magnitude + right._magnitude);
  }
#pragma endregion
#pragma region Multiplicative group
  template<typename D_Left, typename D_Right>
  inline Product<typename Quantity<D_Left>, typename Quantity <D_Right>>
    operator *(Quantity<D_Left> left, Quantity<D_Right> right) {
      return Product<typename Quantity<D_Left>, typename Quantity<D_Right>>(left._magnitude * right._magnitude);
    }
  template<typename D_Left, typename D_Right>
  inline Quotient<typename Quantity<D_Left>, typename Quantity <D_Right>>
    operator /(Quantity<D_Left> left, Quantity<D_Right> right) {
      return Quotient<typename Quantity<D_Left>, typename Quantity<D_Right>>(left._magnitude / right._magnitude);
    }
  template<typename D_Right>
  inline Quantity<D_Right> operator *(double left, Quantity<D_Right> right) {
    return Quantity<D_Right>(left * right._magnitude);
  }
  template<typename D_Left>
  inline Quantity<D_Left> operator *(Quantity<D_Left> left, double right) {
    return Quantity<D_Right>(left._magnitude * right)
  }
  template<typename D_Left>
  inline Quantity<D_Left> operator /(Quantity<D_Left> left, double right) {
    return Quantity<D_Right>(left._magnitude / right)
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
  inline void operator *= (Quantity<D> left, DimensionlessNumber right) {
    left = left * right;
  }
  template<typename D>
  inline void operator /= (Quantity<D> left, DimensionlessNumber right) {
    left = left / right;
  }
#pragma endregion
#pragma region Dimensionless numbers
  inline double Value(DimensionlessNumber number) {
    return number._magnitude;
  }
#pragma endregion
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
#pragma region Common Units
#pragma region Base SI units
  inline Length Metres(double number) { return Length(number); }
  inline Time Seconds(double number) { return Time(number); }
  inline Mass Kilograms(double number) { return Mass(number); }
  inline Temperature Kelvins(double number) { return Temperature(number); }
  const Length Metre = Metres(1.0);
  const Time Second = Seconds(1.0);
  const Mass Kilogram = Kilograms(1.0);
  const Temperature Kelvin = Kelvins(1.0);
#pragma endregion
#pragma region Further units for base quantities
  inline Temperature Celsius(double number) {
    return Kelvins(number) + Kelvins(273.15);
  }
#pragma endregion
#pragma region General mechanics
  const Force Newton = Metre * Kilogram / (Second * Second);
  const Energy Joule = Newton * Metre;
#pragma endregion
#pragma region Thermodynamics
  const Pressure Pascal = Newton / (Metre * Metre);
  const Volume Litre = 1e-3 * Metre * Metre * Metre;
#pragma endregion
#pragma endregion
#pragma region Constants
  const Entropy BoltzmannConstant = 1.3806488e-23 * Joule / Kelvin;
#pragma endregion
  void test() {
    Mass m = Kilograms(5.0);
    Speed v = 1.2 * Metre / Second;
    v += 43 * Metre / Second;
    v *= 5.2;
    v /= 0.7;
    DimensionlessNumber x = 3.0;
    Momentum p = m * v;
    Energy E = .5 * m * v * v;
    Force F = 1000.0 * Newton;
    double numberOfKelvins = Value((9.8 * Kelvin - Celsius(14)) / Kelvin);
    DimensionlessNumber N = 1e23;
    Volume V = 5 * Metre * Metre * Metre + Litre;
    Temperature T = 3 * Kelvin;
    Pressure P = N * BoltzmannConstant * T / V;
  }
}
