// PhysicalQuantities.h
#include <functional>

#pragma once

using namespace System;

namespace PhysicalQuantities {
#pragma region Metaprogramming
  template <typename T> using static_not = std::integral_constant<bool, !T::value>;
#pragma endregion
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
  template<typename D> struct ScalarQuantity;
#pragma region Base quantities
  typedef ScalarQuantity<Dimensions<0, 0, 0, 0>> Dimensionless;
  typedef ScalarQuantity<Dimensions<1, 0, 0, 0>> Length;
  typedef ScalarQuantity<Dimensions<0, 1, 0, 0>> Time;
  typedef ScalarQuantity<Dimensions<0, 0, 1, 0>> Mass;
  typedef ScalarQuantity<Dimensions<0, 0, 0, 1>> Temperature;
#pragma endregion
  template<typename D>
  struct ScalarQuantity {
  public:
    ScalarQuantity() = default;
    typedef typename D Dimensions;
    template <typename Dummy = std::enable_if<std::is_same<ScalarQuantity<D>, Dimensionless>::value, double>::type>
    ScalarQuantity(double magnitude) : magnitude_(magnitude) {}
    friend double Value(Dimensionless);
    friend Length Metres(double);
    friend Time Seconds(double);
    friend Mass Kilograms(double);
    friend Temperature Kelvins(double);
    template<typename D> friend ScalarQuantity<D> operator+ (ScalarQuantity<D>);
    template<typename D> friend ScalarQuantity<D> operator- (ScalarQuantity<D>);
    template<typename D>
    friend ScalarQuantity<D> operator+ (ScalarQuantity<D>, ScalarQuantity<D>);
    template<typename D>
    friend ScalarQuantity<D> operator- (ScalarQuantity<D>, ScalarQuantity<D>);
    template<typename D_Left, typename D_Right>
    friend Product<typename ScalarQuantity<D_Left>, typename ScalarQuantity <D_Right>>
      operator *(ScalarQuantity<D_Left>, ScalarQuantity<D_Right>);
    template<typename D_Right>
    friend ScalarQuantity<D_Right> operator *(double, ScalarQuantity<D_Right>);
    template<typename D_Left>
    friend ScalarQuantity<D_Left> operator *(ScalarQuantity<D_Left>, double);
    template<typename D_Left>
    friend ScalarQuantity<D_Left> operator /(ScalarQuantity<D_Left>, double);
    template<typename D_Left, typename D_Right>
    friend Quotient<typename ScalarQuantity<D_Left>, typename ScalarQuantity <D_Right>>
      operator /(ScalarQuantity<D_Left>, ScalarQuantity<D_Right>);
  private:
    explicit ScalarQuantity(double magnitude) : magnitude_(magnitude) {};
    double   magnitude_;
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
    typedef ScalarQuantity<Dimensions<Length, Time, Mass, Temperature>> ResultType;
  };
  template<typename Left, typename Right>
  struct QuotientGenerator {
    enum {
      Length = Left::Dimensions::Length - Right::Dimensions::Length,
      Time = Left::Dimensions::Time - Right::Dimensions::Time,
      Mass = Left::Dimensions::Mass - Right::Dimensions::Mass,
      Temperature = Left::Dimensions::Temperature - Right::Dimensions::Temperature
    };
    typedef ScalarQuantity<Dimensions<Length, Time, Mass, Temperature>> ResultType;
  };
#pragma endregion
#pragma region Additive group
  template<typename D>
  inline ScalarQuantity<D> operator +(ScalarQuantity<D> right) {
    return ScalarQuantity<D>(+right.magnitude_);
  }
  template<typename D>
  inline ScalarQuantity<D> operator -(ScalarQuantity<D> right) {
    return ScalarQuantity<D>(-right.magnitude_);
  }
  template<typename D>
  inline ScalarQuantity<D> operator +(ScalarQuantity<D> left,
                                      ScalarQuantity<D> right) {
    return ScalarQuantity<D>(left.magnitude_ + right.magnitude_);
  }
  template<typename D>
  inline ScalarQuantity<D> operator -(ScalarQuantity<D> left,
                                      ScalarQuantity<D> right) {
    return ScalarQuantity<D>(left.magnitude_ + right.magnitude_);
  }
#pragma endregion
#pragma region Multiplicative group
  template<typename D_Left, typename D_Right>
  inline Product < typename ScalarQuantity<D_Left>,
                   typename ScalarQuantity <D_Right> >
    operator *(ScalarQuantity<D_Left> left, ScalarQuantity<D_Right> right) {
      return Product<typename ScalarQuantity<D_Left>, typename ScalarQuantity<D_Right>>(left.magnitude_ * right.magnitude_);
    }
  template<typename D_Left, typename D_Right>
  inline Quotient<typename ScalarQuantity<D_Left>, typename ScalarQuantity <D_Right>>
    operator /(ScalarQuantity<D_Left> left, ScalarQuantity<D_Right> right) {
      return Quotient<typename ScalarQuantity<D_Left>, typename ScalarQuantity<D_Right>>(left.magnitude_ / right.magnitude_);
    }
  template<typename D_Right>
  inline ScalarQuantity<D_Right> operator *(double left, ScalarQuantity<D_Right> right) {
    return ScalarQuantity<D_Right>(left * right.magnitude_);
  }
  template<typename D_Left>
  inline ScalarQuantity<D_Left> operator *(ScalarQuantity<D_Left> left, double right) {
    return ScalarQuantity<D_Right>(left.magnitude_ * right)
  }
  template<typename D_Left>
  inline ScalarQuantity<D_Left> operator /(ScalarQuantity<D_Left> left, double right) {
    return ScalarQuantity<D_Right>(left.magnitude_ / right)
  }
#pragma endregion
#pragma region Assigment operators
  template<typename D>
  inline void operator += (ScalarQuantity<D> left, ScalarQuantity<D> right) {
    left = left + right;
  }
  template<typename D>
  inline void operator -= (ScalarQuantity<D> left, ScalarQuantity<D> right) {
    left = left - right;
  }
  template<typename D>
  inline void operator *= (ScalarQuantity<D> left, Dimensionless right) {
    left = left * right;
  }
  template<typename D>
  inline void operator /= (ScalarQuantity<D> left, Dimensionless right) {
    left = left / right;
  }
#pragma endregion
#pragma region Dimensionless numbers
  inline double Value(Dimensionless number) {
    return number.magnitude_;
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
    Dimensionless x = 3.0;
    Momentum p = m * v;
    Energy E = .5 * m * v * v;
    Force F = 1000.0 * Newton;
    double numberOfKelvins = Value((9.8 * Kelvin - Celsius(14)) / Kelvin);
    Dimensionless N = 1e23;
    Volume V = 5 * Metre * Metre * Metre + 2 * Litre;
    Temperature T = 3 * Kelvin;
    Pressure P = N * BoltzmannConstant * T / V;
  }
}
