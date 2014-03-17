// PhysicalQuantities.h

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
  template<typename Left, typename Right>  struct QuotientGenerator;
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
    friend double Value(DimensionlessNumber);
    friend DimensionlessNumber Dimensionless(double);
    friend Length Metres(double);
    friend Time Seconds(double);
    friend Mass Kilograms(double);
    friend Temperature Kelvins(double);
    template<typename D> friend Quantity<D> operator+ (Quantity<D> right);
    template<typename D> friend Quantity<D> operator- (Quantity<D> right);
    template<typename D>
    friend Quantity<D> operator+ (Quantity<D> left, Quantity<D> right);
    template<typename D>
    friend Quantity<D> operator- (Quantity<D> left, Quantity<D> right);
    template<typename D_Left, typename D_Right>
    friend Product<typename Quantity<D_Left>, typename Quantity <D_Right>>
      operator *(Quantity<D_Left> left, Quantity<D_Right> right);
    template<typename D_Left, typename D_Right>
    friend Quotient<typename Quantity<D_Left>, typename Quantity <D_Right>>
      operator /(Quantity<D_Left> left, Quantity<D_Right> right);
  private:
    explicit Quantity(double m) : _magnitude(m) {};
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
  inline DimensionlessNumber Dimensionless(double value) {
    return DimensionlessNumber(value);
  }
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
  inline Speed MetresPerSecond(double number) {
    return Metres(number) / Seconds(1.0);
  }
  inline Acceleration MetresPerSquaredSecond(double number) {
    return MetresPerSecond(number) / Seconds(1.0);
  }
  const Force Newton = Metre * Kilogram / (Second * Second);
  inline Force Newtons(double number) {
    return MetresPerSquaredSecond(number) * Kilograms(1.0);
  }
  const Energy Joule = Newton * Metre;
  inline Energy Joules(double number) { return Dimensionless(number) * Joule; }
#pragma endregion
#pragma region Thermodynamics
  inline Pressure Pascals(double number) {
    return Newtons(number) / (Metre * Metre);
  }
  const Pressure Pascal = Pascals(1.0);
  const Volume Litres = Dimensionless(1e-3) * Metre * Metre * Metre;
#pragma endregion
#pragma endregion
#pragma region Constants
  const Entropy BoltzmannConstant = Joules(1.3806488e-23) / Kelvin;
#pragma endregion
  public ref class foo {
    double test() {
      Mass m = Kilograms(5.0);
      Speed v = Metres(1.2) / Seconds(4.0) - MetresPerSecond(3.14);
      v += Metres(42.0) / Second;
      v *= Dimensionless(5.2);
      v /= Dimensionless(0.7);
      Momentum p = m * v;
      Energy E = Dimensionless(.5) * m * v * v;
      Force F = Newtons(1000.0);
      double numberOfKelvins = Value((Kelvins(3) - Celsius(5)) / Kelvin);
      DimensionlessNumber N = Dimensionless(1e23);
      Volume V = Metres(3) * Metre * Metre;
      Temperature T = Kelvins(3);
      Pressure P = N * BoltzmannConstant * T / V;
      return Value(P / Pascal);
    }
  };
}
