// PhysicalQuantities.h

#pragma once

using namespace System;

namespace PhysicalQuantities {
  template<int LengthExponent, int TimeExponent, int MassExponent, int TemperatureExponent>
  struct Dimensions {
    enum { Length = LengthExponent, Time = TimeExponent, Mass = MassExponent, Temperature = TemperatureExponent };
  };
  template<typename Left, typename Right> struct ProductGenerator;
  template<typename Left, typename Right> struct QuotientGenerator;
  template<typename Left, typename Right> using Quotient = typename QuotientGenerator<Left, Right>::ResultType;
  template<typename Left, typename Right> using Product = typename ProductGenerator<Left, Right>::ResultType;
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
  inline Product <typename Quantity<D_Left>, typename Quantity <D_Right>> operator *(Quantity<D_Left> left,
                                                                                     Quantity<D_Right> right) {
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
  }
#pragma endregion
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
  const Unit<Length> Metre = Unit<Length>(Metres(1.0));
  const Unit<Time> Second = Unit<Time>(Seconds(1.0));
  const Unit<Mass> Kilogram = Unit<Mass>(Kilograms(1.0));
  const Unit<Temperature> Kelvin = Unit<Temperature>(Kelvins(1.0));
#pragma endregion
#pragma region Further units for base quantities
  inline Temperature Celsius(double number) {
    return Kelvins(number) + Kelvins(273.15);
  }
#pragma endregion
#pragma region General mechanics
  const Unit<Force> Newton = Metre * Kilogram / (Second * Second);
  const Unit<Energy> Joule = Newton * Metre;
#pragma endregion
#pragma region Thermodynamics
  const Unit<Pressure> Pascal = Newton / (Metre * Metre);
  const Unit<Volume> Litre = Unit<Volume>(1e-3 * (Metre * Metre * Metre));
#pragma endregion
#pragma endregion
#pragma region Constants
  const Entropy BoltzmannConstant = 1.3806488e-23 * (Joule / Kelvin);
#pragma endregion
  void test() {
    Mass m = 5.0 * Kilogram;
    Speed v = 1.2 * (Metre / Second);
    v += 43 * (Metre / Second);
    v *= Dimensionless(5.2);
    v /= Dimensionless(0.7);
    DimensionlessScalar x = Dimensionless(3.0);
    Momentum p = m * v;
    Energy E = Dimensionless(.5) * m * v * v;
    Force F = 1000.0 * Newton;
    double numberOfKelvins = Value((9.8 * Kelvin - Celsius(14)) / (1.0 * Kelvin));
    DimensionlessScalar N = Dimensionless(1e23);
    Volume V = 5 * (Metre * Metre * Metre) + 2 * Litre;
    Temperature T = 3 * Kelvin;
    Pressure P = N * BoltzmannConstant * T / V;
  }
}
