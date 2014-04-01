#pragma once

#include "SI.hpp"

namespace Principia {
namespace Quantities {

namespace TypeGenerators {
  template<typename Q>
struct SquareRootGenerator<
    Q, Condition<! (Q::Dimensions::Length & 1 ||
                    Q::Dimensions::Mass & 1 ||
                    Q::Dimensions::Time & 1 ||
                    Q::Dimensions::Current & 1 ||
                    Q::Dimensions::Temperature & 1 ||
                    Q::Dimensions::Amount & 1 ||
                    Q::Dimensions::LuminousIntensity & 1 ||
                    Q::Dimensions::Winding & 1 ||
                    Q::Dimensions::Angle & 1 ||
                    Q::Dimensions::SolidAngle & 1)>> {
  enum {
    Length            = Q::Dimensions::Length / 2,
    Mass              = Q::Dimensions::Mass / 2,
    Time              = Q::Dimensions::Time / 2,
    Current           = Q::Dimensions::Current / 2,
    Temperature       = Q::Dimensions::Temperature / 2,
    Amount            = Q::Dimensions::Amount / 2,
    LuminousIntensity = Q::Dimensions::LuminousIntensity / 2,
    Winding           = Q::Dimensions::Winding / 2,
    Angle             = Q::Dimensions::Angle / 2,
    SolidAngle        = Q::Dimensions::SolidAngle / 2
  };
  typedef Quantity<
      Dimensions<Length, Mass, Time, Current, Temperature, Amount,
                 LuminousIntensity, Winding, Angle, SolidAngle>> ResultType;
};
}



inline Dimensionless Sqrt(Dimensionless const& x) {
  return std::sqrt(x.Value());
}
inline Dimensionless Log(Dimensionless const& x) {
  return std::log(x.Value());
}
inline Dimensionless Log(Dimensionless const& base, Dimensionless const& x) {
  if (base == 2) { return std::log2(x.Value()); }
  else if (base == 10) { return std::log10(x.Value()); }
  else { return Log(x) / Log(base); }
}
inline Dimensionless Exp(Dimensionless const& x) {
  return std::exp(x.Value());
}

inline Dimensionless Sin(Angle const& α) {
  return std::sin((α / SI::Radian).Value());
}
inline Dimensionless Cos(Angle const& α) {
  return std::cos((α / SI::Radian).Value());
}
inline Dimensionless Tan(Angle const& α) {
  return std::tan((α / SI::Radian).Value());
}
inline Dimensionless Cot(Angle const& α) {
  return Tan(π / 2 * SI::Radian - α);
}

inline Angle ArcSin(Dimensionless const& x) {
  return std::asin(x.Value()) * SI::Radian;
}
inline Angle ArcCos(Dimensionless const& x) {
  return std::acos(x.Value()) * SI::Radian;
}
inline Angle ArcTan(Dimensionless const& y, Dimensionless const& x) {
  return std::atan2(y.Value(), x.Value()) * SI::Radian;
}
inline Angle ArcCot(Dimensionless const& x, Dimensionless const& y) {
  return ArcTan(y, x);
}

inline Dimensionless Sinh(Angle const& α) {
  return std::sinh((α / SI::Radian).Value());
}
inline Dimensionless Cosh(Angle const& α) {
  return std::cosh((α / SI::Radian).Value());
}
inline Dimensionless Tanh(Angle const& α) {
  return std::tanh((α / SI::Radian).Value());
}
inline Dimensionless Coth(Angle const& α) {
  return 1 / Tanh(α);
}

inline Angle ArcSinh(Dimensionless const& x) {
  return std::asinh(x.Value()) * SI::Radian;
}
inline Angle ArcCosh(Dimensionless const& x) {
  return std::acosh(x.Value()) * SI::Radian;
}
inline Angle ArcTanh(Dimensionless const& x) {
  return std::atanh(x.Value()) * SI::Radian;
}
inline Angle ArcCoth(Dimensionless const& x) {
  return std::atanh(1 / x.Value()) * SI::Radian;
}

template<typename D>
inline SquareRoot<Quantity<D>> Sqrt(Quantity<D> const& x) {
  return SquareRoot<Quantity<D>>(Sqrt(x.magnitude_));
}

}
}
