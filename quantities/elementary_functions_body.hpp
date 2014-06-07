#pragma once

#include <math.h>

#include "quantities/si.hpp"

namespace principia {
namespace quantities {

namespace type_generators {
template<typename Q>
struct SquareRootGenerator<
    Q, Condition<!(Q::Dimensions::Length & 1 ||
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
}  // namespace type_generators



inline double Sqrt(double const x) {
  return std::sqrt(x);
}

inline double Sin(Angle const& α) {
  return std::sin(α / si::Radian);
}
inline double Cos(Angle const& α) {
  return std::cos(α / si::Radian);
}
inline double Tan(Angle const& α) {
  return std::tan(α / si::Radian);
}

inline Angle ArcSin(double const x) {
  return std::asin(x) * si::Radian;
}
inline Angle ArcCos(double const x) {
  return std::acos(x) * si::Radian;
}
inline Angle ArcTan(double const y, double const x) {
  return std::atan2(y, x) * si::Radian;
}
template<typename D>
inline Angle ArcTan(Quantity<D> const& y, Quantity<D> const& x) {
  return ArcTan(y.magnitude_, x.magnitude_);
}

inline double Sinh(Angle const& α) {
  return std::sinh((α / si::Radian));
}
inline double Cosh(Angle const& α) {
  return std::cosh((α / si::Radian));
}
inline double Tanh(Angle const& α) {
  return std::tanh((α / si::Radian));
}

inline Angle ArcSinh(double const x) {
  return std::asinh(x) * si::Radian;
}
inline Angle ArcCosh(double const x) {
  return std::acosh(x) * si::Radian;
}
inline Angle ArcTanh(double const x) {
  return std::atanh(x) * si::Radian;
}

template<typename D>
inline SquareRoot<Quantity<D>> Sqrt(Quantity<D> const& x) {
  return SquareRoot<Quantity<D>>(std::sqrt(x.magnitude_));
}

}  // namespace quantities
}  // namespace principia
