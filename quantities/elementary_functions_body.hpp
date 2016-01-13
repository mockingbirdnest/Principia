#pragma once

#include <cmath>
#include <type_traits>

#include "quantities/si.hpp"

namespace principia {
namespace quantities {

namespace internal {

template<int n, typename Q>
struct NthRootGenerator<
    n,
    Q,
    std::enable_if_t<Q::Dimensions::Length % n == 0 &&
                     Q::Dimensions::Mass % n == 0 &&
                     Q::Dimensions::Time % n == 0 &&
                     Q::Dimensions::Current % n == 0 &&
                     Q::Dimensions::Temperature % n == 0 &&
                     Q::Dimensions::Amount % n == 0 &&
                     Q::Dimensions::LuminousIntensity % n == 0 &&
                     Q::Dimensions::Winding % n == 0 &&
                     Q::Dimensions::Angle % n == 0 &&
                     Q::Dimensions::SolidAngle % n == 0>> {
  enum {
    Length            = Q::Dimensions::Length / n,
    Mass              = Q::Dimensions::Mass / n,
    Time              = Q::Dimensions::Time / n,
    Current           = Q::Dimensions::Current / n,
    Temperature       = Q::Dimensions::Temperature / n,
    Amount            = Q::Dimensions::Amount / n,
    LuminousIntensity = Q::Dimensions::LuminousIntensity / n,
    Winding           = Q::Dimensions::Winding / n,
    Angle             = Q::Dimensions::Angle / n,
    SolidAngle        = Q::Dimensions::SolidAngle / n
  };
  using Type = Quantity<Dimensions<Length, Mass, Time, Current, Temperature,
                                   Amount, LuminousIntensity, Winding, Angle,
                                   SolidAngle>>;
};

}  // namespace internal

inline double Sqrt(double const x) {
  return std::sqrt(x);
}

template<typename D>
inline SquareRoot<Quantity<D>> Sqrt(Quantity<D> const& x) {
  return SquareRoot<Quantity<D>>(std::sqrt(x.magnitude_));
}

inline double Cbrt(double const x) {
  return std::cbrt(x);
}

template<typename D>
inline CubeRoot<Quantity<D>> Cbrt(Quantity<D> const& x) {
  return CubeRoot<Quantity<D>>(std::cbrt(x.magnitude_));
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
Angle ArcTan(Quantity<D> const& y, Quantity<D> const& x) {
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

}  // namespace quantities
}  // namespace principia
