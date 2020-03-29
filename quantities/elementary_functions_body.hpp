
#pragma once

#include "quantities/elementary_functions.hpp"

#include <pmmintrin.h>

#include <cmath>
#include <type_traits>

#include "quantities/si.hpp"
#include "numerics/cbrt.hpp"

namespace principia {
namespace quantities {
namespace internal_elementary_functions {

using si::Radian;

template<typename Q1, typename Q2>
Product<Q1, Q2> FusedMultiplyAdd(Q1 const& x,
                                 Q2 const& y,
                                 Product<Q1, Q2> const& z) {
  return si::Unit<Product<Q1, Q2>> * std::fma(x / si::Unit<Q1>,
                                              y / si::Unit<Q2>,
                                              z / si::Unit<Product<Q1, Q2>>);
}

template<typename Q>
FORCE_INLINE(inline) Q Abs(Q const& quantity) {
  return si::Unit<Q> * std::abs(quantity / si::Unit<Q>);
}

template<typename Q>
Q Mod(Q const& argument, Q const& modulus) {
  double const result =
      std::fmod(argument / si::Unit<Q>, modulus / si::Unit<Q>);
  if (result > 0.0) {
    return result * si::Unit<Q>;
  } else {
    return result * si::Unit<Q> + modulus;
  }
}

template<typename Q>
SquareRoot<Q> Sqrt(Q const& x) {
#if PRINCIPIA_USE_SSE3_INTRINSICS
  auto const x_128d = _mm_set_sd(x / si::Unit<Q>);
  return si::Unit<SquareRoot<Q>> * _mm_cvtsd_f64(_mm_sqrt_sd(x_128d, x_128d));
#else
  return si::Unit<SquareRoot<Q>> * std::sqrt(x / si::Unit<Q>);
#endif
}

template<typename Q>
CubeRoot<Q> Cbrt(Q const& x) {
  return si::Unit<CubeRoot<Q>> * numerics::Cbrt(x / si::Unit<Q>);
}

template<int exponent>
constexpr double Pow(double x) {
  return std::pow(x, exponent);
}

// Static specializations for frequently-used exponents, so that this gets
// turned into multiplications at compile time.

template<>
inline constexpr double Pow<-3>(double x) {
  return 1 / (x * x * x);
}

template<>
inline constexpr double Pow<-2>(double x) {
  return 1 / (x * x);
}

template<>
inline constexpr double Pow<-1>(double x) {
  return 1 / x;
}

template<>
inline constexpr double Pow<0>(double x) {
  return 1;
}

template<>
inline constexpr double Pow<1>(double x) {
  return x;
}

template<>
inline constexpr double Pow<2>(double x) {
  return x * x;
}

template<>
inline constexpr double Pow<3>(double x) {
  return x * x * x;
}

template<int exponent, typename Q>
constexpr Exponentiation<Q, exponent> Pow(Q const& x) {
  return si::Unit<Exponentiation<Q, exponent>> * Pow<exponent>(x / si::Unit<Q>);
}

inline double Sin(Angle const& α) {
  return std::sin(α / Radian);
}

inline double Cos(Angle const& α) {
  return std::cos(α / Radian);
}

inline double Tan(Angle const& α) {
  return std::tan(α / Radian);
}

inline Angle ArcSin(double const x) {
  return std::asin(x) * Radian;
}
inline Angle ArcCos(double const x) {
  return std::acos(x) * Radian;
}
inline Angle ArcTan(double const x) {
  return std::atan(x) * Radian;
}
inline Angle ArcTan(double const y, double const x) {
  return std::atan2(y, x) * Radian;
}
template<typename D>
Angle ArcTan(Quantity<D> const& y, Quantity<D> const& x) {
  return ArcTan(y / si::Unit<Quantity<D>>, x / si::Unit<Quantity<D>>);
}

inline double Sinh(Angle const& α) {
  return std::sinh((α / Radian));
}
inline double Cosh(Angle const& α) {
  return std::cosh((α / Radian));
}
inline double Tanh(Angle const& α) {
  return std::tanh((α / Radian));
}

inline Angle ArcSinh(double const x) {
  return std::asinh(x) * Radian;
}
inline Angle ArcCosh(double const x) {
  return std::acosh(x) * Radian;
}
inline Angle ArcTanh(double const x) {
  return std::atanh(x) * Radian;
}

inline Angle UnwindFrom(Angle const& previous_angle, Angle const& α) {
  return α + std::nearbyint((previous_angle - α) / (2 * π * Radian)) *
                 (2 * π * Radian);
}

}  // namespace internal_elementary_functions
}  // namespace quantities
}  // namespace principia
