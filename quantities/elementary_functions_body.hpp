
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

template<typename Q1, typename Q2>
Product<Q1, Q2> FusedMultiplyAdd(Q1 const& x,
                                 Q2 const& y,
                                 Product<Q1, Q2> const& z) {
  return SIUnit<Product<Q1, Q2>>() * std::fma(x / SIUnit<Q1>(),
                                              y / SIUnit<Q2>(),
                                              z / SIUnit<Product<Q1, Q2>>());
}

template<typename Q>
FORCE_INLINE(inline) Q Abs(Q const& quantity) {
  return SIUnit<Q>() * std::abs(quantity / SIUnit<Q>());
}

template<typename Q>
Q Mod(Q const& argument, Q const& modulus) {
  double const result =
      std::fmod(argument / SIUnit<Q>(), modulus / SIUnit<Q>());
  if (result > 0.0) {
    return result * SIUnit<Q>();
  } else {
    return result * SIUnit<Q>() + modulus;
  }
}

template<typename Q>
SquareRoot<Q> Sqrt(Q const& x) {
#if PRINCIPIA_USE_SSE3_INTRINSICS
  auto const x_128d = _mm_set_sd(x / SIUnit<Q>());
  return SIUnit<SquareRoot<Q>>() * _mm_cvtsd_f64(_mm_sqrt_sd(x_128d, x_128d));
#else
  return SIUnit<SquareRoot<Q>>() * std::sqrt(x / SIUnit<Q>());
#endif
}

template<typename Q>
CubeRoot<Q> Cbrt(Q const& x) {
  return SIUnit<CubeRoot<Q>>() * numerics::Cbrt(x / SIUnit<Q>());
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
  return SIUnit<Exponentiation<Q, exponent>>() * Pow<exponent>(x / SIUnit<Q>());
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
  return ArcTan(y / SIUnit<Quantity<D>>(), x / SIUnit<Quantity<D>>());
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

}  // namespace internal_elementary_functions
}  // namespace quantities
}  // namespace principia
