#pragma once

#include "numerics/elementary_functions.hpp"

#include <immintrin.h>

#include <cmath>

#include "boost/multiprecision/cpp_bin_float.hpp"
#include "glog/logging.h"
#include "numerics/cbrt.hpp"
#include "numerics/fma.hpp"
#include "numerics/next.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace _elementary_functions {
namespace internal {

using namespace principia::numerics::_cbrt;
using namespace principia::numerics::_fma;
using namespace principia::numerics::_next;
using namespace principia::quantities::_si;

// Pointers used for indirect calls, set by `ConfigureElementaryFunctions`.
inline ElementaryFunctionPointer<double> cos = nullptr;
inline ElementaryFunctionPointer<double> sin = nullptr;
inline ElementaryFunctionPointer<SC<double>> sin_cos = nullptr;

inline nullptr_t static_initialization = [](){
  // By default, use the correctly-rounded implementation.  This can be
  // overridden based on the save.
  ConfigureElementaryFunctions(/*uses_correct_sin_cos=*/true);
  return nullptr;
}();

inline SC<double> __cdecl StdSinCos(double const x) {
  return {.sin = std::sin(x), .cos = std::cos(x)};
}

inline ElementaryFunctionsConfigurationSaver::
ElementaryFunctionsConfigurationSaver()
    : cos_(cos), sin_(sin), sin_cos_(sin_cos) {
  CHECK(!active_);
  active_ = true;
}

inline ElementaryFunctionsConfigurationSaver::
~ElementaryFunctionsConfigurationSaver() {
  cos = cos_;
  sin = sin_;
  sin_cos = sin_cos_;
  active_ = false;
}

inline std::atomic_bool ElementaryFunctionsConfigurationSaver::active_ = false;

inline void ConfigureElementaryFunctions(bool const uses_correct_sin_cos) {
  // Beware!  This function cannot use `CanUseHardwareFMA` as it's called from a
  // static declaration and we have no way to know if that variable has already
  // been initialized.
  if (!uses_correct_sin_cos) {
    cos = &std::cos;
    sin = &std::sin;
    sin_cos = &StdSinCos;
  } else if (EarlyCanUseHardwareFMA()) {
    cos = &numerics::_sin_cos::Cos<FMAPresence::Present>;
    sin = &numerics::_sin_cos::Sin<FMAPresence::Present>;
    sin_cos = &numerics::_sin_cos::SinCos<FMAPresence::Present>;
  } else {
    cos = &numerics::_sin_cos::Cos<FMAPresence::Absent>;
    sin = &numerics::_sin_cos::Sin<FMAPresence::Absent>;
    sin_cos = &numerics::_sin_cos::SinCos<FMAPresence::Absent>;
  }
}

template<typename Q1, typename Q2>
  requires((boost_cpp_int<Q1> && boost_cpp_int<Q2>) ||
           (boost_cpp_rational<Q1> && boost_cpp_rational<Q2>))
Product<Q1, Q2> FusedMultiplyAdd(Q1 const& x,
                                 Q2 const& y,
                                 Product<Q1, Q2> const& z) {
  return x * y + z;
}

template<boost_cpp_bin_float Q1, boost_cpp_bin_float Q2>
Product<Q1, Q2> FusedMultiplyAdd(Q1 const& x,
                                 Q2 const& y,
                                 Product<Q1, Q2> const& z) {
  return fma(x, y, z);
}

template<convertible_to_quantity Q1, convertible_to_quantity Q2>
Product<Q1, Q2> FusedMultiplyAdd(Q1 const& x,
                                 Q2 const& y,
                                 Product<Q1, Q2> const& z) {
  return si::Unit<Product<Q1, Q2>> *
          numerics::_fma::FusedMultiplyAdd(x / si::Unit<Q1>,
                                          y / si::Unit<Q2>,
                                          z / si::Unit<Product<Q1, Q2>>);
}


template<typename Q1, typename Q2>
  requires((boost_cpp_int<Q1> && boost_cpp_int<Q2>) ||
           (boost_cpp_rational<Q1> && boost_cpp_rational<Q2>))
Product<Q1, Q2> FusedMultiplySubtract(Q1 const& x,
                                      Q2 const& y,
                                      Product<Q1, Q2> const& z) {
  return x * y - z;
}

template<boost_cpp_bin_float Q1, boost_cpp_bin_float Q2>
Product<Q1, Q2> FusedMultiplySubtract(Q1 const& x,
                                      Q2 const& y,
                                      Product<Q1, Q2> const& z) {
  return fma(x, y, -z);
}

template<convertible_to_quantity Q1, convertible_to_quantity Q2>
Product<Q1, Q2> FusedMultiplySubtract(Q1 const& x,
                                      Q2 const& y,
                                      Product<Q1, Q2> const& z) {
  return si::Unit<Product<Q1, Q2>> *
          numerics::_fma::FusedMultiplySubtract(x / si::Unit<Q1>,
                                                y / si::Unit<Q2>,
                                                z / si::Unit<Product<Q1, Q2>>);
}


template<typename Q1, typename Q2>
  requires((boost_cpp_int<Q1> && boost_cpp_int<Q2>) ||
           (boost_cpp_rational<Q1> && boost_cpp_rational<Q2>))
Product<Q1, Q2> FusedNegatedMultiplyAdd(Q1 const& x,
                                        Q2 const& y,
                                        Product<Q1, Q2> const& z) {
  return -x * y + z;
}

template<boost_cpp_bin_float Q1, boost_cpp_bin_float Q2>
Product<Q1, Q2> FusedNegatedMultiplyAdd(Q1 const& x,
                                        Q2 const& y,
                                        Product<Q1, Q2> const& z) {
  return fma(-x, y, z);
}

template<convertible_to_quantity Q1, convertible_to_quantity Q2>
Product<Q1, Q2> FusedNegatedMultiplyAdd(Q1 const& x,
                                        Q2 const& y,
                                        Product<Q1, Q2> const& z) {
  return si::Unit<Product<Q1, Q2>> * numerics::_fma::FusedNegatedMultiplyAdd(
                                          x / si::Unit<Q1>,
                                          y / si::Unit<Q2>,
                                          z / si::Unit<Product<Q1, Q2>>);
}


template<typename Q1, typename Q2>
  requires((boost_cpp_int<Q1> && boost_cpp_int<Q2>) ||
           (boost_cpp_rational<Q1> && boost_cpp_rational<Q2>))
Product<Q1, Q2> FusedNegatedMultiplySubtract(Q1 const& x,
                                             Q2 const& y,
                                             Product<Q1, Q2> const& z) {
  return -x * y - z;
}

template<boost_cpp_bin_float Q1, boost_cpp_bin_float Q2>
Product<Q1, Q2> FusedNegatedMultiplySubtract(Q1 const& x,
                                             Q2 const& y,
                                             Product<Q1, Q2> const& z) {
  return fma(-x, y, -z);
}

template<convertible_to_quantity Q1, convertible_to_quantity Q2>
Product<Q1, Q2> FusedNegatedMultiplySubtract(Q1 const& x,
                                             Q2 const& y,
                                             Product<Q1, Q2> const& z) {
  return si::Unit<Product<Q1, Q2>> *
          numerics::_fma::FusedNegatedMultiplySubtract(
              x / si::Unit<Q1>,
              y / si::Unit<Q2>,
              z / si::Unit<Product<Q1, Q2>>);
}


template<boost_cpp_number Q>
Q Abs(Q const& x) {
  return abs(x);
}

template<convertible_to_quantity Q>
FORCE_INLINE(inline)
Q Abs(Q const& x) {
  return si::Unit<Q> * std::abs(x / si::Unit<Q>);
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
#if PRINCIPIA_USE_SSE3_INTRINSICS()
  auto const x_128d = _mm_set_sd(x / si::Unit<Q>);
  return si::Unit<SquareRoot<Q>> * _mm_cvtsd_f64(_mm_sqrt_sd(x_128d, x_128d));
#else
  return si::Unit<SquareRoot<Q>> * std::sqrt(x / si::Unit<Q>);
#endif
}

template<typename Q>
CubeRoot<Q> Cbrt(Q const& x) {
  return si::Unit<CubeRoot<Q>> * numerics::_cbrt::Cbrt(x / si::Unit<Q>);
}

template<int N, typename Q>
NthRoot<Q, N> Root(Q const& x) {
  return si::Unit<NthRoot<Q, N>> * Root(N, x / si::Unit<Q>);
}

inline double Root(int const n, double const x) {
  return std::pow(x, 1.0 / n);
}

template<typename Q>
constexpr Q NextUp(Q const& x) {
  return si::Unit<Q> * numerics::_next::NextUp(x / si::Unit<Q>);
}

template<typename Q>
constexpr Q NextDown(Q const& x) {
  return si::Unit<Q> * numerics::_next::NextDown(x / si::Unit<Q>);
}

template<int exponent>
constexpr double Pow(double x) {
  // Use the Russian peasant algorithm for small exponents.
  if constexpr (exponent > 0 && exponent < 32) {
    // The end of the recursion is handled by the specializations below.
    auto const y = Pow<exponent / 2>(x);
    auto const y² = y * y;
    if constexpr (exponent % 2 == 1) {
      return y² * x;
    } else {
      return y²;
    }
  } else if constexpr (exponent < 0 && exponent > -32) {
    return 1 / Pow<-exponent>(x);
  } else {
    return std::pow(x, exponent);
  }
}

// Static specializations for frequently-used exponents, so that this gets
// turned into multiplications at compile time.

template<>
inline constexpr double Pow<0>(double) {
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
  if constexpr (boost_cpp_rational<Q>) {
    // It seems that Boost does not define `pow` for `cpp_rational`.
    return cpp_rational(pow(numerator(x), exponent),
                        pow(denominator(x), exponent));
  } else if constexpr (boost_cpp_number<Q>) {
    return pow(x, exponent);
  } else {
    return si::Unit<Exponentiation<Q, exponent>> *
           Pow<exponent>(x / si::Unit<Q>);
  }
}

inline double Sin(Angle const& α) {
  return sin(α / Radian);
}

inline double Cos(Angle const& α) {
  return cos(α / Radian);
}

inline double Tan(Angle const& α) {
  return std::tan(α / Radian);
}

inline SC<double> SinCos(Angle const& α) {
  return sin_cos(α / Radian);
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

template<typename Q>
  requires boost_cpp_bin_float<Q> || std::floating_point<Q>
Q Round(Q const& x) {
  if constexpr (boost_cpp_bin_float<Q>) {
    return static_cast<Q>(round(x));
  } else {
    return std::round(x);
  }
}

inline cpp_int Round(cpp_rational const& x) {
  // TODO(phl): This is clunky.  Use `divide_qr` or something.
  return static_cast<cpp_int>(round(static_cast<cpp_bin_float_50>(x)));
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

}  // namespace internal
}  // namespace _elementary_functions
}  // namespace numerics
}  // namespace principia
