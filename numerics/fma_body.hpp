#pragma once

#include "numerics/fma.hpp"

#include "glog/logging.h"

namespace principia {
namespace numerics {
namespace _fma {
namespace internal {

inline double FusedMultiplyAdd(double const a, double const b, double const c) {
  return _mm_cvtsd_f64(
      _mm_fmadd_sd(_mm_set_sd(a), _mm_set_sd(b), _mm_set_sd(c)));
}

inline double FusedMultiplySubtract(double const a,
                                    double const b,
                                    double const c) {
  return _mm_cvtsd_f64(
      _mm_fmsub_sd(_mm_set_sd(a), _mm_set_sd(b), _mm_set_sd(c)));
}
inline double FusedNegatedMultiplyAdd(double const a,
                                      double const b,
                                      double const c) {
    return _mm_cvtsd_f64(
        _mm_fnmadd_sd(_mm_set_sd(a), _mm_set_sd(b), _mm_set_sd(c)));
}

inline double FusedNegatedMultiplySubtract(double const a,
                                           double const b,
                                           double const c) {
  return _mm_cvtsd_f64(
      _mm_fnmsub_sd(_mm_set_sd(a), _mm_set_sd(b), _mm_set_sd(c)));
}

}  // namespace internal
}  // namespace _fma
}  // namespace numerics
}  // namespace principia
