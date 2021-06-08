#pragma once

#include "numerics/fma.hpp"

#include "glog/logging.h"

namespace principia {
namespace numerics {
namespace internal_fma {

inline double FusedMultiplyAdd(double const a, double const b, double const c) {
  if constexpr (CanEmitFMAInstructions) {
    return _mm_cvtsd_f64(
        _mm_fmadd_sd(_mm_set_sd(a), _mm_set_sd(b), _mm_set_sd(c)));
  } else {
    LOG(FATAL) << "Clang cannot use FMA without VEX-encoding everything";
  }
}

inline double FusedMultiplySubtract(double const a,
                                    double const b,
                                    double const c) {
  if constexpr (CanEmitFMAInstructions) {
    return _mm_cvtsd_f64(
        _mm_fmsub_sd(_mm_set_sd(a), _mm_set_sd(b), _mm_set_sd(c)));
  } else {
    LOG(FATAL) << "Clang cannot use FMA without VEX-encoding everything";
  }
}
inline double FusedNegatedMultiplyAdd(double const a,
                                      double const b,
                                      double const c) {
  if constexpr (CanEmitFMAInstructions) {
    return _mm_cvtsd_f64(
        _mm_fnmadd_sd(_mm_set_sd(a), _mm_set_sd(b), _mm_set_sd(c)));
  } else {
    LOG(FATAL) << "Clang cannot use FMA without VEX-encoding everything";
  }
}

inline double FusedNegatedMultiplySubtract(double const a,
                                           double const b,
                                           double const c) {
  if constexpr (CanEmitFMAInstructions) {
    return _mm_cvtsd_f64(
        _mm_fnmsub_sd(_mm_set_sd(a), _mm_set_sd(b), _mm_set_sd(c)));
  } else {
    LOG(FATAL) << "Clang cannot use FMA without VEX-encoding everything";
  }
}

}  // namespace internal_fma
}  // namespace numerics
}  // namespace principia
