#pragma once

#include <immintrin.h>

#include "base/cpuid.hpp"
#include "base/macros.hpp"

namespace principia {
namespace numerics {
namespace internal_fma {

using base::CPUFeatureFlags;

// The functions in this file unconditionally wrap the appropriate intrinsics.
// The caller may only use them if |UseHardwareFMA| is true.
#if PRINCIPIA_USE_FMA_IF_AVAILABLE
inline bool const UseHardwareFMA = HasCPUFeatures(CPUFeatureFlags::FMA);
#else
inline bool const UseHardwareFMA = false;
#endif

// ⟦ab + c⟧.
inline double FusedMultiplyAdd(double a, double b, double c);

// ⟦ab - c⟧.
inline double FusedMultiplySubtract(double a, double b, double c);

// ⟦-ab + c⟧.
inline double FusedNegatedMultiplyAdd(double a, double b, double c);

// ⟦-ab - c⟧.
inline double FusedNegatedMultiplySubtract(double a, double b, double c);

}  // namespace internal_fma

using internal_fma::FusedMultiplyAdd;
using internal_fma::FusedMultiplySubtract;
using internal_fma::FusedNegatedMultiplyAdd;
using internal_fma::FusedNegatedMultiplySubtract;
using internal_fma::UseHardwareFMA;

}  // namespace numerics
}  // namespace principia

#include "numerics/fma_body.hpp"
