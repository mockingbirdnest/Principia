#pragma once

#include <cmath>
#include <immintrin.h>

#include "base/cpuid.hpp"
#include "base/macros.hpp"

namespace principia {
namespace numerics {
namespace internal_fma {

using base::FeatureFlags;

#if PRINCIPIA_USE_FMA_IF_AVAILABLE
inline bool UseHardwareFMA = HasCPUFeatures(FeatureFlags::FMA);
#else
inline bool UseHardwareFMA = false;
#endif

// ⟦ab + c⟧.
double FusedMultiplyAdd(double a, double b, double c);

// ⟦ab - c⟧.
double FusedMultiplySubtract(double a, double b, double c);

// ⟦-ab + c⟧.
double FusedNegatedMultiplyAdd(double a, double b, double c);

// ⟦-ab - c⟧.
double FusedNegatedMultiplySubtract(double a, double b, double c);

}  // namespace internal_fma

using internal_fma::UseHardwareFMA;
using internal_fma::FusedMultiplyAdd;
using internal_fma::FusedMultiplySubtract;
using internal_fma::FusedNegatedMultiplyAdd;
using internal_fma::FusedNegatedMultiplySubtract;

}  // namespace numerics
}  // namespace principia

#include "numerics/fma_body.hpp"
