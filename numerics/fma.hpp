#pragma once

#include <immintrin.h>

#include "base/cpuid.hpp"
#include "base/macros.hpp"  // 🧙 For PRINCIPIA_USE_FMA_IF_AVAILABLE.

namespace principia {
namespace numerics {
namespace _fma {
namespace internal {

using namespace principia::base::_cpuid;

// With clang, using FMA requires VEX-encoding everything; see #3019.
#if PRINCIPIA_COMPILER_MSVC
constexpr bool CanEmitFMAInstructions = true;
#else
constexpr bool CanEmitFMAInstructions = false;
#endif

#if PRINCIPIA_USE_FMA_IF_AVAILABLE()
#define PRINCIPIA_USE_HARDWARE_FMA_DEFAULT \
  (CanEmitFMAInstructions && HasCPUFeatures(CPUFeatureFlags::FMA));
#else
#define PRINCIPIA_USE_HARDWARE_FMA_DEFAULT false;
#endif

#if PRINCIPIA_CAN_OVERRIDE_FMA_USAGE_AT_RUNTIME
class FMAPreventer {
 public:
  FMAPreventer() {
    use_hardware_fma = false;
  }
  ~FMAPreventer() {
    use_hardware_fma = PRINCIPIA_USE_HARDWARE_FMA_DEFAULT;
  }
  static bool use_hardware_fma;
};
bool FMAPreventer::use_hardware_fma = PRINCIPIA_USE_HARDWARE_FMA_DEFAULT;
inline bool const& UseHardwareFMA = FMAPreventer::use_hardware_fma;
#else
// The functions in this file unconditionally wrap the appropriate intrinsics.
// The caller may only use them if |UseHardwareFMA| is true.
inline bool const UseHardwareFMA = PRINCIPIA_USE_HARDWARE_FMA_DEFAULT;
class FMAPreventer;  // Undefined.
#endif

#undef PRINCIPIA_USE_HARDWARE_FMA_DEFAULT


// ⟦ab + c⟧.
inline double FusedMultiplyAdd(double a, double b, double c);

// ⟦ab - c⟧.
inline double FusedMultiplySubtract(double a, double b, double c);

// ⟦-ab + c⟧.
inline double FusedNegatedMultiplyAdd(double a, double b, double c);

// ⟦-ab - c⟧.
inline double FusedNegatedMultiplySubtract(double a, double b, double c);

}  // namespace internal

using internal::CanEmitFMAInstructions;
using internal::FMAPreventer;
using internal::FusedMultiplyAdd;
using internal::FusedMultiplySubtract;
using internal::FusedNegatedMultiplyAdd;
using internal::FusedNegatedMultiplySubtract;
using internal::UseHardwareFMA;

}  // namespace _fma
}  // namespace numerics
}  // namespace principia

#include "numerics/fma_body.hpp"
