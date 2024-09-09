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
inline bool const UseHardwareFMA =
    (CanEmitFMAInstructions && CPUIDFeatureFlag::FMA.IsSet());
#else
inline bool const UseHardwareFMA = false;
#endif

// The policy used for emitting FMA instructions.  This type is not used by this
// file, but is declared here for the convenience of the clients.  The intended
// semantics are:
// * `Auto`: FMA is used if supported by the processor, the decision must be
//   made dynamically by calling `UseHardwareFMA`.
// * `Disallow`: FMA is never used.
// * `Force`: FMA is always used.  The caller is expected to determine upstream
//   if FMA is supported by the processor by calling `UseHardwareFMA`.
enum class FMAPolicy {
  Auto = 0,
  Disallow = 1,
  Force = 2,
};

// The functions in this file unconditionally wrap the appropriate intrinsics.
// The caller may only use them if `UseHardwareFMA` is true.

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
using internal::FMAPolicy;
using internal::FusedMultiplyAdd;
using internal::FusedMultiplySubtract;
using internal::FusedNegatedMultiplyAdd;
using internal::FusedNegatedMultiplySubtract;
using internal::UseHardwareFMA;

}  // namespace _fma
}  // namespace numerics
}  // namespace principia

#include "numerics/fma_body.hpp"
