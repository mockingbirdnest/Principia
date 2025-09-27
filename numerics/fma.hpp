#pragma once

#include <immintrin.h>

#include "base/cpuid.hpp"
#include "base/macros.hpp"  // 🧙 For PRINCIPIA_USE_FMA_IF_AVAILABLE.

namespace principia {
namespace numerics {
namespace _fma {
namespace internal {

using namespace principia::base::_cpuid;

// Same as `CanUseHardwareFMA`, but may be used in static contexts where we
// don't know if `CanUseHardwareFMA` has already been initialized.
#if PRINCIPIA_HAS_FMA
constexpr bool EarlyCanUseHardwareFMA() { return true; }
#else
bool EarlyCanUseHardwareFMA();
#endif

#if PRINCIPIA_HAS_FMA || PRINCIPIA_COMPILER_MSVC
constexpr bool CanEmitFMAInstructions = true;
#else
constexpr bool CanEmitFMAInstructions = false;
#endif

#if PRINCIPIA_USE_FMA_IF_AVAILABLE()
#if PRINCIPIA_HAS_FMA
constexpr bool CanUseHardwareFMA = true;
#else
inline bool const CanUseHardwareFMA = EarlyCanUseHardwareFMA();
#endif
#else
inline bool const CanUseHardwareFMA = false;
#endif

// Whether FMA support is present.  This type is not used by this file, but is
// declared here for the convenience of clients.  The intended semantics are:
// * `Unknown`: FMA may or may not be present, the bit `CanUseHardwareFMA` must
//   be tested because trying to execute an FMA instruction.
// * `Absent`: Hardware FMA is not availabile, an attempt to use an FMA
//   instruction will cause a compile-time or run-time error.
// * `Present`: Hardware FMA support is available.
enum class FMAPresence {
  Unknown = 0,
  Absent = 1,
  Present = 2,
};

// The policy used for emitting FMA instructions.  This type is not used by this
// file, but is declared here for the convenience of the clients.  The intended
// semantics are:
// * `Auto`: FMA is used if supported by the processor, the decision must be
//   made dynamically by calling `CanUseHardwareFMA`.
// * `Disallow`: FMA is never used.
enum class FMAPolicy {
  Auto = 0,
  Disallow = 1,
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
using internal::CanUseHardwareFMA;
using internal::EarlyCanUseHardwareFMA;
using internal::FMAPolicy;
using internal::FMAPresence;
using internal::FusedMultiplyAdd;
using internal::FusedMultiplySubtract;
using internal::FusedNegatedMultiplyAdd;
using internal::FusedNegatedMultiplySubtract;

}  // namespace _fma
}  // namespace numerics
}  // namespace principia

#include "numerics/fma_body.hpp"
