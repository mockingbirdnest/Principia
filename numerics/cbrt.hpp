#pragma once

namespace principia {
namespace numerics {
namespace _cbrt {
namespace internal {

// Computes ∛y, correctly rounded to nearest.
double Cbrt(double y);

// Specific methods, uncorrected rounding, and the corrector function exposed
// for testing.
// See documentation/cbrt.pdf, appendix F, for the naming of the methods and the
// misrounding rates of the unfaithful versions.

enum class Rounding {
  Faithful,
  Correct,
};

namespace method_3²ᴄZ5¹ {
template<Rounding rounding>
double Cbrt(double y);
extern template double Cbrt<Rounding::Faithful>(double y);
extern template double Cbrt<Rounding::Correct>(double y);
}  // namespace method_3²ᴄZ5¹

namespace method_5²Z4¹FMA {
template<Rounding rounding>
double Cbrt(double y);
extern template double Cbrt<Rounding::Faithful>(double y);
extern template double Cbrt<Rounding::Correct>(double y);
}  // namespace method_5²Z4¹FMA

// Computes one additional bit of ∛y, rounded toward 0.
// All arguments must be positive.
// a is the already-computed approximation of the cube root.
// b is the bit being computed; it must be a power of 2.
// The least significant bit of a must be greater than b.
// ∛y must lie in [a, a + 2b[.
// The result is the value of ∛y ≥ a + b, i.e., the value of the bit b in the
// binary expansion of ∛y.
bool CbrtOneBit(double y, double a, double b);

}  // namespace internal

using internal::Cbrt;

}  // namespace _cbrt
}  // namespace numerics
}  // namespace principia

namespace principia::numerics {
using namespace principia::numerics::_cbrt;
}  // namespace principia::numerics
