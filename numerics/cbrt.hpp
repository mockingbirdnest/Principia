#pragma once

namespace principia {
namespace numerics {
// Computes ∛y, correctly rounded to nearest.
double Cbrt(double y);

// Specific methods and uncorrected rounding exposed for testing.

enum class Rounding {
  Faifthful,
  Correct,
};

namespace method_3²ᴄZ5¹ {
template<Rounding rounding = Rounding::Correct>
double Cbrt(double y);
extern template double Cbrt<Rounding::Faifthful>(double y);
extern template double Cbrt<Rounding::Correct>(double y);
}  // namespace method_3²ᴄZ5¹

namespace method_5²Z4¹FMA {
template<Rounding rounding = Rounding::Correct>
double Cbrt(double y);
extern template double Cbrt<Rounding::Faifthful>(double y);
extern template double Cbrt<Rounding::Correct>(double y);
}  // namespace method_5²Z4¹FMA

}  // namespace numerics
}  // namespace principia
