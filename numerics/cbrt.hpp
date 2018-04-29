#pragma once

namespace principia {
namespace numerics {

// Computes ∛y with a maximal error in [0.50005, 0.50022] ULPs; the result is
// incorrectly rounded for approximately 5 inputs per million.
double cbrt(double y);

}  // namespace numerics
}  // namespace principia
