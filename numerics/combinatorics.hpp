
#pragma once

namespace principia {
namespace numerics {

// n choose k.
constexpr int Binomial(int n, int k);

constexpr int FallingFactorial(int n, int k);

}  // namespace numerics
}  // namespace principia

#include "numerics/combinatorics_body.hpp"
