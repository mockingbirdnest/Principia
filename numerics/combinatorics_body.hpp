
#pragma once

#include "numerics/combinatorics.hpp"

namespace principia {
namespace numerics {

constexpr int Binomial(int n, int k) {
  return k == 0 ? 1 : n == k ? 1 : Binomial(n - 1, k - 1) + Binomial(n - 1, k);
}

constexpr int FallingFactorial(int n, int k) {
  return k == 0 ? 1 : n * FallingFactorial(n - 1, k - 1);
}

}  // namespace numerics
}  // namespace principia
