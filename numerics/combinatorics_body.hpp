
#pragma once

#include "numerics/combinatorics.hpp"

#include <cstdint>

namespace principia {
namespace numerics {

constexpr std::int64_t Binomial(std::int64_t const n, std::int64_t const k) {
  return FallingFactorial(n, k) / Factorial(k);
}

constexpr std::int64_t DoubleFactorial(std::int64_t n) {
  std::int64_t result = 1;
  for (std::int64_t i = n; i >= 1; i -= 2) {
    result *= i;
  }
  return result;
}

constexpr std::int64_t Factorial(std::int64_t n) {
  return FallingFactorial(n, n);
}

constexpr std::int64_t FallingFactorial(std::int64_t const n,
                                        std::int64_t const k) {
  std::int64_t result = 1;
  for (std::int64_t i = 0; i < k; ++i) {
    result *= n - i;
  }
  return result;
}

}  // namespace numerics
}  // namespace principia
