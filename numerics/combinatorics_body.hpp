
#pragma once

#include "numerics/combinatorics.hpp"

#include <cstdint>

namespace principia {
namespace numerics {

constexpr std::int64_t Binomial(std::int64_t n, std::int64_t k) {
  std::int64_t result = FallingFactorial(n, k);
  for (std::int64_t i = 1; i <= k; ++i) {
    result /= i;
  }
  return result;
}

constexpr std::int64_t FallingFactorial(std::int64_t n, std::int64_t k) {
  std::int64_t result = 1;
  for (std::int64_t i = 0; i < k; ++i) {
    result *= n - i;
  }
  return result;
}

}  // namespace numerics
}  // namespace principia
