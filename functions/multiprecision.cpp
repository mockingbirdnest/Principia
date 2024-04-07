#include "functions/multiprecision.hpp"

#include <cstdint>

#include "numerics/combinatorics.hpp"

namespace principia {
namespace functions {
namespace _multiprecision {
namespace internal {

using namespace principia::numerics::_combinatorics;

cpp_bin_float_50 Sin(cpp_rational const& α) {
  cpp_rational sum;
  std::int64_t sign = 1;
  cpp_rational α²ⁿ⁺¹ = α;
  for (std::int64_t n = 0; n < 10; ++n) {
    sum += cpp_rational(sign, Factorial(2 * n + 1)) * α²ⁿ⁺¹;
    sign *= -1;
    α²ⁿ⁺¹ *= α * α;
  }
  return static_cast<cpp_bin_float_50>(sum);
}

}  // namespace internal
}  // namespace _multiprecision
}  // namespace functions
}  // namespace principia
