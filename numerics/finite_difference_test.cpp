
#include "numerics/finite_difference.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "mathematica/mathematica.hpp"
#include "numerics/cbrt.hpp"

namespace principia {
namespace numerics {
namespace internal_finite_difference {

using mathematica::ToMathematica;
using quantities::Sqrt;
using quantities::Pow;

class FiniteDifferenceTest : public ::testing::Test {};

TEST_F(FiniteDifferenceTest, Meow) {
  std::vector<std::vector<double>> estimates;
  auto const f = [](double x) { return 1e12 - Pow<2>(x); };
  for (IMPL = 1; IMPL <= 4; ++IMPL) {
    estimates.emplace_back();
    for (double step = 0x1p64; 1 + step != 1; step /= 2) {
      FixedVector<double, 6> values = {f(1),
                                       f(1 + step),
                                       f(1 + 2 * step),
                                       f(1 + 3 * step),
                                       f(1 + 4 * step),
                                       f(1 + 5 * step)};
      estimates.back().push_back(FiniteDifference(values, step, 0));
    }
  }
  LOG(ERROR)<<ToMathematica(estimates);
}

}  // namespace internal_finite_difference
}  // namespace numerics
}  // namespace principia
