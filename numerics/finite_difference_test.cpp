
#include "numerics/finite_difference.hpp"

#include "base/file.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "mathematica/mathematica.hpp"
#include "numerics/cbrt.hpp"

namespace principia {
namespace numerics {
namespace internal_finite_difference {

using quantities::Sqrt;
using quantities::Pow;

class FiniteDifferenceTest : public ::testing::Test {
 protected:
  template<int n>
  std::vector<std::vector<double>> estimate() const {
    std::vector<std::vector<double>> estimates;
    auto const f = [](double x) {
      return 1e9 + Pow<5>(x);
    };
      estimates.emplace_back();
      for (double step = 0x1p64; 1 + step != 1; step /= 2) {
        FixedVector<double, n> values;
        for (int i = 0; i < n; ++i) {
          values[i] = f(1 + i * step);
        }
        estimates.back().push_back(FiniteDifference(values, step, 0));
      }
    return estimates;
  }
};

TEST_F(FiniteDifferenceTest, Meow) {
  base::OFStream file(SOLUTION_DIR / "finite_difference");
  file << mathematica::Assign("estimates2", estimate<2>());
  file << mathematica::Assign("estimates5", estimate<5>());
  file << mathematica::Assign("estimates8", estimate<8>());
  file << mathematica::Assign("estimates16", estimate<16>());
}

}  // namespace internal_finite_difference
}  // namespace numerics
}  // namespace principia
