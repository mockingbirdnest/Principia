#include "numerics/sin_cos.hpp"

#include <algorithm>
#include <limits>
#include <random>

#include "boost/multiprecision/cpp_int.hpp"
#include "functions/multiprecision.hpp"
#include "glog/logging.h"
#include "gtest/gtest.h"
#include "numerics/next.hpp"
#include "quantities/numbers.hpp"

// This test lives in `functions` to avoid pulling `boost` into `numerics`.
namespace principia {
namespace numerics {
namespace _sin_cos {

using namespace boost::multiprecision;
using namespace principia::numerics::_next;

class SinCosTest : public ::testing::Test {};

TEST_F(SinCosTest, Random) {
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> uniformly_at(-π / 4, π / 4);

  cpp_bin_float_50 max_sin_ulps_error = 0;
  cpp_bin_float_50 max_cos_ulps_error = 0;
  double worst_sin_argument = 0;
  double worst_cos_argument = 0;

#if _DEBUG
  static constexpr std::int64_t iterations = 100;
#else
  static constexpr std::int64_t iterations = 300'000;
#endif

  for (std::int64_t i = 0; i < iterations; ++i) {
    double const principia_argument = uniformly_at(random);
    auto const boost_argument = cpp_rational(principia_argument);
    {
      auto const boost_sin =
          functions::_multiprecision::Sin(boost_argument);
      double const principia_sin = Sin(principia_argument);
      auto const sin_error =
          abs(boost_sin - static_cast<cpp_bin_float_50>(principia_sin));
      auto const ulp = NextUp(principia_sin) - principia_sin;
      auto const sin_ulps_error = sin_error / ulp;
      if (sin_ulps_error > max_sin_ulps_error) {
        max_sin_ulps_error = sin_ulps_error;
        worst_sin_argument = principia_argument;
      }
    }
    {
      auto const boost_cos =
          functions::_multiprecision::Cos(boost_argument);
      double const principia_cos = Cos(principia_argument);
      auto const cos_error =
          abs(boost_cos - static_cast<cpp_bin_float_50>(principia_cos));
      auto const ulp = NextUp(principia_cos) - principia_cos;
      auto const cos_ulps_error = cos_error / ulp;
      if (cos_ulps_error > max_cos_ulps_error) {
        max_cos_ulps_error = cos_ulps_error;
        worst_cos_argument = principia_argument;
      }
    }
  }

  // This implementation is not quite correctly rounded, but not far from it.
  EXPECT_LE(max_sin_ulps_error, 0.500003);
  EXPECT_LE(max_cos_ulps_error, 0.500001);

  LOG(ERROR) << "Sin error: " << max_sin_ulps_error << std::setprecision(25)
             << " ulps for argument: " << worst_sin_argument
             << " value: " << Sin(worst_sin_argument);
  LOG(ERROR) << "Cos error: " << max_cos_ulps_error << std::setprecision(25)
             << " ulps for argument: " << worst_cos_argument
             << " value: " << Cos(worst_cos_argument);
}

}  // namespace _sin_cos
}  // namespace numerics
}  // namespace principia
