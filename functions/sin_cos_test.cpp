#include "numerics/sin_cos.hpp"

#include <algorithm>
#include <limits>
#include <random>

#include "boost/multiprecision/cpp_int.hpp"
#include "functions/multiprecision.hpp"
#include "glog/logging.h"
#include "gtest/gtest.h"
#include "quantities/numbers.hpp"

// This test lives in `functions` to avoid pulling `boost` into `numerics`.
namespace principia {
namespace numerics {
namespace _sin_cos {

using namespace boost::multiprecision;

class SinCosTest : public ::testing::Test {};

TEST_F(SinCosTest, Random) {
  std::mt19937_64 random(42);
  // TODO(phl): Negative angles.
  std::uniform_real_distribution<> uniformly_at(0, π / 4);

  cpp_bin_float_50 max_sin_ulps_error = 0;
  cpp_bin_float_50 max_cos_ulps_error = 0;
  double worst_sin_argument = 0;
  double worst_cos_argument = 0;

#if _DEBUG
  static constexpr std::int64_t iterations = 100;
#else
  static constexpr std::int64_t iterations = 100'000;
#endif

  for (std::int64_t i = 0; i < iterations; ++i) {
    double const argument = uniformly_at(random);
    {
      auto const boost_sin =
          functions::_multiprecision::Sin(cpp_rational(argument));
      double const principia_sin = Sin(argument);
      int sin_exponent;
      std::frexp(principia_sin, &sin_exponent);
      auto const sin_error =
          abs(boost_sin - static_cast<cpp_bin_float_50>(principia_sin));
      auto const sin_ulps_error =
          ldexp(sin_error, sin_exponent + std::numeric_limits<double>::digits);
      if (sin_ulps_error > max_sin_ulps_error) {
        max_sin_ulps_error = sin_ulps_error;
        worst_sin_argument = argument;
      }
    }
    {
      auto const boost_cos =
          functions::_multiprecision::Cos(cpp_rational(argument));
      double const principia_cos = Cos(argument);
      int cos_exponent;
      std::frexp(principia_cos, &cos_exponent);
      auto const cos_error =
          abs(boost_cos - static_cast<cpp_bin_float_50>(principia_cos));
      auto const cos_ulps_error =
          ldexp(cos_error, cos_exponent + std::numeric_limits<double>::digits);
      if (cos_ulps_error > max_cos_ulps_error) {
        max_cos_ulps_error = cos_ulps_error;
        worst_cos_argument = argument;
      }
    }
  }

  LOG(ERROR) << "SIN " << max_sin_ulps_error << std::setprecision(20) << " "
             << worst_sin_argument << " " << Sin(worst_sin_argument);
  LOG(ERROR) << "COS " << max_cos_ulps_error << std::setprecision(20) << " "
             << worst_cos_argument << " " << Cos(worst_cos_argument);
}

}  // namespace _sin_cos
}  // namespace numerics
}  // namespace principia
