#include "numerics/sin_cos.hpp"

#include <algorithm>
#include <limits>
#include <random>
#include <vector>

#include "base/bundle.hpp"
#include "boost/multiprecision/cpp_int.hpp"
#include "functions/multiprecision.hpp"
#include "glog/logging.h"
#include "gtest/gtest.h"
#include "numerics/next.hpp"
#include "quantities/numbers.hpp"  // 🧙 For π.
#include "testing_utilities/almost_equals.hpp"

// This test lives in `functions` to avoid pulling `boost` into `numerics`.
// It uses neither the `functions` nor the `numerics` namespace so that the Sin
// and Cos from both (`principia::functions::_multiprecision` and
// `principia::numerics::_sin_cos`) are made visible by the using directives
// below.

namespace principia {
namespace functions_test {

using namespace boost::multiprecision;
using namespace principia::base::_bundle;
using namespace principia::functions::_multiprecision;
using namespace principia::numerics::_next;
using namespace principia::numerics::_sin_cos;
using namespace principia::testing_utilities::_almost_equals;
namespace sin_cos = principia::numerics::_sin_cos;

class SinCosTest : public ::testing::Test {
 protected:
  struct FunctionStatistics {
    cpp_bin_float_50 max_ulps_error = 0;
    double worst_argument = 0;
    cpp_bin_float_50 boost_fn_worst_argument = 0;
    double principia_fn_worst_argument = 0;
    std::int64_t incorrectly_rounded = 0;
  };

  struct SinCosStatistics {
    FunctionStatistics sin;
    FunctionStatistics cos;
  };

  static void SetUpTestCase() {
    sin_cos::StaticInitialization();
  }

  template<std::int64_t iterations_quantum>
  static SinCosStatistics RandomArgumentTest(double const lower_bound,
                                             double const upper_bound,
                                             std::int64_t const seed) {
    std::mt19937_64 random(seed);
    std::uniform_real_distribution<> uniformly_at(lower_bound, upper_bound);
    std::uniform_int_distribution<> uniform_sign(0, 1);

    SinCosStatistics s;

    for (std::int64_t i = 0; i < iterations_quantum; ++i) {
      double const principia_argument =
          uniformly_at(random) * ((uniform_sign(random) << 1) - 1);
      auto const boost_argument = cpp_rational(principia_argument);
      {
        auto const boost_sin = Sin(boost_argument);
        double const principia_sin = Sin(principia_argument);
        auto const sin_error =
            abs(boost_sin - static_cast<cpp_bin_float_50>(principia_sin));
        auto const ulp = NextUp(principia_sin) - principia_sin;
        auto const sin_ulps_error = sin_error / ulp;
        if (sin_ulps_error > s.sin.max_ulps_error) {
          s.sin.max_ulps_error = sin_ulps_error;
          s.sin.worst_argument = principia_argument;
          s.sin.boost_fn_worst_argument = boost_sin;
          s.sin.principia_fn_worst_argument = principia_sin;
        }
        if (sin_ulps_error > 0.5) {
          ++s.sin.incorrectly_rounded;
          LOG(ERROR) << "Sin: " << sin_ulps_error << " ulps at "
                     << std::setprecision(25) << principia_argument;
        }
      }
      {
        auto const boost_cos = Cos(boost_argument);
        double const principia_cos = Cos(principia_argument);
        auto const cos_error =
            abs(boost_cos - static_cast<cpp_bin_float_50>(principia_cos));
        auto const ulp = NextUp(principia_cos) - principia_cos;
        auto const cos_ulps_error = cos_error / ulp;
        if (cos_ulps_error > s.cos.max_ulps_error) {
          s.cos.max_ulps_error = cos_ulps_error;
          s.cos.worst_argument = principia_argument;
          s.cos.boost_fn_worst_argument = boost_cos;
          s.cos.principia_fn_worst_argument = principia_cos;
        }
        if (cos_ulps_error > 0.5) {
          ++s.cos.incorrectly_rounded;
          LOG(ERROR) << "Cos: " << cos_ulps_error << " ulps at "
                     << std::setprecision(25) << principia_argument;
        }
      }
    }

    EXPECT_LE(s.sin.max_ulps_error, 0.5);
    EXPECT_LE(s.cos.max_ulps_error, 0.5);
    EXPECT_EQ(s.sin.incorrectly_rounded, 0);
    EXPECT_EQ(s.cos.incorrectly_rounded, 0);

    return s;
  }

  void ParallelRandomArgumentTest(double const lower_bound,
                                  double const upper_bound) {
#if _DEBUG
    static constexpr std::int64_t iterations = 30'000;
    static constexpr std::int64_t iterations_quantum = 100;
#else
    static constexpr std::int64_t iterations = 10'000'000;
    static constexpr std::int64_t iterations_quantum = 10'000;
#endif
    static_assert(iterations % iterations_quantum == 0);

    Bundle bundle;

    std::vector<SinCosStatistics> statistics(iterations / iterations_quantum);
    for (std::int64_t i = 0; i < statistics.size(); ++i) {
      bundle.Add([i, lower_bound, upper_bound, &statistics]() {
        statistics[i] = RandomArgumentTest<iterations_quantum>(
            lower_bound, upper_bound, /*seed=*/i);
        return absl::OkStatus();
      });
    }
    absl::Status const status = bundle.Join();

    SinCosStatistics final_statistics;
    for (auto const& s : statistics) {
      if (s.sin.max_ulps_error > final_statistics.sin.max_ulps_error) {
        final_statistics.sin.max_ulps_error = s.sin.max_ulps_error;
        final_statistics.sin.worst_argument = s.sin.worst_argument;
        final_statistics.sin.boost_fn_worst_argument =
            s.sin.boost_fn_worst_argument;
        final_statistics.sin.principia_fn_worst_argument =
            s.sin.principia_fn_worst_argument;
      }
      final_statistics.sin.incorrectly_rounded += s.sin.incorrectly_rounded;
      if (s.cos.max_ulps_error > final_statistics.cos.max_ulps_error) {
        final_statistics.cos.max_ulps_error = s.cos.max_ulps_error;
        final_statistics.cos.worst_argument = s.cos.worst_argument;
        final_statistics.cos.boost_fn_worst_argument =
            s.cos.boost_fn_worst_argument;
        final_statistics.cos.principia_fn_worst_argument =
            s.cos.principia_fn_worst_argument;
      }
      final_statistics.cos.incorrectly_rounded += s.cos.incorrectly_rounded;
    }

    auto log_statistics = [](std::string_view const fn_name,
                             FunctionStatistics const& s) {
      LOG(ERROR) << fn_name << " error: " << s.max_ulps_error
                 << std::setprecision(25)
                 << " ulps for argument: " << s.worst_argument << " ("
                 << std::hexfloat << s.worst_argument << std::defaultfloat
                 << ") value: " << s.principia_fn_worst_argument << " ("
                 << std::hexfloat << s.principia_fn_worst_argument
                 << std::defaultfloat << ") vs. expected "
                 << s.boost_fn_worst_argument << " (" << std::hexfloat
                 << static_cast<double>(s.boost_fn_worst_argument)
                 << std::defaultfloat << "); incorrectly rounded probability: "
                 << std::setprecision(3)
                 << s.incorrectly_rounded / static_cast<double>(iterations);
    };

    log_statistics("Sin", final_statistics.sin);
    log_statistics("Cos", final_statistics.cos);
  }
};

//TODO(phl)update this and the next
TEST_F(SinCosTest, AccurateTableIndex) {
  static constexpr std::int64_t iterations = 100;

  constexpr std::int64_t table_spacing_bits = 9;
  constexpr double table_spacing_reciprocal = 1 << table_spacing_bits;
  static const __m128d mantissa_index_bits =
      _mm_castsi128_pd(_mm_cvtsi64_si128(0x0000'0000'0000'01ff));
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> uniformly_at(0, π / 4);

  for (std::int64_t i = 0; i < iterations; ++i) {
    double const x = uniformly_at(random);

    std::int64_t const n =
        _mm_cvtsd_si64(_mm_set_sd(x * table_spacing_reciprocal));

    std::int64_t const m = _mm_cvtsi128_si64(_mm_castpd_si128(
        _mm_and_pd(mantissa_index_bits,
                   _mm_set_sd(x + (1LL << (std::numeric_limits<double>::digits -
                                           table_spacing_bits - 1))))));

    EXPECT_EQ(n, m);
  }
}

TEST_F(SinCosTest, ReduceIndex) {
  static constexpr std::int64_t iterations = 100;

  static const __m128d sign_bit =
      _mm_castsi128_pd(_mm_cvtsi64_si128(0x8000'0000'0000'0000));
  static constexpr double mantissa_reduce_shifter =
      1LL << (std::numeric_limits<double>::digits - 1);
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> uniformly_at(-1000.0, 1000.0);

  for (std::int64_t i = 0; i < iterations; ++i) {
    double const θ = uniformly_at(random);

    __m128d const n_128d = _mm_round_sd(
        _mm_setzero_pd(), _mm_set_sd(θ * (2 / π)), _MM_FROUND_RINT);
    std::int64_t const n = _mm_cvtsd_si64(n_128d);

    double const abs_θ = std::abs(θ);
    __m128d const sign = _mm_and_pd(sign_bit, _mm_set_sd(θ));
    double const m_shifted = abs_θ * (2 / π) + mantissa_reduce_shifter;
    double const m_double = _mm_cvtsd_f64(
        _mm_xor_pd(_mm_set_sd(m_shifted - mantissa_reduce_shifter), sign));
    std::int64_t const m = _mm_cvtsd_si64(_mm_set_sd(m_double));

    EXPECT_EQ(n, m);
    EXPECT_EQ(m, m_double);
  }
}

TEST_F(SinCosTest, RandomSmall) {
  ParallelRandomArgumentTest(0, π / 4);
}

TEST_F(SinCosTest, RandomTwoTerms) {
  ParallelRandomArgumentTest(π / 4, 1 << 8);
}

TEST_F(SinCosTest, RandomThreeTerms) {
  ParallelRandomArgumentTest(1 << 8, 1 << 18);
}

TEST_F(SinCosTest, RandomLarge) {
  ParallelRandomArgumentTest(1 << 18,
                             std::numeric_limits<double>::max() / 1.0e30);
}

// Values for which the base algorithm gives an error of 1 ULP.
TEST_F(SinCosTest, HardRounding) {
  EXPECT_THAT(Sin(1.777288458404935767021016),
              AlmostEquals(0.9787561457198967196367773, 0));
  EXPECT_THAT(Cos(3.912491942337291916942377),
              AlmostEquals(-0.7172843528140595004137653, 0));
  EXPECT_THAT(Sin(5.528810471911395296729097),
              AlmostEquals(-0.6848332450871304488693046, 0));
  EXPECT_THAT(Sin(2.670333644894535396474566),
              AlmostEquals(0.4540084183741445456039384, 0));
  EXPECT_THAT(Cos(1.486604973422413600303571),
              AlmostEquals(0.0840919279825555407437241, 0));
  EXPECT_THAT(Sin(-2.496680544289484160458414),
              AlmostEquals(-0.6011282027544306294509797, 0));
  EXPECT_THAT(Sin(3.348980952786005715893225),
              AlmostEquals(-0.2059048676737040700634683, 0));
  EXPECT_THAT(Sin(3.523452575387961971387085),
              AlmostEquals(-0.3726470704519433130297035, 0));
  EXPECT_THAT(Cos(-6.265702600230396157598989),
              AlmostEquals(0.99984718137127853720984932, 0));
  EXPECT_THAT(Sin(1.881458091523454001503524),
              AlmostEquals(0.9521314843257784876761001, 0));
  EXPECT_THAT(Sin(-1.763163156774038675678185),
              AlmostEquals(-0.9815544881044536151825223, 0));
  EXPECT_THAT(Cos(-3.885819786017697730073905),
              AlmostEquals(-0.7356116652133562472394118, 0));
  EXPECT_THAT(Sin(-2.58105062034143273308473),
              AlmostEquals(-0.5316453603071467637339815, 0));
  EXPECT_THAT(Sin(1.657419885978818285821035),
              AlmostEquals(0.99625052493662308306103561, 0));
  EXPECT_THAT(Sin(5.094301519947547873812255),
              AlmostEquals(-0.9279535374988051033005616, 0));
  EXPECT_THAT(Sin(5.262137362438826571064965),
              AlmostEquals(-0.8526560125576488347044409, 0));
  EXPECT_THAT(Cos(-5.026994177012682030181168),
              AlmostEquals(0.3094410694753661206223057, 0));
  EXPECT_THAT(Cos(0.2388111698570396512764091),
              AlmostEquals(0.9716198764286143041422587, 0));
}

TEST_F(SinCosTest, HardReduction) {
  // A previous implementation was incorrectly reducing this angle.
  EXPECT_THAT(Sin(0x5ad7bcep-8), AlmostEquals(0.99999999999999999818, 0));
  EXPECT_THAT(Cos(0x5ad7bcep-8),
              AlmostEquals(1.9060601347136373148e-9, 0));

  // Muller's bad angle.
  EXPECT_THAT(Sin(0x16ac5b262ca1ffp797), AlmostEquals(1.0, 0));
  EXPECT_THAT(Cos(0x16ac5b262ca1ffp797),
              AlmostEquals(-4.687165924254627611122582801963884e-19, 0));
}

}  // namespace functions_test
}  // namespace principia
