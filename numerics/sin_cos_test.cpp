#include "numerics/sin_cos.hpp"

#include <algorithm>
#include <limits>
#include <random>
#include <vector>

#include "absl/strings/str_cat.h"
#include "absl/synchronization/mutex.h"
#include "base/bundle.hpp"
#include "base/multiprecision.hpp"
#include "functions/multiprecision.hpp"
#include "glog/logging.h"
#include "gtest/gtest.h"
#include "numerics/fma.hpp"
#include "numerics/m128d.hpp"
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

using namespace principia::base::_bundle;
using namespace principia::base::_multiprecision;
using namespace principia::functions::_multiprecision;
using namespace principia::numerics::_fma;
using namespace principia::numerics::_m128d;
using namespace principia::numerics::_next;
using namespace principia::numerics::_sin_cos;
using namespace principia::testing_utilities::_almost_equals;

constexpr std::int64_t table_spacing_bits = 9;
constexpr double table_spacing_reciprocal = 1 << table_spacing_bits;

namespace m128d {

M128D const mantissa_index_bits = M128D::MakeFromBits(0x0000'0000'0000'01ffull);
M128D const accurate_table_index_addend(static_cast<double>(
    1LL << (std::numeric_limits<double>::digits - table_spacing_bits - 1)));

M128D const mantissa_reduce_shifter(
    static_cast<double>(1LL << (std::numeric_limits<double>::digits - 1)));
M128D const two_over_π(2.0 / π);

}  // namespace m128d

class SinCosTest : public ::testing::Test {
 protected:
  struct FunctionStatistics {
    cpp_bin_float_50 max_ulps_error = 0;
    double worst_argument = 0;
    cpp_bin_float_50 boost_fn_worst_argument = 0;
    double principia_fn_worst_argument = 0;
    std::int64_t incorrectly_rounded = 0;
    std::int64_t fallbacks = 0;
  };

  struct SinCosStatistics {
    FunctionStatistics sin;
    FunctionStatistics cos;
  };

  static void SetUpTestCase() {
    SetSlowPathsCallbacks(&CountSinFallbacks, &CountCosFallbacks);
  }

  static void CountSinFallbacks(double const θ) {
    absl::MutexLock lock(&lock_);
    ++sin_fallbacks_;
  }

  static void CountCosFallbacks(double const θ) {
    absl::MutexLock lock(&lock_);
    ++cos_fallbacks_;
  }

  static double Sin(double const θ) {
    return CanUseHardwareFMA ? numerics::_sin_cos::Sin<FMAPresence::Present>(θ)
                             : numerics::_sin_cos::Sin<FMAPresence::Absent>(θ);
  }

  static double Cos(double const θ) {
    return CanUseHardwareFMA ? numerics::_sin_cos::Cos<FMAPresence::Present>(θ)
                             : numerics::_sin_cos::Cos<FMAPresence::Absent>(θ);
  }

  static auto SinCos(double const θ) {
    return CanUseHardwareFMA
               ? numerics::_sin_cos::SinCos<FMAPresence::Present>(θ)
               : numerics::_sin_cos::SinCos<FMAPresence::Absent>(θ);
  }

  template<std::int64_t iterations_quantum>
  static SinCosStatistics RandomArgumentTest(double const lower_bound,
                                             double const upper_bound,
                                             bool const sin_cos,
                                             std::int64_t const seed) {
    std::mt19937_64 random(seed);
    std::uniform_real_distribution<> uniformly_at(lower_bound, upper_bound);
    std::uniform_int_distribution<> uniform_sign(0, 1);

    SinCosStatistics s;

    for (std::int64_t i = 0; i < iterations_quantum; ++i) {
      double const principia_argument =
          uniformly_at(random) * ((uniform_sign(random) << 1) - 1);
      auto const boost_argument = cpp_rational(principia_argument);
      double principia_sin;
      double principia_cos;
      if (sin_cos) {
        const auto [sin, cos] = SinCos(principia_argument);
        principia_sin = sin;
        principia_cos = cos;
      } else {
        principia_sin = Sin(principia_argument);
        principia_cos = Cos(principia_argument);
      }
      {
        auto const boost_sin = functions::_multiprecision::Sin(boost_argument);
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
        auto const boost_cos = functions::_multiprecision::Cos(boost_argument);
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
                                  double const upper_bound,
                                  bool const sin_cos) {
#if _DEBUG
    static constexpr std::int64_t iterations = 30'000;
    static constexpr std::int64_t iterations_quantum = 100;
#else
    static constexpr std::int64_t iterations = 10'000'000;
    static constexpr std::int64_t iterations_quantum = 10'000;
#endif
    static_assert(iterations % iterations_quantum == 0);

    {
      absl::MutexLock l(&lock_);
      sin_fallbacks_ = 0;
      cos_fallbacks_ = 0;
    }

    Bundle bundle;

    std::vector<SinCosStatistics> statistics(iterations / iterations_quantum);
    for (std::int64_t i = 0; i < statistics.size(); ++i) {
      bundle.Add([i, lower_bound, upper_bound, sin_cos, &statistics]() {
        statistics[i] = RandomArgumentTest<iterations_quantum>(
            lower_bound, upper_bound, sin_cos, /*seed=*/i);
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

    {
      absl::ReaderMutexLock l(&lock_);
      final_statistics.sin.fallbacks = sin_fallbacks_;
      final_statistics.cos.fallbacks = cos_fallbacks_;
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
                 << s.incorrectly_rounded / static_cast<double>(iterations)
                 << "; fallback probability: "
                 << (s.fallbacks == 0
                         ? "0"
                         : absl::StrCat("1/", iterations / s.fallbacks));
    };

    log_statistics("Sin", final_statistics.sin);
    log_statistics("Cos", final_statistics.cos);
  }

  static absl::Mutex lock_;
  static std::int64_t sin_fallbacks_;
  static std::int64_t cos_fallbacks_;
};

ABSL_CONST_INIT absl::Mutex SinCosTest::lock_(absl::kConstInit);
std::int64_t SinCosTest::sin_fallbacks_ = 0;
std::int64_t SinCosTest::cos_fallbacks_ = 0;

TEST_F(SinCosTest, AccurateTableIndex) {
  static constexpr std::int64_t iterations = 100;

  static const __m128d mantissa_index_bits =
      _mm_castsi128_pd(_mm_cvtsi64_si128(0x0000'0000'0000'01ff));
  std::mt19937_64 random(42);
  std::uniform_real_distribution<> uniformly_at(0, π / 4);

  for (std::int64_t i = 0; i < iterations; ++i) {
    double const x = uniformly_at(random);
    M128D const x_m128d(x);

    std::int64_t const n =
        _mm_cvtsd_si64(_mm_set_sd(x * table_spacing_reciprocal));

    std::int64_t const m = _mm_cvtsi128_si64(_mm_castpd_si128(
        _mm_and_pd(mantissa_index_bits,
                   _mm_set_sd(x + (1LL << (std::numeric_limits<double>::digits -
                                           table_spacing_bits - 1))))));

    std::int64_t const p = (m128d::mantissa_index_bits &
          (x_m128d + m128d::accurate_table_index_addend))
      .Bits<std::int64_t>();

    EXPECT_EQ(n, m);
    EXPECT_EQ(n, p);
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
    M128D const θ_m128d(θ);

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

    if constexpr (CanEmitFMAInstructions) {
      M128D const abs_θ_m128d = Abs(θ_m128d);
      M128D const sign_m128d = Sign(θ_m128d);
      M128D p_m128d =
          FusedMultiplyAdd(
              abs_θ_m128d, m128d::two_over_π, m128d::mantissa_reduce_shifter) -
          m128d::mantissa_reduce_shifter;
      p_m128d = p_m128d ^ sign_m128d;

      EXPECT_EQ(m_double, static_cast<double>(p_m128d));
    }
  }
}

TEST_F(SinCosTest, RandomSmall) {
  ParallelRandomArgumentTest(0, π / 4, /*sin_cos=*/false);
}

TEST_F(SinCosTest, RandomTwoTerms) {
  ParallelRandomArgumentTest(π / 4, 1 << 8, /*sin_cos=*/false);
}

TEST_F(SinCosTest, RandomThreeTerms) {
  ParallelRandomArgumentTest(1 << 8, 1 << 18, /*sin_cos=*/false);
}

TEST_F(SinCosTest, RandomLarge) {
  ParallelRandomArgumentTest(
      1 << 18, std::numeric_limits<double>::max() / 1.0e30, /*sin_cos=*/false);
}

TEST_F(SinCosTest, RandomSinCos) {
  ParallelRandomArgumentTest(0, 1 << 18, /*sin_cos=*/true);
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
