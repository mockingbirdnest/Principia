#include "functions/accurate_table_generator.hpp"

#include <cmath>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

#include "absl/strings/strip.h"
#include "boost/multiprecision/cpp_int.hpp"
#include "functions/multiprecision.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "mathematica/logger.hpp"
#include "mathematica/mathematica.hpp"
#include "quantities/numbers.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics_matchers.hpp"

namespace principia {
namespace functions {
namespace _accurate_table_generator {

using ::testing::AnyOf;
using ::testing::Eq;
using ::testing::Lt;
using ::testing::SizeIs;
using namespace boost::multiprecision;
using namespace principia::functions::_multiprecision;
using namespace principia::mathematica::_logger;
using namespace principia::mathematica::_mathematica;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_matchers;
using namespace principia::testing_utilities::_numerics_matchers;

class AccurateTableGeneratorTest : public ::testing::Test {
 protected:
  AccurateTableGeneratorTest() {
    FLAGS_v = 0;
    google::LogToStderr();
  }
};

#if !_DEBUG

TEST_F(AccurateTableGeneratorTest, GalSin5) {
  auto const x = GalExhaustiveSearch<5>({Sin}, 5.0 / 128.0);
  EXPECT_EQ(x, cpp_rational(2814749767106647, 72057594037927936));
  EXPECT_THAT(static_cast<double>(x),
              RelativeErrorFrom(5.0 / 128.0, IsNear(3.1e-14_(1))));

  // The stupid language doesn't allow printing a float in binary.  So to verify
  // that the |Sin| has zeroes in the right place we fiddle with the Mathematica
  // output.  We check the 7 bits starting at the last bit of the mantissa
  // (i.e., the bits starting at index 52) and verify that they are made of 5
  // zeroes surrounded by ones.  This validates that we actually shoot zeroes at
  // the right place.
  std::string const mathematica = ToMathematica(Sin(x),
                                                /*express_in=*/std::nullopt,
                                                /*base=*/2);
  std::string_view mantissa = mathematica;
  CHECK(absl::ConsumePrefix(&mantissa, "Times[2^^"));
  EXPECT_EQ("1000001", mantissa.substr(52, 7));
}

TEST_F(AccurateTableGeneratorTest, GalSinCos5) {
  auto const x = GalExhaustiveSearch<5>({Sin, Cos}, 95.0 / 128.0);
  EXPECT_EQ(x, cpp_rational(6685030696878177, 9007199254740992));
  EXPECT_THAT(static_cast<double>(x),
              RelativeErrorFrom(95.0 / 128.0, IsNear(1.5e-14_(1))));
  {
    std::string const mathematica = ToMathematica(Sin(x),
                                                  /*express_in=*/std::nullopt,
                                                  /*base=*/2);
    std::string_view mantissa = mathematica;
    CHECK(absl::ConsumePrefix(&mantissa, "Times[2^^"));
    EXPECT_EQ("1000001", mantissa.substr(52, 7));
  }
  {
    std::string const mathematica = ToMathematica(Cos(x),
                                                  /*express_in=*/std::nullopt,
                                                  /*base=*/2);
    std::string_view mantissa = mathematica;
    CHECK(absl::ConsumePrefix(&mantissa, "Times[2^^"));
    EXPECT_EQ("1000001", mantissa.substr(52, 7));
  }
}

TEST_F(AccurateTableGeneratorTest, GalMultisearchSinCos5) {
  static constexpr std::int64_t index_begin = 17;
  static constexpr std::int64_t index_end = 100;
  std::vector<cpp_rational> starting_arguments;
  for (std::int64_t i = index_begin; i < index_end; ++i) {
    starting_arguments.push_back(i / 128.0);
  }
  auto const xs = GalExhaustiveMultisearch<5>({Sin, Cos}, starting_arguments);
  EXPECT_THAT(xs, SizeIs(index_end - index_begin));
  for (std::int64_t i = 0; i < xs.size(); ++i) {
    auto const& x = xs[i];
    EXPECT_THAT(static_cast<double>(x),
                RelativeErrorFrom((i + index_begin) / 128.0, Lt(8.6e-13)));
    {
      std::string const mathematica = ToMathematica(Sin(x),
                                                    /*express_in=*/std::nullopt,
                                                    /*base=*/2);
      std::string_view mantissa = mathematica;
      CHECK(absl::ConsumePrefix(&mantissa, "Times[2^^"));
      EXPECT_EQ("00000", mantissa.substr(53, 5));
    }
    {
      std::string const mathematica = ToMathematica(Cos(x),
                                                    /*express_in=*/std::nullopt,
                                                    /*base=*/2);
      std::string_view mantissa = mathematica;
      CHECK(absl::ConsumePrefix(&mantissa, "Times[2^^"));
      EXPECT_EQ("00000", mantissa.substr(53, 5));
    }
  }
}

TEST_F(AccurateTableGeneratorTest, StehléZimmermannSinCos15) {
  double const x₀ = 1140850681.0 / 8589934592.0;
  double const u₀ = 4 * x₀;
  auto const sin = [](cpp_rational const& u) { return 4 * Sin(u / 4); };
  auto const cos = [](cpp_rational const& u) { return Cos(u / 4); };
  CHECK(0.5 <= u₀ & u₀ < 1.0) << u₀;
  CHECK(0.5 <= sin(u₀) && sin(u₀) < 1.0);
  CHECK(0.5 <= cos(u₀) && cos(u₀) < 1.0);
  AccuratePolynomial<cpp_rational, 2> sin_taylor2(
      {4 * cpp_rational(Sin(u₀ / 4)),
       cpp_rational(Cos(u₀ / 4)),
       -cpp_rational(Sin(u₀ / 4)) / 8},
      u₀);
  AccuratePolynomial<cpp_rational, 2> cos_taylor2(
      {cpp_rational(Cos(u₀ / 4)),
       -cpp_rational(Sin(u₀ / 4) / 4),
       -cpp_rational(Cos(u₀ / 4) / 32)},
      u₀);

  auto const u = StehléZimmermannSimultaneousSearch<15>(
      {sin, cos},
      {sin_taylor2, cos_taylor2},
      u₀,
      /*N=*/1ll << 53,
      /*T=*/1ll << 21);
  EXPECT_THAT(u,
              IsOkAndHolds(cpp_rational(4785074575333183, 9007199254740992)));
  EXPECT_THAT(static_cast<double>(*u),
              RelativeErrorFrom(u₀, Lt(1.3e-10)));
  {
    std::string const mathematica = ToMathematica(sin(*u),
                                                  /*express_in=*/std::nullopt,
                                                  /*base=*/2);
    std::string_view mantissa = mathematica;
    CHECK(absl::ConsumePrefix(&mantissa, "Times[2^^"));
    EXPECT_EQ("00000""00000""00000", mantissa.substr(53, 15));
  }
  {
    std::string const mathematica = ToMathematica(cos(*u),
                                                  /*express_in=*/std::nullopt,
                                                  /*base=*/2);
    std::string_view mantissa = mathematica;
    CHECK(absl::ConsumePrefix(&mantissa, "Times[2^^"));
    EXPECT_EQ("00000""00000""00000", mantissa.substr(53, 15));
  }
}

TEST_F(AccurateTableGeneratorTest, StehléZimmermannFullSinCos5NoScaling) {
  double const x₀ = 17.0 / 128;
  double const u₀ = 4 * x₀;
  auto const sin = [](cpp_rational const& u) { return 4 * Sin(u / 4); };
  auto const cos = [](cpp_rational const& u) { return Cos(u / 4); };
  CHECK(0.5 <= u₀ & u₀ < 1.0) << u₀;
  CHECK(0.5 <= sin(u₀) && sin(u₀) < 1.0);
  CHECK(0.5 <= cos(u₀) && cos(u₀) < 1.0);
  AccuratePolynomial<cpp_rational, 2> sin_taylor2(
      {4 * cpp_rational(Sin(u₀ / 4)),
       cpp_rational(Cos(u₀ / 4)),
       -cpp_rational(Sin(u₀ / 4)) / 8},
      u₀);
  AccuratePolynomial<cpp_rational, 2> cos_taylor2(
      {cpp_rational(Cos(u₀ / 4)),
       -cpp_rational(Sin(u₀ / 4) / 4),
       -cpp_rational(Cos(u₀ / 4) / 32)},
      u₀);

  auto const u = StehléZimmermannSimultaneousFullSearch<5>(
      {sin, cos},
      {sin_taylor2, cos_taylor2},
      u₀);
  EXPECT_THAT(u,
              IsOkAndHolds(cpp_rational(4785074604081885, 9007199254740992)));
  EXPECT_THAT(static_cast<double>(*u),
              RelativeErrorFrom(u₀, Lt(1.6e-13)));
  {
    std::string const mathematica = ToMathematica(sin(*u),
                                                  /*express_in=*/std::nullopt,
                                                  /*base=*/2);
    std::string_view mantissa = mathematica;
    CHECK(absl::ConsumePrefix(&mantissa, "Times[2^^"));
    EXPECT_THAT(mantissa.substr(53, 5), Eq("00000"));
  }
  {
    std::string const mathematica = ToMathematica(cos(*u),
                                                  /*express_in=*/std::nullopt,
                                                  /*base=*/2);
    std::string_view mantissa = mathematica;
    CHECK(absl::ConsumePrefix(&mantissa, "Times[2^^"));
    EXPECT_THAT(mantissa.substr(53, 5), Eq("00000"));
  }
}

TEST_F(AccurateTableGeneratorTest, StehléZimmermannFullSinCos15NoScaling) {
  double const x₀ = 17.0 / 128;
  double const u₀ = 4 * x₀;
  auto const sin = [](cpp_rational const& u) { return 4 * Sin(u / 4); };
  auto const cos = [](cpp_rational const& u) { return Cos(u / 4); };
  CHECK(0.5 <= u₀ & u₀ < 1.0) << u₀;
  CHECK(0.5 <= sin(u₀) && sin(u₀) < 1.0);
  CHECK(0.5 <= cos(u₀) && cos(u₀) < 1.0);
  AccuratePolynomial<cpp_rational, 2> sin_taylor2(
      {4 * cpp_rational(Sin(u₀ / 4)),
       cpp_rational(Cos(u₀ / 4)),
       -cpp_rational(Sin(u₀ / 4)) / 8},
      u₀);
  AccuratePolynomial<cpp_rational, 2> cos_taylor2(
      {cpp_rational(Cos(u₀ / 4)),
       -cpp_rational(Sin(u₀ / 4) / 4),
       -cpp_rational(Cos(u₀ / 4) / 32)},
      u₀);

  auto const u = StehléZimmermannSimultaneousFullSearch<15>(
      {sin, cos},
      {sin_taylor2, cos_taylor2},
      u₀);
  EXPECT_THAT(u,
              IsOkAndHolds(cpp_rational(4785074575333183, 9007199254740992)));
  EXPECT_THAT(static_cast<double>(*u),
              RelativeErrorFrom(u₀, Lt(6.1e-9)));
  {
    std::string const mathematica = ToMathematica(sin(*u),
                                                  /*express_in=*/std::nullopt,
                                                  /*base=*/2);
    std::string_view mantissa = mathematica;
    CHECK(absl::ConsumePrefix(&mantissa, "Times[2^^"));
    EXPECT_THAT(mantissa.substr(53, 15), Eq("00000""00000""00000"));
  }
  {
    std::string const mathematica = ToMathematica(cos(*u),
                                                  /*express_in=*/std::nullopt,
                                                  /*base=*/2);
    std::string_view mantissa = mathematica;
    CHECK(absl::ConsumePrefix(&mantissa, "Times[2^^"));
    EXPECT_THAT(mantissa.substr(53, 15), Eq("00000""00000""00000"));
  }
}

TEST_F(AccurateTableGeneratorTest, StehléZimmermannFullSinCos15WithScaling) {
  double const x₀ = 17.0 / 128;
  AccuratePolynomial<cpp_rational, 2> sin_taylor2({cpp_rational(Sin(x₀)),
                                                   cpp_rational(Cos(x₀)),
                                                   -cpp_rational(Sin(x₀)) / 2},
                                                  x₀);
  AccuratePolynomial<cpp_rational, 2> cos_taylor2({cpp_rational(Cos(x₀)),
                                                   -cpp_rational(Sin(x₀)),
                                                   -cpp_rational(Cos(x₀) / 2)},
                                                  x₀);

  auto const x = StehléZimmermannSimultaneousFullSearch<15>(
      {Sin, Cos},
      {sin_taylor2, cos_taylor2},
      x₀);
  EXPECT_THAT(x,
              IsOkAndHolds(cpp_rational(4785074575333183, 36028797018963968)));
  EXPECT_THAT(static_cast<double>(*x),
              RelativeErrorFrom(x₀, Lt(6.1e-9)));
  {
    std::string const mathematica = ToMathematica(Sin(*x),
                                                  /*express_in=*/std::nullopt,
                                                  /*base=*/2);
    std::string_view mantissa = mathematica;
    CHECK(absl::ConsumePrefix(&mantissa, "Times[2^^"));
    EXPECT_THAT(mantissa.substr(53, 15), Eq("00000""00000""00000"));
  }
  {
    std::string const mathematica = ToMathematica(Cos(*x),
                                                  /*express_in=*/std::nullopt,
                                                  /*base=*/2);
    std::string_view mantissa = mathematica;
    CHECK(absl::ConsumePrefix(&mantissa, "Times[2^^"));
    EXPECT_THAT(mantissa.substr(53, 15), Eq("00000""00000""00000"));
  }
}

TEST_F(AccurateTableGeneratorTest, StehléZimmermannMultisearchSinCos15) {
  static constexpr std::int64_t index_begin = 17;
  static constexpr std::int64_t index_end = 100;
  std::vector<cpp_rational> starting_arguments;
  std::vector<std::array<AccuratePolynomial<cpp_rational, 2>, 2>> polynomials;
  for (std::int64_t i = index_begin; i < index_end; ++i) {
    auto const x₀ = i / 128.0;
    AccuratePolynomial<cpp_rational, 2> const sin_taylor2(
        {cpp_rational(Sin(x₀)),
         cpp_rational(Cos(x₀)),
         -cpp_rational(Sin(x₀)) / 2},
        x₀);
    AccuratePolynomial<cpp_rational, 2> const cos_taylor2(
        {cpp_rational(Cos(x₀)),
         -cpp_rational(Sin(x₀)),
         -cpp_rational(Cos(x₀) / 2)},
        x₀);
    starting_arguments.push_back(x₀);
    polynomials.push_back({sin_taylor2, cos_taylor2});
  }
  auto const xs = StehléZimmermannSimultaneousMultisearch<15>(
      {Sin, Cos}, polynomials, starting_arguments);
  EXPECT_THAT(xs, SizeIs(index_end - index_begin));
  for (std::int64_t i = 0; i < xs.size(); ++i) {
    CHECK_OK(xs[i].status());
    auto const& x = *xs[i];
    EXPECT_THAT(static_cast<double>(x),
                RelativeErrorFrom((i + index_begin) / 128.0, Lt(1.3e-7)));
    {
      std::string const mathematica = ToMathematica(Sin(x),
                                                    /*express_in=*/std::nullopt,
                                                    /*base=*/2);
      std::string_view mantissa = mathematica;
      CHECK(absl::ConsumePrefix(&mantissa, "Times[2^^"));
      EXPECT_THAT(mantissa.substr(53, 15),
                  AnyOf(Eq("00000""00000""00000"),
                        Eq("11111""11111""11111")));
    }
    {
      std::string const mathematica = ToMathematica(Cos(x),
                                                    /*express_in=*/std::nullopt,
                                                    /*base=*/2);
      std::string_view mantissa = mathematica;
      CHECK(absl::ConsumePrefix(&mantissa, "Times[2^^"));
      EXPECT_THAT(mantissa.substr(53, 15),
                  AnyOf(Eq("00000""00000""00000"),
                        Eq("11111""11111""11111")));
    }
  }
}

TEST_F(AccurateTableGeneratorTest, DISABLED_SECULAR_SinCos18) {
  for (std::int64_t n = 0; n < 10; ++n) {
    Logger logger(TEMP_DIR / absl::StrCat("sin_cos_18_", n, ".wl"));

    // Process the binade [1 / 2^(n + 1), 1 / 2^n[ (except that for n = 0 the
    // upper bound is π / 4).
    double const lower_bound = std::asin(1.0 / (1 << (n + 1)));
    double const upper_bound = n == 0 ? π / 4 : std::asin(1.0 / (1 << n));

    double const h = 1.0 / (1 << (n + 10));
    double const h_over_2 = h / 2.0;

    std::vector<cpp_rational> starting_arguments;
    std::vector<std::array<AccuratePolynomial<cpp_rational, 2>, 2>> polynomials;
    for (std::int64_t i = std::floor(lower_bound / h_over_2);
         i <= std::ceil(upper_bound / h_over_2);
         ++i) {
      // The arguments are odd multiples of h/2.
      if (i % 2 == 1 && i == 1077) {
        double const x₀ = i * h_over_2;
        if (lower_bound <= x₀ && x₀ < upper_bound) {
          AccuratePolynomial<cpp_rational, 2> const sin_taylor2(
              {cpp_rational(Sin(x₀)),
               cpp_rational(Cos(x₀)),
               -cpp_rational(Sin(x₀)) / 2},
              x₀);
          AccuratePolynomial<cpp_rational, 2> const cos_taylor2(
              {cpp_rational(Cos(x₀)),
               -cpp_rational(Sin(x₀)),
               -cpp_rational(Cos(x₀) / 2)},
              x₀);
          starting_arguments.push_back(x₀);
          polynomials.push_back({sin_taylor2, cos_taylor2});
        }
      }
      if (i >= 1080) break;
    }

    auto const xs = StehléZimmermannSimultaneousMultisearch<18>(
        {Sin, Cos}, polynomials, starting_arguments);

    for (auto const& status_or_x : xs) {
      CHECK_OK(status_or_x.status());
      auto const& x = status_or_x.value();
      logger.Append(
          absl::StrCat("accurateTables[", n, "]"),
          std::tuple{static_cast<cpp_bin_float_50>(x), Sin(x), Cos(x)});
    }
  }
}

#endif

}  // namespace _accurate_table_generator
}  // namespace functions
}  // namespace principia
