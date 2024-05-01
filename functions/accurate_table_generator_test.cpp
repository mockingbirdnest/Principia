#include "functions/accurate_table_generator.hpp"

#include <string>
#include <string_view>
#include <vector>

#include "absl/strings/strip.h"
#include "boost/multiprecision/cpp_int.hpp"
#include "functions/multiprecision.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "mathematica/mathematica.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/numerics_matchers.hpp"

namespace principia {
namespace functions {
namespace _accurate_table_generator {

using ::testing::Lt;
using ::testing::SizeIs;
using namespace boost::multiprecision;
using namespace principia::functions::_multiprecision;
using namespace principia::mathematica::_mathematica;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_numerics_matchers;

class AccurateTableGeneratorTest : public ::testing::Test {};

#if 1

TEST_F(AccurateTableGeneratorTest, Sin5) {
  auto const x = ExhaustiveSearch<5>({Sin}, 5.0 / 128.0);
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

TEST_F(AccurateTableGeneratorTest, SinCos5) {
  auto const x = ExhaustiveSearch<5>({Sin, Cos}, 95.0 / 128.0);
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

TEST_F(AccurateTableGeneratorTest, SinCos5Multisearch) {
  static constexpr std::int64_t index_begin = 17;
  static constexpr std::int64_t index_end = 100;
  std::vector<cpp_rational> starting_arguments;
  for (std::int64_t i = index_begin; i < index_end; ++i) {
    starting_arguments.push_back(i / 128.0);
  }
  auto const xs = ExhaustiveMultisearch<5>({Sin, Cos}, starting_arguments);
  EXPECT_THAT(xs, SizeIs(index_end - index_begin));
  for (std::int64_t i = 0; i < xs.size(); ++i) {
    auto const& x = xs[i];
LOG(ERROR)<<i + index_begin<<" "<<x;
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

TEST_F(AccurateTableGeneratorTest, SinCos5BadCase) {
  double const x₀ = 17.0 / 128.0;
  AccuratePolynomial<2> sin_taylor2(
      AccuratePolynomial<2>::Coefficients{cpp_rational(4 * Sin(x₀)),
                                          cpp_rational(4 * Cos(x₀)),
                                          -0.5 * cpp_rational(4 * Sin(x₀))},
      x₀);
  AccuratePolynomial<2> cos_taylor2(
      AccuratePolynomial<2>::Coefficients{cpp_rational(Cos(x₀)),
                                          -cpp_rational(Sin(x₀)),
                                          -0.5 * cpp_rational(Cos(x₀))},
      x₀);
  auto const x = SimultaneousBadCaseSearch<5>(
      {[](cpp_rational const& x) { return 4 * Sin(x); }, Cos},
      {sin_taylor2, cos_taylor2},
      x₀,
      /*M=*/1ll << 5,
      /*N=*/1ll << 53,
      /*T=*/1ll << 10);
  LOG(ERROR)<<x.status();
}

#endif

}  // namespace _accurate_table_generator
}  // namespace functions
}  // namespace principia
