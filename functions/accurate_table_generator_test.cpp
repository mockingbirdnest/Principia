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
#include "testing_utilities/matchers.hpp"
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
using namespace principia::testing_utilities::_matchers;
using namespace principia::testing_utilities::_numerics_matchers;

class AccurateTableGeneratorTest : public ::testing::Test {};

#if 1

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

TEST_F(AccurateTableGeneratorTest, GalSinCos5Multisearch) {
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

#if 0
TEST_F(AccurateTableGeneratorTest, StehléZimmermannSinCos15NotFound) {
  auto const u₀ =
      cpp_rational("160560481864624295660215/302231454903657293676544");
  auto const u_expected =
      cpp_rational("1284483861309301654012917/604462909807314587353088");
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
      /*T=*/1ll << 22);
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
#endif

TEST_F(AccurateTableGeneratorTest, IterativeStehléZimmermannSinCos15) {
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

  auto const u = IterativeStehléZimmermannSimultaneousSearch<15>(
      {sin, cos},
      {sin_taylor2, cos_taylor2},
      u₀,
      /*N=*/1ll << 53);
  //EXPECT_THAT(u,
  //            IsOkAndHolds(cpp_rational(4785074575333183, 9007199254740992)));
  //EXPECT_THAT(static_cast<double>(*u),
  //            RelativeErrorFrom(u₀, Lt(1.3e-10)));
  //{
  //  std::string const mathematica = ToMathematica(sin(*u),
  //                                                /*express_in=*/std::nullopt,
  //                                                /*base=*/2);
  //  std::string_view mantissa = mathematica;
  //  CHECK(absl::ConsumePrefix(&mantissa, "Times[2^^"));
  //  EXPECT_EQ("00000""00000""00000", mantissa.substr(53, 15));
  //}
  //{
  //  std::string const mathematica = ToMathematica(cos(*u),
  //                                                /*express_in=*/std::nullopt,
  //                                                /*base=*/2);
  //  std::string_view mantissa = mathematica;
  //  CHECK(absl::ConsumePrefix(&mantissa, "Times[2^^"));
  //  EXPECT_EQ("00000""00000""00000", mantissa.substr(53, 15));
  //}
}
#endif

}  // namespace _accurate_table_generator
}  // namespace functions
}  // namespace principia
