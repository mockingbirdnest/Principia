#include "functions/accurate_table_generator.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <string_view>
#include <vector>

#include "absl/strings/strip.h"
#include "boost/multiprecision/cpp_int.hpp"
#include "functions/multiprecision.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "mathematica/logger.hpp"
#include "mathematica/mathematica.hpp"
#include "numerics/combinatorics.hpp"
#include "quantities/numbers.hpp"  // 🧙 For π.
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
using namespace principia::numerics::_combinatorics;
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
  // that the `Sin` has zeroes in the right place we fiddle with the Mathematica
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

  // Use the Lagrange form of the remainder.  u may be above or below u₀.  Note
  // that we use the fact that the functions are monotonic.
  auto const remainder_sin_taylor2 = [u₀](cpp_rational const& u) {
    auto const Δu = static_cast<double>(u) - u₀;
    auto const Δu³ = Δu * Δu * Δu;
    return Δu³ * (-std::cos(std::min(u₀ + Δu, u₀) / 4) / 16) / Factorial(3);
  };
  auto const remainder_cos_taylor2 = [u₀](cpp_rational const& u) {
    auto const Δu = static_cast<double>(u) - u₀;
    auto const Δu³ = Δu * Δu * Δu;
    return Δu³ * (std::sin(std::max(u₀ + Δu, u₀) / 4) / 64) / Factorial(3);
  };

  auto const u = StehléZimmermannSimultaneousSearch<15>(
      {sin, cos},
      {sin_taylor2, cos_taylor2},
      {remainder_sin_taylor2, remainder_cos_taylor2},
      u₀,
      /*N=*/1ll << 53,
      /*T=*/1ll << 21);
  EXPECT_THAT(u,
              IsOkAndHolds(cpp_rational(4785074575333183, 9007199254740992)));
  EXPECT_EQ(*u, cpp_rational(static_cast<double>(*u)));
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
  auto const sin_taylor2 = [](cpp_rational const& u₀) {
    return AccuratePolynomial<cpp_rational, 2>(
        {4 * cpp_rational(Sin(u₀ / 4)),
         cpp_rational(Cos(u₀ / 4)),
         -cpp_rational(Sin(u₀ / 4)) / 8},
        u₀);
  };
  auto const cos_taylor2 = [](cpp_rational const& u₀) {
    return AccuratePolynomial<cpp_rational, 2>(
        {cpp_rational(Cos(u₀ / 4)),
         -cpp_rational(Sin(u₀ / 4) / 4),
         -cpp_rational(Cos(u₀ / 4) / 32)},
        u₀);
  };

  auto const remainder_sin_taylor2 = [](cpp_rational const& u₀) {
    return [u₀ = static_cast<double>(u₀)](cpp_rational const& u) {
      auto const Δu = static_cast<double>(u) - u₀;
      auto const Δu³ = Δu * Δu * Δu;
      return Δu³ * (-std::cos(std::min(u₀ + Δu, u₀) / 4) / 16) / Factorial(3);
    };
  };
  auto const remainder_cos_taylor2 = [](cpp_rational const& u₀) {
    return [u₀ = static_cast<double>(u₀)](cpp_rational const& u) {
      auto const Δu = static_cast<double>(u) - u₀;
      auto const Δu³ = Δu * Δu * Δu;
      return Δu³ * (std::sin(std::max(u₀ + Δu, u₀) / 4) / 64) / Factorial(3);
    };
  };

  auto const u = StehléZimmermannSimultaneousFullSearch<5>(
      {sin, cos},
      {sin_taylor2, cos_taylor2},
      {remainder_sin_taylor2, remainder_cos_taylor2},
      u₀);
  EXPECT_THAT(u,
              IsOkAndHolds(cpp_rational(1196268651020245, 2251799813685248)));
  EXPECT_EQ(*u, cpp_rational(static_cast<double>(*u)));
  EXPECT_THAT(static_cast<double>(*u),
              RelativeErrorFrom(u₀, Lt(3.7e-14)));
  {
    std::string const mathematica = ToMathematica(sin(*u),
                                                  /*express_in=*/std::nullopt,
                                                  /*base=*/2);
    std::string_view mantissa = mathematica;
    CHECK(absl::ConsumePrefix(&mantissa, "Times[2^^"));
    EXPECT_THAT(mantissa.substr(53, 5), Eq("11111"));
  }
  {
    std::string const mathematica = ToMathematica(cos(*u),
                                                  /*express_in=*/std::nullopt,
                                                  /*base=*/2);
    std::string_view mantissa = mathematica;
    CHECK(absl::ConsumePrefix(&mantissa, "Times[2^^"));
    EXPECT_THAT(mantissa.substr(53, 5), Eq("11111"));
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
  auto const sin_taylor2 = [](cpp_rational const& u₀) {
    return AccuratePolynomial<cpp_rational, 2>(
        {4 * cpp_rational(Sin(u₀ / 4)),
         cpp_rational(Cos(u₀ / 4)),
         -cpp_rational(Sin(u₀ / 4)) / 8},
        u₀);
  };
  auto const cos_taylor2 = [](cpp_rational const& u₀) {
    return AccuratePolynomial<cpp_rational, 2>(
        {cpp_rational(Cos(u₀ / 4)),
         -cpp_rational(Sin(u₀ / 4) / 4),
         -cpp_rational(Cos(u₀ / 4) / 32)},
        u₀);
  };

  auto const remainder_sin_taylor2 = [](cpp_rational const& u₀) {
    return [u₀ = static_cast<double>(u₀)](cpp_rational const& u) {
      auto const Δu = static_cast<double>(u) - u₀;
      auto const Δu³ = Δu * Δu * Δu;
      return Δu³ * (-std::cos(std::min(u₀ + Δu, u₀) / 4) / 16) / Factorial(3);
    };
  };
  auto const remainder_cos_taylor2 = [](cpp_rational const& u₀) {
    return [u₀ = static_cast<double>(u₀)](cpp_rational const& u) {
      auto const Δu = static_cast<double>(u) - u₀;
      auto const Δu³ = Δu * Δu * Δu;
      return Δu³ * (std::sin(std::max(u₀ + Δu, u₀) / 4) / 64) / Factorial(3);
    };
  };

  auto const u = StehléZimmermannSimultaneousFullSearch<15>(
      {sin, cos},
      {sin_taylor2, cos_taylor2},
      {remainder_sin_taylor2, remainder_cos_taylor2},
      u₀);
  EXPECT_THAT(u,
              IsOkAndHolds(cpp_rational(4785074575333183, 9007199254740992)));
  EXPECT_EQ(*u, cpp_rational(static_cast<double>(*u)));
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
  auto const sin_taylor2 = [](cpp_rational const& x₀) {
    return AccuratePolynomial<cpp_rational, 2>({cpp_rational(Sin(x₀)),
                                                cpp_rational(Cos(x₀)),
                                                -cpp_rational(Sin(x₀)) / 2},
                                               x₀);
  };
  auto const cos_taylor2 = [](cpp_rational const& x₀) {
    return AccuratePolynomial<cpp_rational, 2>({cpp_rational(Cos(x₀)),
                                                -cpp_rational(Sin(x₀)),
                                                -cpp_rational(Cos(x₀) / 2)},
                                               x₀);
  };

  auto const remainder_sin_taylor2 = [](cpp_rational const& x₀) {
    return [x₀ = static_cast<double>(x₀)](cpp_rational const& x) {
      auto const Δx = static_cast<double>(x) - x₀;
      auto const Δx³ = Δx * Δx * Δx;
      return -Δx³ * -std::cos(std::min(x₀ + Δx, x₀)) / Factorial(3);
    };
  };
  auto const remainder_cos_taylor2 = [](cpp_rational const& x₀) {
    return [x₀ = static_cast<double>(x₀)](cpp_rational const& x) {
      auto const Δx = static_cast<double>(x) - x₀;
      auto const Δx³ = Δx * Δx * Δx;
      return Δx³ * std::sin(std::max(x₀ + Δx, x₀)) / Factorial(3);
    };
  };

  auto const x = StehléZimmermannSimultaneousFullSearch<15>(
      {Sin, Cos},
      {sin_taylor2, cos_taylor2},
      {remainder_sin_taylor2, remainder_cos_taylor2},
      x₀);
  EXPECT_THAT(x,
              IsOkAndHolds(cpp_rational(4785074575333183, 36028797018963968)));
  EXPECT_EQ(*x, cpp_rational(static_cast<double>(*x)));
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
  std::vector<std::array<AccuratePolynomialFactory<cpp_rational, 2>, 2>>
      polynomials;
  std::vector<std::array<ApproximateFunctionFactory, 2>> remainders;
  for (std::int64_t i = index_begin; i < index_end; ++i) {
    auto const x₀ = i / 128.0;
    auto const sin_taylor2 = [](cpp_rational const& x₀) {
      return AccuratePolynomial<cpp_rational, 2>({cpp_rational(Sin(x₀)),
                                                  cpp_rational(Cos(x₀)),
                                                  -cpp_rational(Sin(x₀)) / 2},
                                                 x₀);
    };
    auto const cos_taylor2 = [](cpp_rational const& x₀) {
      return AccuratePolynomial<cpp_rational, 2>({cpp_rational(Cos(x₀)),
                                                  -cpp_rational(Sin(x₀)),
                                                  -cpp_rational(Cos(x₀) / 2)},
                                                 x₀);
    };

    auto const remainder_sin_taylor2 = [](cpp_rational const& x₀) {
      return [x₀ = static_cast<double>(x₀)](cpp_rational const& x) {
        auto const Δx = static_cast<double>(x) - x₀;
        auto const Δx³ = Δx * Δx * Δx;
        return -Δx³ * -std::cos(std::min(x₀ + Δx, x₀)) / Factorial(3);
      };
    };
    auto const remainder_cos_taylor2 = [](cpp_rational const& x₀) {
      return [x₀ = static_cast<double>(x₀)](cpp_rational const& x) {
        auto const Δx = static_cast<double>(x) - x₀;
        auto const Δx³ = Δx * Δx * Δx;
        return Δx³ * std::sin(std::max(x₀ + Δx, x₀)) / Factorial(3);
      };
    };

    starting_arguments.push_back(x₀);
    polynomials.push_back({sin_taylor2, cos_taylor2});
    remainders.push_back({remainder_sin_taylor2, remainder_cos_taylor2});
  }
  auto const xs = StehléZimmermannSimultaneousMultisearch<15>(
      {Sin, Cos}, polynomials, remainders, starting_arguments);
  EXPECT_THAT(xs, SizeIs(index_end - index_begin));
  for (std::int64_t i = 0; i < xs.size(); ++i) {
    CHECK_OK(xs[i].status());
    auto const& x = *xs[i];
    EXPECT_EQ(x, cpp_rational(static_cast<double>(x)));
    EXPECT_THAT(static_cast<double>(x),
                RelativeErrorFrom((i + index_begin) / 128.0, Lt(1.3e-7)))
        << starting_arguments[i];
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
  static constexpr std::int64_t bits = 18;
  Logger logger(TEMP_DIR / absl::StrCat("sin_cos_", bits, ".wl"),
                /*make_unique=*/false);

  // The radius of each interval.
  double const h = 1.0 / (1 << 10);

  // The centre of the interval with index `i`.
  auto const centre = [h](std::int64_t const i) { return 2 * i * h; };

  // The index of the first interval, which starts at `h` with a centre at
  // `2 * h`.
  std::int64_t const i_min = 1;

  // The index of the last interval, which goes a bit beyond π / 4.
  std::int64_t i_max = std::ceil(π / (8 * h) - 0.5);

  // Check that the last interval straddles π / 4.
  CHECK_LT(centre(i_max) - h, π / 4);
  CHECK_LT(π / 4, centre(i_max) + h);

  // No need for fancy angle reduction as the angles are small.
  AccurateFunction const accurate_sin = [](cpp_rational const& x) {
    return sin(static_cast<cpp_bin_float_50>(x));
  };
  AccurateFunction const accurate_cos = [](cpp_rational const& x) {
    return cos(static_cast<cpp_bin_float_50>(x));
  };

  std::vector<cpp_rational> starting_arguments;
  std::vector<std::array<AccuratePolynomialFactory<cpp_rational, 2>, 2>>
      polynomials;
  std::vector<std::array<ApproximateFunctionFactory, 2>> remainders;
  for (std::int64_t i = i_min; i <= i_max; ++i) {
    double const x₀ = centre(i);
    auto const sin_taylor2 = [](cpp_rational const& x₀) {
      return AccuratePolynomial<cpp_rational, 2>({cpp_rational(Sin(x₀)),
                                                  cpp_rational(Cos(x₀)),
                                                  -cpp_rational(Sin(x₀)) / 2},
                                                 x₀);
    };
    auto const cos_taylor2 = [](cpp_rational const& x₀) {
      return AccuratePolynomial<cpp_rational, 2>({cpp_rational(Cos(x₀)),
                                                  -cpp_rational(Sin(x₀)),
                                                  -cpp_rational(Cos(x₀) / 2)},
                                                 x₀);
    };

    auto const remainder_sin_taylor2 = [](cpp_rational const& x₀) {
      return [x₀ = static_cast<double>(x₀)](cpp_rational const& x) {
        auto const Δx = static_cast<double>(x) - x₀;
        auto const Δx³ = Δx * Δx * Δx;
        return -Δx³ * -std::cos(std::min(x₀ + Δx, x₀)) / Factorial(3);
      };
    };
    auto const remainder_cos_taylor2 = [](cpp_rational const& x₀) {
      return [x₀ = static_cast<double>(x₀)](cpp_rational const& x) {
        auto const Δx = static_cast<double>(x) - x₀;
        auto const Δx³ = Δx * Δx * Δx;
        return Δx³ * std::sin(std::max(x₀ + Δx, x₀)) / Factorial(3);
      };
    };

    starting_arguments.push_back(x₀);
    polynomials.push_back({sin_taylor2, cos_taylor2});
    remainders.push_back({remainder_sin_taylor2, remainder_cos_taylor2});
  }

  StehléZimmermannSimultaneousStreamingMultisearch<bits>(
      {accurate_sin, accurate_cos},
      polynomials,
      remainders,
      starting_arguments,
      [i_min, &logger](std::int64_t const index,
                       absl::StatusOr<cpp_rational> status_or_x) {
        auto const& x = status_or_x.value();
        auto const sin_x = Sin(x);
        auto const cos_x = Cos(x);
        {
          std::string const mathematica =
              ToMathematica(sin_x,
                            /*express_in=*/std::nullopt,
                            /*base=*/2);
          std::string_view mantissa = mathematica;
          CHECK(absl::ConsumePrefix(&mantissa, "Times[2^^"));
          EXPECT_THAT(mantissa.substr(53, bits),
                      AnyOf(Eq("00000""00000""00000""000"),
                            Eq("11111""11111""11111""111")));
        }
        {
          std::string const mathematica =
              ToMathematica(cos_x,
                            /*express_in=*/std::nullopt,
                            /*base=*/2);
          std::string_view mantissa = mathematica;
          CHECK(absl::ConsumePrefix(&mantissa, "Times[2^^"));
          EXPECT_THAT(mantissa.substr(53, bits),
                      AnyOf(Eq("00000""00000""00000""000"),
                            Eq("11111""11111""11111""111")));
        }
        logger.Set(absl::StrCat("accurateTables[", index + i_min, "]"),
                   std::tuple{static_cast<cpp_bin_float_50>(x), sin_x, cos_x});
        logger.FlushAndClear();
      });
}

#endif

}  // namespace _accurate_table_generator
}  // namespace functions
}  // namespace principia
