#include "numerics/approximation.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "numerics/чебышёв_series.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/numerics_matchers.hpp"

namespace principia {
namespace numerics {
namespace _approximation {

using ::testing::AllOf;
using ::testing::ElementsAre;
using ::testing::Ge;
using ::testing::Lt;
using ::testing::Property;
using ::testing::SizeIs;
using namespace principia::numerics::_чебышёв_series;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_numerics_matchers;

TEST(ApproximationTest, SinInverse) {
  auto const f = [](double const x) { return Sin(1 * Radian / x); };
  TerminationPredicate<double, double> const done =
      [](auto const& _, double const& error_estimate) -> bool {
    return error_estimate < 1e-6;
  };
  double error_estimate;
  auto const interpolant =
      ЧебышёвPolynomialInterpolant<128>(f,
                                        /*lower_bound=*/0.1,
                                        /*upper_bound=*/10.0,
                                        done,
                                        &error_estimate);
  EXPECT_EQ(128, interpolant.degree());
  // Didn't reach the desired accuracy.
  EXPECT_THAT(error_estimate, IsNear(4.9e-6_(1)));
  for (double x = 0.1; x < 10.0; x += 0.01) {
    EXPECT_THAT(interpolant.Evaluate(x),
                AbsoluteErrorFrom(f(x), AllOf(Lt(3.0e-6), Ge(0))));
  }
}

TEST(ApproximationTest, Exp) {
  auto const f = [](double const x) { return std::exp(x); };
  TerminationPredicate<double, double> const done =
      [](auto const& _, double const& error_estimate) -> bool {
    return error_estimate < 1e-6;
  };
  double error_estimate;
  auto const interpolant =
      ЧебышёвPolynomialInterpolant<128>(f,
                                        /*lower_bound=*/0.01,
                                        /*upper_bound=*/3.0,
                                        done,
                                        &error_estimate);
  EXPECT_EQ(16, interpolant.degree());
  EXPECT_THAT(error_estimate, IsNear(4.7e-14_(1)));
  for (double x = 0.01; x < 3; x += 0.01) {
    EXPECT_THAT(interpolant.Evaluate(x),
                AbsoluteErrorFrom(f(x), AllOf(Lt(7.2e-15), Ge(0))));
  }
}

TEST(ApproximationTest, AdaptiveSinInverse) {
  auto const f = [](double const x) { return Sin(1 * Radian / x); };
  TerminationPredicate<double, double> const done =
      [](auto const& _, double const& error_estimate) -> bool {
    return error_estimate < 1e-6;
  };
  double error_estimate;
  auto const interpolants =
      AdaptiveЧебышёвPolynomialInterpolant<8>(f,
                                              /*lower_bound=*/0.1,
                                              /*upper_bound=*/10.0,
                                              done,
                                              &error_estimate);
  EXPECT_THAT(error_estimate, IsNear(7.1e-7_(1)));
  EXPECT_THAT(interpolants, SizeIs(11));
  EXPECT_THAT(
      interpolants,
      ElementsAre(AllOf(Property(&ЧебышёвSeries<double, double>::lower_bound,
                                 AlmostEquals(0.1, 0)),
                        Property(&ЧебышёвSeries<double, double>::upper_bound,
                                 AlmostEquals(0.1 + 9.9 / 512, 0)),
                        Property(&ЧебышёвSeries<double, double>::degree, 8)),
                  AllOf(Property(&ЧебышёвSeries<double, double>::lower_bound,
                                 AlmostEquals(0.1 + 9.9 / 512, 0)),
                        Property(&ЧебышёвSeries<double, double>::upper_bound,
                                 AlmostEquals(0.1 + 9.9 / 256, 0)),
                        Property(&ЧебышёвSeries<double, double>::degree, 8)),
                  AllOf(Property(&ЧебышёвSeries<double, double>::lower_bound,
                                 AlmostEquals(0.1 + 9.9 / 256, 0)),
                        Property(&ЧебышёвSeries<double, double>::upper_bound,
                                 AlmostEquals(0.1 + 9.9 / 128, 0)),
                        Property(&ЧебышёвSeries<double, double>::degree, 8)),
                  AllOf(Property(&ЧебышёвSeries<double, double>::lower_bound,
                                 AlmostEquals(0.1 + 9.9 / 128, 0)),
                        Property(&ЧебышёвSeries<double, double>::upper_bound,
                                 AlmostEquals(0.1 + 9.9 / 64, 0)),
                        Property(&ЧебышёвSeries<double, double>::degree, 8)),
                  AllOf(Property(&ЧебышёвSeries<double, double>::lower_bound,
                                 AlmostEquals(0.1 + 9.9 / 64, 0)),
                        Property(&ЧебышёвSeries<double, double>::upper_bound,
                                 AlmostEquals(0.1 + 9.9 / 32, 1)),
                        Property(&ЧебышёвSeries<double, double>::degree, 8)),
                  AllOf(Property(&ЧебышёвSeries<double, double>::lower_bound,
                                 AlmostEquals(0.1 + 9.9 / 32, 1)),
                        Property(&ЧебышёвSeries<double, double>::upper_bound,
                                 AlmostEquals(0.1 + 9.9 * 3 / 64, 0)),
                        Property(&ЧебышёвSeries<double, double>::degree, 8)),
                  AllOf(Property(&ЧебышёвSeries<double, double>::lower_bound,
                                 AlmostEquals(0.1 + 9.9 * 3 / 64, 0)),
                        Property(&ЧебышёвSeries<double, double>::upper_bound,
                                 AlmostEquals(0.1 + 9.9 / 16, 0)),
                        Property(&ЧебышёвSeries<double, double>::degree, 8)),
                  AllOf(Property(&ЧебышёвSeries<double, double>::lower_bound,
                                 AlmostEquals(0.1 + 9.9 / 16, 0)),
                        Property(&ЧебышёвSeries<double, double>::upper_bound,
                                 AlmostEquals(0.1 + 9.9 / 8, 1)),
                        Property(&ЧебышёвSeries<double, double>::degree, 8)),
                  AllOf(Property(&ЧебышёвSeries<double, double>::lower_bound,
                                 AlmostEquals(0.1 + 9.9 / 8, 1)),
                        Property(&ЧебышёвSeries<double, double>::upper_bound,
                                 AlmostEquals(0.1 + 9.9 / 4, 1)),
                        Property(&ЧебышёвSeries<double, double>::degree, 8)),
                  AllOf(Property(&ЧебышёвSeries<double, double>::lower_bound,
                                 AlmostEquals(0.1 + 9.9 / 4, 1)),
                        Property(&ЧебышёвSeries<double, double>::upper_bound,
                                 AlmostEquals(0.1 + 9.9 / 2, 0)),
                        Property(&ЧебышёвSeries<double, double>::degree, 8)),
                  AllOf(Property(&ЧебышёвSeries<double, double>::lower_bound,
                                 AlmostEquals(0.1 + 9.9 / 2, 0)),
                        Property(&ЧебышёвSeries<double, double>::upper_bound,
                                 AlmostEquals(10.0, 0)),
                        Property(&ЧебышёвSeries<double, double>::degree, 8))));
  for (auto const& interpolant : interpolants) {
    for (double x = interpolant.lower_bound();
         x < interpolant.upper_bound();
         x += 0.01) {
      EXPECT_THAT(interpolant.Evaluate(x),
                  AbsoluteErrorFrom(f(x), AllOf(Lt(6.2e-7), Ge(2.7e-17))));
    }
  }
}

}  // namespace _approximation
}  // namespace numerics
}  // namespace principia
