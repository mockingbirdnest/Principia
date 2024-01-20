#include "numerics/approximation.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/numerics_matchers.hpp"

namespace principia {
namespace numerics {
namespace _approximation {

using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_numerics_matchers;
using ::testing::AllOf;
using ::testing::Ge;
using ::testing::Lt;

TEST(ApproximationTest, SinInverse) {
  auto const f = [](double const x) { return Sin(1 * Radian / x); };
  auto const approximation =
      ЧебышёвPolynomialInterpolant<128>(f,
                                        /*lower_bound=*/0.1,
                                        /*upper_bound=*/10.0,
                                        /*max_error=*/1e-6);
  EXPECT_EQ(128, approximation.degree());
  for (double x = 0.1; x < 10.0; x += 0.01) {
    EXPECT_THAT(approximation.Evaluate(x),
                AbsoluteErrorFrom(f(x), AllOf(Lt(3.0e-6), Ge(0))));
  }
}

TEST(ApproximationTest, Exp) {
  auto const f = [](double const x) { return std::exp(x); };
  auto const approximation =
      ЧебышёвPolynomialInterpolant<128>(f,
                                        /*lower_bound=*/0.01,
                                        /*upper_bound=*/3.0,
                                        /*max_error=*/1e-6);
  EXPECT_EQ(16, approximation.degree());
  for (double x = 0.01; x < 3; x += 0.01) {
    EXPECT_THAT(approximation.Evaluate(x),
                AbsoluteErrorFrom(f(x), AllOf(Lt(7.2e-15), Ge(0))));
  }
}

}  // namespace _approximation
}  // namespace numerics
}  // namespace principia
