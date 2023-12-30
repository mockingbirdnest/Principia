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
using ::testing::Gt;
using ::testing::Lt;

TEST(ApproximationTest, SinInverse) {
  auto const f = [](double const x) { return Sin(1 * Radian / x); };
  auto const approximation =
      ЧебышёвPolynomialInterpolant<double>(f,
                                           /*lower_bound=*/0.01,
                                           /*upper_bound=*/10,
                                           /*max_error=*/1e-6);
  EXPECT_EQ(32, approximation.degree());
  for (double x = 0.01; x < 10.0; x += 0.01) {
    EXPECT_THAT(approximation.Evaluate(x),
                AbsoluteErrorFrom(f(x), AllOf(Lt(1.9), Gt(3.7e-14))));
  }
}

TEST(ApproximationTest, Exp) {
  auto const f = [](double const x) { return std::exp(x); };
  auto const approximation =
      ЧебышёвPolynomialInterpolant<double>(f,
                                           /*lower_bound=*/0.01,
                                           /*upper_bound=*/3,
                                           /*max_error=*/1e-6);
  EXPECT_EQ(32, approximation.degree());
  for (double x = 0.01; x < 3; x += 0.01) {
    EXPECT_THAT(approximation.Evaluate(x),
                AbsoluteErrorFrom(f(x), AllOf(Lt(1.9), Gt(3.7e-14))));
  }
}

}  // namespace _approximation
}  // namespace numerics
}  // namespace principia
