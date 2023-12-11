#include "numerics/approximation.hpp"

#include "gtest/gtest.h"
#include "quantities/elementary_functions.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {
namespace _approximation {

using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_si;

TEST(ApproximationTest, SinInverse) {
  std::function<double(double)> const f = [](double const x) -> double {
    return Sin(1 * Radian / x);
  };
  auto const approximation = ЧебышёвPolynomialInterpolant(f,
                                                          /*lower_bound=*/0.01,
                                                          /*upper_bound=*/10.0,
                                                          /*max_error=*/1e-6);
}

}  // namespace _approximation
}  // namespace numerics
}  // namespace principia
