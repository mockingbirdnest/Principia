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
//  using T1 = Difference<internal::Value<decltype(f), double>>;
  using T2 = Difference<std::invoke_result_t<decltype(f), double>>;
  // auto const approximation =
  //     ЧебышёвPolynomialInterpolant<double>(f,
  //                                          /*lower_bound=*/0.01,
  //                                          /*upper_bound=*/10,
  //                                          /*max_error=*/1e-6);
}

}  // namespace _approximation
}  // namespace numerics
}  // namespace principia
