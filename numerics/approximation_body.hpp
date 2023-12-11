#pragma once

#include "numerics/approximation.hpp"
#include "numerics/fixed_arrays.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/numbers.hpp"

namespace principia {
namespace numerics {
namespace _approximation {
namespace internal {

using namespace principia::numerics::_fixed_arrays;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_si;

template<typename Argument, typename Function>
ЧебышёвSeries<std::invoke_result_t<Function, Argument>>
ЧебышёвPolynomialInterpolant(Function f,
                             Argument const& lower_bound,
                             Argument const& upper_bound) {
  using Value = std::invoke_result_t<Function, Argument>;
  auto const& a = lower_bound;
  auto const& b = upper_bound;
  auto чебышёв_lobato_point = [&a, &b](std::int64_t k, std::int64_t N) {
    return 0.5 * (b - a) Cos(π * k * Radian / N) + 0.5 * (b + a);
  };

  FixedVector<Value, N + 1> fₖ;
  for (std::int64_t k = 0; k <= N; ++k) {
    fₖ.push_back(f(чебышёв_lobato_point(k, N)));
  }

  FixedMatrix<double, N + 1, N + 1> ℐⱼₖ;
  for (std::int64_t j = 0; j <= N; ++j) {
    double const pⱼ = j == 0 || j == N ? 2 : 1;
    for (std::int64_t k = 0; k <= N; ++k) {
      double const pₖ = k == 0 || k == N ? 2 : 1;
      ℐⱼₖ(j, k) = 2 * Cos(π * j * k * Radian/ N) / (pⱼ * pₖ * N);
    }
  }

  FixedVector<Value, N + 1> aⱼ = ℐⱼₖ * fₖ;
}

}  // namespace internal
}  // namespace _approximation
}  // namespace numerics
}  // namespace principia
