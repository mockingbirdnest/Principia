
#include <random>

#include "benchmark/benchmark.h"
#include "numerics/polynomial.hpp"

namespace principia {
namespace numerics {

namespace {
int const evaluations_per_iteration = 1000;
}  // namespace

void BM_EvaluatePolynomial(benchmark::State& state) {
  using P4 = PolynomialInMonomialBasis<double, double, 4>;
  std::mt19937_64 random(42);
  P4::Coefficients const coefficients = {static_cast<double>(random()),
                                         static_cast<double>(random()),
                                         static_cast<double>(random()),
                                         static_cast<double>(random()),
                                         static_cast<double>(random())};
  double const min = static_cast<double>(random());
  double const max = static_cast<double>(random());
  P4 const p4(coefficients, min, max);

  double argument = min;
  double const Δargument = (max - min) * 1e-9;
  double result = 0.0;

  while (state.KeepRunning()) {
    for (int i = 0; i < evaluations_per_iteration; ++i) {
      result += p4.Evaluate(argument);
      argument += Δargument;
    }
  }

  // This weird call to |SetLabel| has no effect except that it uses |result|
  // and therefore prevents the loop from being optimized away.
  state.SetLabel(std::to_string(result).substr(0, 0));
}

BENCHMARK(BM_EvaluatePolynomial);

}  // namespace numerics
}  // namespace principia
