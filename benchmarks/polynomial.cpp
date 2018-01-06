
#include <random>
#include <tuple>

#include "benchmark/benchmark.h"
#include "numerics/polynomial.hpp"

namespace principia {
namespace numerics {

namespace {
int const evaluations_per_iteration = 1000;
}  // namespace

template<typename T>
T GenerateElement(std::mt19937_64& random);

template<>
double GenerateElement(std::mt19937_64& random) {
  return static_cast<double>(random());
}

template<typename Tuple, int k, int size = std::tuple_size_v<Tuple>>
struct RandomTupleGenerator {
  static void Fill(Tuple& t, std::mt19937_64& random) {
    std::get<k>(t) = GenerateElement<std::tuple_element_t<k, Tuple>>(random);
    RandomTupleGenerator<Tuple, k + 1, size>::Fill(t, random);
  }
};

template<typename Tuple, int size>
struct RandomTupleGenerator<Tuple, size, size> {
  static void Fill(Tuple& t, std::mt19937_64& random) {}
};

template<typename Value, typename Argument, typename degree>
void EvaluatePolynomialInMonomialBasis(benchmark::State& state) {
  using P = PolynomialInMonomialBasis<Value, Argument, degree>;
  std::mt19937_64 random(42);
  P::Coefficients coefficients;
  RandomTupleGenerator<P::Coefficients, 0>::Fill(coefficients, random);
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

void BM_EvaluatePolynomial(benchmark::State& state) {
  using P4 = PolynomialInMonomialBasis<double, double, 4>;
  std::mt19937_64 random(42);
  P4::Coefficients coefficients;
  RandomTupleGenerator<P4::Coefficients, 0>::Fill(coefficients, random);
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
