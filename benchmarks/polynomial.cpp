
#include <random>
#include <tuple>

#include "benchmark/benchmark.h"
#include "numerics/polynomial.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {

using quantities::Length;
using quantities::Quantity;
using quantities::SIUnit;
using quantities::Time;

namespace numerics {

namespace {
constexpr int evaluations_per_iteration = 1000;
}  // namespace

template<typename T>
struct ValueGenerator;

template<>
struct ValueGenerator<double> {
  static double Get(std::mt19937_64& random) {
    return static_cast<double>(random());
  }
};

template<typename D>
struct ValueGenerator<Quantity<D>> {
  static Quantity<D> Get(std::mt19937_64& random) {
    return static_cast<double>(random()) * SIUnit<Quantity<D>>();
  }
};

template<typename Tuple, int k, int size = std::tuple_size_v<Tuple>>
struct RandomTupleGenerator {
  static void Fill(Tuple& t, std::mt19937_64& random) {
    std::get<k>(t) =
        ValueGenerator<std::tuple_element_t<k, Tuple>>::Get(random);
    RandomTupleGenerator<Tuple, k + 1, size>::Fill(t, random);
  }
};

template<typename Tuple, int size>
struct RandomTupleGenerator<Tuple, size, size> {
  static void Fill(Tuple& t, std::mt19937_64& random) {}
};

template<typename Value, typename Argument, int degree>
void EvaluatePolynomialInMonomialBasis(benchmark::State& state) {
  using P = PolynomialInMonomialBasis<Value, Argument, degree>;
  std::mt19937_64 random(42);
  P::Coefficients coefficients;
  RandomTupleGenerator<P::Coefficients, 0>::Fill(coefficients, random);
  P const p(coefficients);

  auto const min = ValueGenerator<Argument>::Get(random);
  auto const max = ValueGenerator<Argument>::Get(random);
  auto argument = min;
  auto const Δargument = (max - min) * 1e-9;
  auto result = Value{};

  while (state.KeepRunning()) {
    for (int i = 0; i < evaluations_per_iteration; ++i) {
      result += p.Evaluate(argument);
      argument += Δargument;
    }
  }

  // This weird call to |SetLabel| has no effect except that it uses |result|
  // and therefore prevents the loop from being optimized away.
  std::stringstream ss;
  ss << result;
  state.SetLabel(ss.str().substr(0, 0));
}

void BM_EvaluatePolynomialInMonomialBasisDouble(benchmark::State& state) {
  int const degree = state.range_x();
  switch (degree) {
    case 4:
      EvaluatePolynomialInMonomialBasis<double, Time, 4>(state);
      break;
    case 8:
      EvaluatePolynomialInMonomialBasis<double, Time, 8>(state);
      break;
    case 16:
      EvaluatePolynomialInMonomialBasis<double, Time, 16>(state);
      break;
    default:
      LOG(FATAL) << "Degree " << degree << " in BM_EvaluatePolynomialDouble";
  }
}

void BM_EvaluatePolynomialInMonomialBasisQuantity(benchmark::State& state) {
  int const degree = state.range_x();
  switch (degree) {
    case 4:
      EvaluatePolynomialInMonomialBasis<Length, Time, 4>(state);
      break;
    case 8:
      EvaluatePolynomialInMonomialBasis<Length, Time, 8>(state);
      break;
    case 16:
      EvaluatePolynomialInMonomialBasis<Length, Time, 16>(state);
      break;
    default:
      LOG(FATAL) << "Degree " << degree << " in BM_EvaluatePolynomialQuantity";
  }
}

BENCHMARK(BM_EvaluatePolynomialInMonomialBasisDouble)->Arg(4)->Arg(8)->Arg(16);
BENCHMARK(BM_EvaluatePolynomialInMonomialBasisQuantity)
    ->Arg(4)->Arg(8)->Arg(16);

}  // namespace numerics
}  // namespace principia
