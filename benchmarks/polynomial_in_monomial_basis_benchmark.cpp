// .\Release\x64\benchmarks.exe --benchmark_repetitions=3 --benchmark_filter=EvaluatePolynomialInMonomialBasis  // NOLINT(whitespace/line_length)

#include <random>
#include <tuple>

#include "astronomy/frames.hpp"
#include "benchmark/benchmark.h"
#include "benchmarks/metric.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/space.hpp"
#include "glog/logging.h"
#include "numerics/polynomial_evaluators.hpp"
#include "numerics/polynomial_in_monomial_basis.hpp"
#include "quantities/concepts.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {

using namespace principia::astronomy::_frames;
using namespace principia::benchmarks::_metric;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_r3_element;
using namespace principia::geometry::_space;
using namespace principia::numerics::_polynomial_evaluators;
using namespace principia::numerics::_polynomial_in_monomial_basis;
using namespace principia::quantities::_concepts;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

static constexpr std::int64_t number_of_iterations = 100;

template<typename T>
struct ValueGenerator;

template<>
struct ValueGenerator<double> {
  static double Get(std::mt19937_64& random) {
    std::uniform_real_distribution<> uniformly_at(-1.0, 1.0);
    return static_cast<double>(uniformly_at(random));
  }
};

template<typename D>
struct ValueGenerator<Quantity<D>> {
  static Quantity<D> Get(std::mt19937_64& random) {
    return ValueGenerator<double>::Get(random) * si::Unit<Quantity<D>>;
  }
};

template<typename S>
struct ValueGenerator<R3Element<S>> {
  static R3Element<S> Get(std::mt19937_64& random) {
    return {ValueGenerator<S>::Get(random),
            ValueGenerator<S>::Get(random),
            ValueGenerator<S>::Get(random)};
  }
};

template<typename S, typename F, int r>
struct ValueGenerator<Multivector<S, F, r>> {
  static Multivector<S, F, r> Get(std::mt19937_64& random) {
    return Multivector<S, F, r>(ValueGenerator<R3Element<S>>::Get(random));
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

// Constructs a polynomial and an argument so that the iteration performed for
// the latency benchmark does not generate nonfinite values.  Note that the
// iteration here must be identical to the one used in the benchmark.
template<typename Value, typename Argument, int degree,
         template<typename, typename, int> class Evaluator>
void PickPolynomialAndArgument(
    std::mt19937_64& random,
    PolynomialInMonomialBasis<Value, Argument, degree>& polynomial,
    Argument& initial_argument) {
  using P = PolynomialInMonomialBasis<Value, Argument, degree>;

  initial_argument = ValueGenerator<Argument>::Get(random);
  for (;;) {
    typename P::Coefficients coefficients;
    RandomTupleGenerator<typename P::Coefficients, 0>::Fill(coefficients,
                                                            random);
    polynomial = P(coefficients, with_evaluator<Evaluator>);
    auto argument = initial_argument;
    Value v;
    for (std::int64_t i = 0; i < number_of_iterations; ++i) {
      v = polynomial(argument);
      if constexpr (convertible_to_quantity<Value>) {
        argument = v * (Second / si::Unit<Value>);
      } else {
        argument = v.coordinates().x * (Second / Metre);
      }
    }
    if (IsFinite(argument)) {
      break;
    }
  }
}

template<typename Value, typename Argument, int degree, Metric metric,
         template<typename, typename, int> class Evaluator>
void EvaluatePolynomialInMonomialBasis(benchmark::State& state) {
  using P = PolynomialInMonomialBasis<Value, Argument, degree>;

  std::mt19937_64 random(42);
  P polynomial(typename P::Coefficients{});
  Argument argument;

  // We do this even for the throughput benchmark to make sure that we use the
  // same polynomial and argument in both cases.
  PickPolynomialAndArgument<Value, Argument, degree, Evaluator>(
      random, polynomial, argument);

  if constexpr (metric == Metric::Throughput) {
    Value v[number_of_iterations];
    while (state.KeepRunningBatch(number_of_iterations)) {
      for (std::int64_t i = 0; i < number_of_iterations; ++i) {
        v[i] = polynomial(argument);
      }
      benchmark::DoNotOptimize(v);
    }
  } else {
    static_assert(metric == Metric::Latency);
    Value v;
    while (state.KeepRunningBatch(number_of_iterations)) {
      for (std::int64_t i = 0; i < number_of_iterations ; ++i) {
        v = polynomial(argument);
        if constexpr (convertible_to_quantity<Value>) {
          argument = v * (Second / si::Unit<Value>);
        } else {
          argument = v.coordinates().x * (Second / Metre);
        }
      }
    }
  }
}

template<typename Value, Metric metric,
         template<typename, typename, int> typename Evaluator>
void BM_EvaluatePolynomialInMonomialBasis(benchmark::State& state) {
  int const degree = state.range(0);
  switch (degree) {
    case 2:
      EvaluatePolynomialInMonomialBasis<Value, Time, 2, metric, Evaluator>(
          state);
      break;
    case 4:
      EvaluatePolynomialInMonomialBasis<Value, Time, 4, metric, Evaluator>(
          state);
      break;
    case 6:
      EvaluatePolynomialInMonomialBasis<Value, Time, 6, metric, Evaluator>(
          state);
      break;
    case 8:
      EvaluatePolynomialInMonomialBasis<Value, Time, 8, metric, Evaluator>(
          state);
      break;
    case 10:
      EvaluatePolynomialInMonomialBasis<Value, Time, 10, metric, Evaluator>(
          state);
      break;
    case 12:
      EvaluatePolynomialInMonomialBasis<Value, Time, 12, metric, Evaluator>(
          state);
      break;
    case 14:
      EvaluatePolynomialInMonomialBasis<Value, Time, 14, metric, Evaluator>(
          state);
      break;
    case 16:
      EvaluatePolynomialInMonomialBasis<Value, Time, 16, metric, Evaluator>(
          state);
      break;
    default:
      LOG(FATAL) << "Degree " << degree
                 << " in BM_EvaluatePolynomialInMonomialBasis";
  }
}

BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   double, Metric::Latency, Estrin)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   double, Metric::Throughput, Estrin)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   double, Metric::Latency, Horner)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   double, Metric::Throughput, Horner)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   double, Metric::Latency, EstrinWithoutFMA)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   double, Metric::Throughput, EstrinWithoutFMA)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   double, Metric::Latency, HornerWithoutFMA)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   double, Metric::Throughput, HornerWithoutFMA)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   Length, Metric::Latency, Estrin)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   Length, Metric::Throughput, Estrin)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   Length, Metric::Latency, Horner)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   Length, Metric::Throughput, Horner)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   Length, Metric::Latency, EstrinWithoutFMA)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   Length, Metric::Throughput, EstrinWithoutFMA)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   Length, Metric::Latency, HornerWithoutFMA)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   Length, Metric::Throughput, HornerWithoutFMA)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   Displacement<ICRS>, Metric::Latency, Estrin)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   Displacement<ICRS>, Metric::Throughput, Estrin)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   Displacement<ICRS>, Metric::Latency, Horner)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   Displacement<ICRS>, Metric::Throughput, Horner)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   Displacement<ICRS>, Metric::Latency, EstrinWithoutFMA)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   Displacement<ICRS>, Metric::Throughput, EstrinWithoutFMA)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   Displacement<ICRS>, Metric::Latency, HornerWithoutFMA)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   Displacement<ICRS>, Metric::Throughput, HornerWithoutFMA)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);

}  // namespace numerics
}  // namespace principia
