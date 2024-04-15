// .\Release\x64\benchmarks.exe --benchmark_repetitions=3 --benchmark_filter=EvaluatePolynomialInMonomialBasis  // NOLINT(whitespace/line_length)

#include <random>
#include <tuple>

#include "astronomy/frames.hpp"
#include "benchmark/benchmark.h"
#include "geometry/grassmann.hpp"
#include "geometry/hilbert.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/space.hpp"
#include "glog/logging.h"
#include "numerics/polynomial_evaluators.hpp"
#include "numerics/polynomial_in_monomial_basis.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {

using namespace principia::astronomy::_frames;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_hilbert;
using namespace principia::geometry::_r3_element;
using namespace principia::geometry::_space;
using namespace principia::numerics::_polynomial_evaluators;
using namespace principia::numerics::_polynomial_in_monomial_basis;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

enum class Metric {
  Latency,
  Throughput
};

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

template<typename Value, typename Argument, int degree, Metric metric,
         template<typename, typename, int> class Evaluator>
void EvaluatePolynomialInMonomialBasis(benchmark::State& state) {
  using P = PolynomialInMonomialBasis<Value, Argument, degree>;
  std::mt19937_64 random(42);
  typename P::Coefficients coefficients;
  RandomTupleGenerator<typename P::Coefficients, 0>::Fill(coefficients, random);
  P p(coefficients, with_evaluator<Evaluator>);
  auto const initial_argument = ValueGenerator<Argument>::Get(random);
  auto argument = initial_argument;
  do {
    RandomTupleGenerator<typename P::Coefficients, 0>::Fill(coefficients,
                                                            random);
    p = P(coefficients, with_evaluator<Evaluator>);
    argument = initial_argument;
    for (std::int64_t i = 0; i < 100; ++i) {
      Value v;
      v = p(argument);
      if constexpr (std::is_same_v<Value, double>) {
        argument = v * Second;
      } else {
        argument = v.coordinates().x * (Second / Metre);
      }
    }
  } while (!IsFinite(argument));

  argument = initial_argument;

  if constexpr (metric == Metric::Throughput) {
    Value v[100];
    while (state.KeepRunningBatch(100)) {
      for (std::int64_t i = 0; i < 100; ++i) {
        v[i] = p(argument);
      }
      benchmark::DoNotOptimize(v);
    }
  } else {
    Value v;
    while (state.KeepRunningBatch(100)) {
      for (std::int64_t i = 0; i < 100 ; ++i) {
        v = p(argument);
        if constexpr (std::is_same_v<Value, double>) {
          argument = v * Second;
        } else {
          argument = v.coordinates().x * (Second / Metre);
        }
      }
    }
  }
}

template<typename Value,
         Metric metric,
         template<typename, typename, int>
         class Evaluator>
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
                 << " in BM_EvaluatePolynomialInMonomialBasisDouble";
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

//BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
//                   Length, Estrin)
//    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
//    ->Unit(benchmark::kNanosecond);

//BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
//                   Length, Horner)
//    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
//    ->Unit(benchmark::kNanosecond);

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

//BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
//                   double, EstrinWithoutFMA)
//    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
//    ->Unit(benchmark::kNanosecond);
//BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
//                   Length, EstrinWithoutFMA)
//    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
//    ->Unit(benchmark::kNanosecond);
//BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
//                   Displacement<ICRS>, EstrinWithoutFMA)
//    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
//    ->Unit(benchmark::kNanosecond);
//
//BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
//                   double, HornerWithoutFMA)
//    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
//    ->Unit(benchmark::kNanosecond);
//BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
//                   Length, HornerWithoutFMA)
//    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
//    ->Unit(benchmark::kNanosecond);
//BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
//                   Displacement<ICRS>, HornerWithoutFMA)
//    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
//    ->Unit(benchmark::kNanosecond);

}  // namespace numerics
}  // namespace principia
