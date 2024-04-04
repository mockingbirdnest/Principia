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
    return static_cast<double>(random()) * si::Unit<Quantity<D>>;
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

template<typename Value, typename Argument, int degree,
         template<typename, typename, int> class Evaluator>
void EvaluatePolynomialInMonomialBasis(benchmark::State& state) {
  using P = PolynomialInMonomialBasis<Value, Argument, degree>;
  std::mt19937_64 random(42);
  typename P::Coefficients coefficients;
  RandomTupleGenerator<typename P::Coefficients, 0>::Fill(coefficients, random);
  P const p(coefficients, with_evaluator<Evaluator>);

  auto const min = ValueGenerator<Argument>::Get(random);
  auto const max = ValueGenerator<Argument>::Get(random);
  auto argument = min;
  auto const Δargument = (max - min) * 1e-9;

#if 0
  while (state.KeepRunningBatch(10)) {
    auto const a0 = p(argument);
    //argument += Δargument;
    auto const a1 = p(argument);
    //argument += Δargument;
    auto const a2 = p(argument);
    //argument += Δargument;
    auto const a3 = p(argument);
    //argument += Δargument;
    auto const a4 = p(argument);
    //argument += Δargument;
    auto const a5 = p(argument);
    //argument += Δargument;
    auto const a6 = p(argument);
    //argument += Δargument;
    auto const a7 = p(argument);
    //argument += Δargument;
    auto const a8 = p(argument);
    //argument += Δargument;
    auto const a9 = p(argument);
    auto const s = a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9;
    benchmark::DoNotOptimize(s);
  }
#else
  using H = Hilbert<Value>;
  auto const scale = 1 / si::Unit<typename H::Norm²Type>;
  while (state.KeepRunningBatch(10)) {
#if 0
    auto const a0 = p(argument);
    argument += H::Norm²(a0) * scale * Δargument;
    auto const a1 = p(argument);
    argument += H::Norm²(a1) * scale * Δargument;
    auto const a2 = p(argument);
    argument += H::Norm²(a2) * scale * Δargument;
    auto const a3 = p(argument);
    argument += H::Norm²(a3) * scale * Δargument;
    auto const a4 = p(argument);
    argument += H::Norm²(a4) * scale * Δargument;
    auto const a5 = p(argument);
    argument += H::Norm²(a5) * scale * Δargument;
    auto const a6 = p(argument);
    argument += H::Norm²(a6) * scale * Δargument;
    auto const a7 = p(argument);
    argument += H::Norm²(a7) * scale * Δargument;
    auto const a8 = p(argument);
    argument += H::Norm²(a8) * scale * Δargument;
    auto const a9 = p(argument);
#elif 1
    auto const a0 = p(argument);
    benchmark::DoNotOptimize(a0);
    //argument += Δargument;
    auto const a1 = p(argument);
    benchmark::DoNotOptimize(a1);
    //argument += Δargument;
    auto const a2 = p(argument);
    benchmark::DoNotOptimize(a2);
    //argument += Δargument;
    auto const a3 = p(argument);
    benchmark::DoNotOptimize(a3);
    //argument += Δargument;
    auto const a4 = p(argument);
    benchmark::DoNotOptimize(a4);
    //argument += Δargument;
    auto const a5 = p(argument);
    benchmark::DoNotOptimize(a5);
    //argument += Δargument;
    auto const a6 = p(argument);
    benchmark::DoNotOptimize(a6);
    //argument += Δargument;
    auto const a7 = p(argument);
    benchmark::DoNotOptimize(a7);
    //argument += Δargument;
    auto const a8 = p(argument);
    benchmark::DoNotOptimize(a8);
    //argument += Δargument;
    auto const a9 = p(argument);
    benchmark::DoNotOptimize(a9);
#else
    auto const a0 = p(argument);
    argument += *reinterpret_cast<double const*>(&a0) * Δargument;
    auto const a1 = p(argument);
    argument += *reinterpret_cast<double const*>(&a1) * Δargument;
    auto const a2 = p(argument);
    argument += *reinterpret_cast<double const*>(&a2) * Δargument;
    auto const a3 = p(argument);
    argument += *reinterpret_cast<double const*>(&a3) * Δargument;
    auto const a4 = p(argument);
    argument += *reinterpret_cast<double const*>(&a4) * Δargument;
    auto const a5 = p(argument);
    argument += *reinterpret_cast<double const*>(&a5) * Δargument;
    auto const a6 = p(argument);
    argument += *reinterpret_cast<double const*>(&a6) * Δargument;
    auto const a7 = p(argument);
    argument += *reinterpret_cast<double const*>(&a7) * Δargument;
    auto const a8 = p(argument);
    argument += *reinterpret_cast<double const*>(&a8) * Δargument;
    auto const a9 = p(argument);
    benchmark::DoNotOptimize(a9);
#endif
  }

#endif
}

template<typename Value,
         template<typename, typename, int> class Evaluator>
void BM_EvaluatePolynomialInMonomialBasis(benchmark::State& state) {
  int const degree = state.range(0);
  switch (degree) {
    case 2:
      EvaluatePolynomialInMonomialBasis<Value, Time, 2, Evaluator>(state);
      break;
    case 4:
      EvaluatePolynomialInMonomialBasis<Value, Time, 4, Evaluator>(state);
      break;
    case 6:
      EvaluatePolynomialInMonomialBasis<Value, Time, 6, Evaluator>(state);
      break;
    case 8:
      EvaluatePolynomialInMonomialBasis<Value, Time, 8, Evaluator>(state);
      break;
    case 10:
      EvaluatePolynomialInMonomialBasis<Value, Time, 10, Evaluator>(state);
      break;
    case 12:
      EvaluatePolynomialInMonomialBasis<Value, Time, 12, Evaluator>(state);
      break;
    case 14:
      EvaluatePolynomialInMonomialBasis<Value, Time, 14, Evaluator>(state);
      break;
    case 16:
      EvaluatePolynomialInMonomialBasis<Value, Time, 16, Evaluator>(state);
      break;
    default:
      LOG(FATAL) << "Degree " << degree
                 << " in BM_EvaluatePolynomialInMonomialBasisDouble";
  }
}

BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   double, Estrin)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   Length, Estrin)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   Displacement<ICRS>, Estrin)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   double, Horner)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   Length, Horner)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   Displacement<ICRS>, Horner)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   double, EstrinWithoutFMA)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   Length, EstrinWithoutFMA)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   Displacement<ICRS>, EstrinWithoutFMA)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);

BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   double, HornerWithoutFMA)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   Length, HornerWithoutFMA)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);
BENCHMARK_TEMPLATE(BM_EvaluatePolynomialInMonomialBasis,
                   Displacement<ICRS>, HornerWithoutFMA)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kNanosecond);

}  // namespace numerics
}  // namespace principia
