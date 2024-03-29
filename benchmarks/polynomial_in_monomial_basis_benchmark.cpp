// .\Release\x64\benchmarks.exe --benchmark_repetitions=3 --benchmark_filter=EvaluatePolynomial  // NOLINT(whitespace/line_length)

#include <random>
#include <tuple>

#include "astronomy/frames.hpp"
#include "benchmark/benchmark.h"
#include "geometry/grassmann.hpp"
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
using namespace principia::geometry::_r3_element;
using namespace principia::geometry::_space;
using namespace principia::numerics::_polynomial_evaluators;
using namespace principia::numerics::_polynomial_in_monomial_basis;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

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

  for (auto _ : state) {
    for (int i = 0; i < evaluations_per_iteration; ++i) {
      benchmark::DoNotOptimize(p(argument));
      argument += Δargument;
    }
  }
}

template<template<typename, typename, int> class Evaluator>
void BM_EvaluatePolynomialInMonomialBasisDouble(benchmark::State& state) {
  int const degree = state.range(0);
  switch (degree) {
    case 2:
      EvaluatePolynomialInMonomialBasis<double, Time, 2, Evaluator>(state);
      break;
    case 4:
      EvaluatePolynomialInMonomialBasis<double, Time, 4, Evaluator>(state);
      break;
    case 6:
      EvaluatePolynomialInMonomialBasis<double, Time, 6, Evaluator>(state);
      break;
    case 8:
      EvaluatePolynomialInMonomialBasis<double, Time, 8, Evaluator>(state);
      break;
    case 10:
      EvaluatePolynomialInMonomialBasis<double, Time, 10, Evaluator>(state);
      break;
    case 12:
      EvaluatePolynomialInMonomialBasis<double, Time, 12, Evaluator>(state);
      break;
    case 14:
      EvaluatePolynomialInMonomialBasis<double, Time, 14, Evaluator>(state);
      break;
    case 16:
      EvaluatePolynomialInMonomialBasis<double, Time, 16, Evaluator>(state);
      break;
    default:
      LOG(FATAL) << "Degree " << degree
                 << " in BM_EvaluatePolynomialInMonomialBasisDouble";
  }
}

template<template<typename, typename, int> class Evaluator>
void BM_EvaluatePolynomialInMonomialBasisQuantity(benchmark::State& state) {
  int const degree = state.range(0);
  switch (degree) {
    case 2:
      EvaluatePolynomialInMonomialBasis<Length, Time, 2, Evaluator>(state);
      break;
    case 4:
      EvaluatePolynomialInMonomialBasis<Length, Time, 4, Evaluator>(state);
      break;
    case 6:
      EvaluatePolynomialInMonomialBasis<Length, Time, 6, Evaluator>(state);
      break;
    case 8:
      EvaluatePolynomialInMonomialBasis<Length, Time, 8, Evaluator>(state);
      break;
    case 10:
      EvaluatePolynomialInMonomialBasis<Length, Time, 10, Evaluator>(state);
      break;
    case 12:
      EvaluatePolynomialInMonomialBasis<Length, Time, 12, Evaluator>(state);
      break;
    case 14:
      EvaluatePolynomialInMonomialBasis<Length, Time, 14, Evaluator>(state);
      break;
    case 16:
      EvaluatePolynomialInMonomialBasis<Length, Time, 16, Evaluator>(state);
      break;
    default:
      LOG(FATAL) << "Degree " << degree
                 << " in BM_EvaluatePolynomialInMonomialBasisQuantity";
  }
}

template<template<typename, typename, int> class Evaluator>
void BM_EvaluatePolynomialInMonomialBasisVectorDouble(benchmark::State& state) {
  int const degree = state.range(0);
  switch (degree) {
    case 2:
      EvaluatePolynomialInMonomialBasis<Multivector<double, ICRS, 1>,
                                        Time,
                                        2,
                                        Evaluator>(state);
      break;
    case 4:
      EvaluatePolynomialInMonomialBasis<Multivector<double, ICRS, 1>,
                                        Time,
                                        4,
                                        Evaluator>(state);
      break;
    case 6:
      EvaluatePolynomialInMonomialBasis<Multivector<double, ICRS, 1>,
                                        Time,
                                        6,
                                        Evaluator>(state);
      break;
    case 8:
      EvaluatePolynomialInMonomialBasis<Multivector<double, ICRS, 1>,
                                        Time,
                                        8,
                                        Evaluator>(state);
      break;
    case 10:
      EvaluatePolynomialInMonomialBasis<Multivector<double, ICRS, 1>,
                                        Time,
                                        10,
                                        Evaluator>(state);
      break;
    case 12:
      EvaluatePolynomialInMonomialBasis<Multivector<double, ICRS, 1>,
                                        Time,
                                        12,
                                        Evaluator>(state);
      break;
    case 14:
      EvaluatePolynomialInMonomialBasis<Multivector<double, ICRS, 1>,
                                        Time,
                                        14,
                                        Evaluator>(state);
      break;
    case 16:
      EvaluatePolynomialInMonomialBasis<Multivector<double, ICRS, 1>,
                                        Time,
                                        16,
                                        Evaluator>(state);
      break;
    default:
      LOG(FATAL) << "Degree " << degree
                 << " in BM_EvaluatePolynomialInMonomialBasisVectorDouble";
  }
}

template<template<typename, typename, int> class Evaluator>
void BM_EvaluatePolynomialInMonomialBasisDisplacement(benchmark::State& state) {
  int const degree = state.range(0);
  switch (degree) {
    case 2:
      EvaluatePolynomialInMonomialBasis<Displacement<ICRS>,
                                        Time,
                                        2,
                                        Evaluator>(state);
      break;
    case 4:
      EvaluatePolynomialInMonomialBasis<Displacement<ICRS>,
                                        Time,
                                        4,
                                        Evaluator>(state);
      break;
    case 6:
      EvaluatePolynomialInMonomialBasis<Displacement<ICRS>,
                                        Time,
                                        6,
                                        Evaluator>(state);
      break;
    case 8:
      EvaluatePolynomialInMonomialBasis<Displacement<ICRS>,
                                        Time,
                                        8,
                                        Evaluator>(state);
      break;
    case 10:
      EvaluatePolynomialInMonomialBasis<Displacement<ICRS>,
                                        Time,
                                        10,
                                        Evaluator>(state);
      break;
    case 12:
      EvaluatePolynomialInMonomialBasis<Displacement<ICRS>,
                                        Time,
                                        12,
                                        Evaluator>(state);
      break;
    case 14:
      EvaluatePolynomialInMonomialBasis<Displacement<ICRS>,
                                        Time,
                                        14,
                                        Evaluator>(state);
      break;
    case 16:
      EvaluatePolynomialInMonomialBasis<Displacement<ICRS>,
                                        Time,
                                        16,
                                        Evaluator>(state);
      break;
    default:
      LOG(FATAL) << "Degree " << degree
                 << " in BM_EvaluatePolynomialInMonomialBasisDisplacement";
  }
}

BENCHMARK_TEMPLATE1(BM_EvaluatePolynomialInMonomialBasisDouble, Estrin)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE1(BM_EvaluatePolynomialInMonomialBasisQuantity, Estrin)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE1(BM_EvaluatePolynomialInMonomialBasisVectorDouble, Estrin)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE1(BM_EvaluatePolynomialInMonomialBasisDisplacement, Estrin)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE1(BM_EvaluatePolynomialInMonomialBasisDouble, Horner)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE1(BM_EvaluatePolynomialInMonomialBasisQuantity, Horner)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE1(BM_EvaluatePolynomialInMonomialBasisVectorDouble, Horner)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE1(BM_EvaluatePolynomialInMonomialBasisDisplacement, Horner)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE1(BM_EvaluatePolynomialInMonomialBasisDouble,
                    EstrinWithoutFMA)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE1(BM_EvaluatePolynomialInMonomialBasisQuantity,
                    EstrinWithoutFMA)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE1(BM_EvaluatePolynomialInMonomialBasisVectorDouble,
                    EstrinWithoutFMA)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE1(BM_EvaluatePolynomialInMonomialBasisDisplacement,
                    EstrinWithoutFMA)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kMicrosecond);

BENCHMARK_TEMPLATE1(BM_EvaluatePolynomialInMonomialBasisDouble,
                    HornerWithoutFMA)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE1(BM_EvaluatePolynomialInMonomialBasisQuantity,
                    HornerWithoutFMA)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE1(BM_EvaluatePolynomialInMonomialBasisVectorDouble,
                    HornerWithoutFMA)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE1(BM_EvaluatePolynomialInMonomialBasisDisplacement,
                    HornerWithoutFMA)
    ->Arg(2)->Arg(4)->Arg(6)->Arg(8)->Arg(10)->Arg(12)->Arg(14)->Arg(16)
    ->Unit(benchmark::kMicrosecond);

}  // namespace numerics
}  // namespace principia
