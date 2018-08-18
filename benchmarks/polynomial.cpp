
#include <random>
#include <tuple>

#include "astronomy/frames.hpp"
#include "benchmark/benchmark.h"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/r3_element.hpp"
#include "numerics/polynomial.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

namespace principia {

using astronomy::ICRS;
using geometry::Displacement;
using geometry::Multivector;
using geometry::R3Element;
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
  using P = PolynomialInMonomialBasis<Value, Argument, degree, Evaluator>;
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

template<template<typename, typename, int> class Evaluator>
void BM_EvaluatePolynomialInMonomialBasisDouble(benchmark::State& state) {
  int const degree = state.range_x();
  switch (degree) {
    case 4:
      EvaluatePolynomialInMonomialBasis<double, Time, 4, Evaluator>(state);
      break;
    case 8:
      EvaluatePolynomialInMonomialBasis<double, Time, 8, Evaluator>(state);
      break;
    case 12:
      EvaluatePolynomialInMonomialBasis<double, Time, 12, Evaluator>(state);
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
  int const degree = state.range_x();
  switch (degree) {
    case 4:
      EvaluatePolynomialInMonomialBasis<Length, Time, 4, Evaluator>(state);
      break;
    case 8:
      EvaluatePolynomialInMonomialBasis<Length, Time, 8, Evaluator>(state);
      break;
    case 12:
      EvaluatePolynomialInMonomialBasis<Length, Time, 12, Evaluator>(state);
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
  int const degree = state.range_x();
  switch (degree) {
    case 4:
      EvaluatePolynomialInMonomialBasis<Multivector<double, ICRS, 1>,
                                        Time,
                                        4,
                                        Evaluator>(state);
      break;
    case 8:
      EvaluatePolynomialInMonomialBasis<Multivector<double, ICRS, 1>,
                                        Time,
                                        8,
                                        Evaluator>(state);
      break;
    case 12:
      EvaluatePolynomialInMonomialBasis<Multivector<double, ICRS, 1>,
                                        Time,
                                        12,
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
  int const degree = state.range_x();
  switch (degree) {
    case 4:
      EvaluatePolynomialInMonomialBasis<Displacement<ICRS>,
                                        Time,
                                        4,
                                        Evaluator>(state);
      break;
    case 8:
      EvaluatePolynomialInMonomialBasis<Displacement<ICRS>,
                                        Time,
                                        8,
                                        Evaluator>(state);
      break;
    case 12:
      EvaluatePolynomialInMonomialBasis<Displacement<ICRS>,
                                        Time,
                                        12,
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

BENCHMARK_TEMPLATE1(BM_EvaluatePolynomialInMonomialBasisDouble,
                    EstrinEvaluator)
    ->Arg(4)->Arg(8)->Arg(12)->Arg(16);
BENCHMARK_TEMPLATE1(BM_EvaluatePolynomialInMonomialBasisQuantity,
                    EstrinEvaluator)
    ->Arg(4)->Arg(8)->Arg(12)->Arg(16);
BENCHMARK_TEMPLATE1(BM_EvaluatePolynomialInMonomialBasisVectorDouble,
                    EstrinEvaluator)
    ->Arg(4)->Arg(8)->Arg(12)->Arg(16);
BENCHMARK_TEMPLATE1(BM_EvaluatePolynomialInMonomialBasisDisplacement,
                    EstrinEvaluator)
    ->Arg(4)->Arg(8)->Arg(12)->Arg(16);
BENCHMARK_TEMPLATE1(BM_EvaluatePolynomialInMonomialBasisDouble,
                    HornerEvaluator)
    ->Arg(4)->Arg(8)->Arg(12)->Arg(16);
BENCHMARK_TEMPLATE1(BM_EvaluatePolynomialInMonomialBasisQuantity,
                    HornerEvaluator)
    ->Arg(4)->Arg(8)->Arg(12)->Arg(16);
BENCHMARK_TEMPLATE1(BM_EvaluatePolynomialInMonomialBasisVectorDouble,
                    HornerEvaluator)
    ->Arg(4)->Arg(8)->Arg(12)->Arg(16);
BENCHMARK_TEMPLATE1(BM_EvaluatePolynomialInMonomialBasisDisplacement,
                    HornerEvaluator)
    ->Arg(4)->Arg(8)->Arg(12)->Arg(16);

}  // namespace numerics
}  // namespace principia
