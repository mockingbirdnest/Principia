
#include "numerics/polynomial_evaluators.hpp"

#include <utility>

#include "gtest/gtest.h"
#include "numerics/combinatorics.hpp"

namespace principia {
namespace numerics {

class PolynomialEvaluatorTest : public ::testing::Test {
 public:
  template<typename Tuple, int n, int... k>
  Tuple MakeBinomialTuple(std::integer_sequence<int, k...>) {
    return {Binomial(n, k)...};
  }

  // This test builds the binomial polynomial (1 + x)^degree and evaluates it
  // using the |Evaluator| and directly using |std::pow|.
  template<template<typename Value, typename Argument, int degree>
           class Evaluator,
           int degree>
  void Test() {
    using E = Evaluator<double, double, degree>;
    auto const binomial_coefficients =
        MakeBinomialTuple<typename E::Coefficients, degree>(
            std::make_integer_sequence<int, degree + 1>());
    for (int argument = -degree; argument <= degree; ++argument) {
      EXPECT_EQ(E::Evaluate(binomial_coefficients, argument),
                std::pow(argument + 1, degree)) << argument << " " << degree;
      EXPECT_EQ(E::EvaluateDerivative(binomial_coefficients, argument),
                degree * std::pow(argument + 1, degree - 1))
          << argument << " " << degree;
    }
  }
};

TEST_F(PolynomialEvaluatorTest, Estrin) {
  Test<EstrinEvaluator, 1>();
  Test<EstrinEvaluator, 2>();
  Test<EstrinEvaluator, 3>();
  Test<EstrinEvaluator, 4>();
  Test<EstrinEvaluator, 5>();
  Test<EstrinEvaluator, 6>();
  Test<EstrinEvaluator, 7>();
  Test<EstrinEvaluator, 8>();
  Test<EstrinEvaluator, 9>();
  Test<EstrinEvaluator, 10>();
  Test<EstrinEvaluator, 11>();
  Test<EstrinEvaluator, 12>();
  Test<EstrinEvaluator, 13>();
  Test<EstrinEvaluator, 14>();
}

TEST_F(PolynomialEvaluatorTest, Horner) {
  Test<HornerEvaluator, 1>();
  Test<HornerEvaluator, 2>();
  Test<HornerEvaluator, 3>();
  Test<HornerEvaluator, 4>();
  Test<HornerEvaluator, 5>();
  Test<HornerEvaluator, 6>();
  Test<HornerEvaluator, 7>();
  Test<HornerEvaluator, 8>();
  Test<HornerEvaluator, 9>();
  Test<HornerEvaluator, 10>();
  Test<HornerEvaluator, 11>();
  Test<HornerEvaluator, 12>();
  Test<HornerEvaluator, 13>();
  Test<HornerEvaluator, 14>();
}

}  // namespace numerics
}  // namespace principia
