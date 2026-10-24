#include "numerics/polynomial_evaluators.hpp"

#include <cmath>
#include <cstddef>
#include <utility>

#include "gtest/gtest.h"
#include "numerics/combinatorics.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace numerics {

using namespace principia::numerics::_combinatorics;
using namespace principia::numerics::_polynomial_evaluators;
using namespace principia::testing_utilities::_almost_equals;

class PolynomialEvaluatorTest : public ::testing::Test {
 public:
  template<typename Tuple, int n, std::size_t... k>
  Tuple MakeBinomialTuple(std::index_sequence<k...>) {
    return {Binomial(n, k)...};
  }

  // This test builds the binomial polynomial (1 + x)^degree and evaluates it
  // using the `Evaluator` and directly using `std::pow`.
  template<template<typename Value, typename Argument, int degree>
           typename Evaluator,
           int degree,
           int min_ulps,
           int max_ulps>
  void Test() {
    using C =
        _polynomial_evaluators::internal::Coefficients<double, double, degree>;
    using E = Evaluator<double, double, degree>;
    auto const binomial_coefficients =
        MakeBinomialTuple<C, degree>(std::make_index_sequence<degree + 1>());
    double value;
    double derivative;
    for (int argument = -degree; argument <= degree; ++argument) {
      EXPECT_THAT(
          E::Evaluate(binomial_coefficients, argument),
          AlmostEquals(std::pow(argument + 1, degree), min_ulps, max_ulps))
          << argument << " " << degree;
      EXPECT_THAT(
          E::EvaluateDerivative(binomial_coefficients, argument),
          AlmostEquals(
              degree * std::pow(argument + 1, degree - 1), min_ulps, max_ulps))
          << argument << " " << degree;
      E::EvaluateWithDerivative(
          binomial_coefficients, argument, value, derivative);
      EXPECT_THAT(
          value,
          AlmostEquals(std::pow(argument + 1, degree), min_ulps, max_ulps))
          << argument << " " << degree;
      EXPECT_THAT(
          derivative,
          AlmostEquals(
              degree * std::pow(argument + 1, degree - 1), min_ulps, max_ulps))
          << argument << " " << degree;
    }
  }
};

TEST_F(PolynomialEvaluatorTest, Estrin) {
  Test<Estrin, 0, 0, 0>();
  Test<Estrin, 1, 0, 0>();
  Test<Estrin, 2, 0, 0>();
  Test<Estrin, 3, 0, 0>();
  Test<Estrin, 4, 0, 0>();
  Test<Estrin, 5, 0, 0>();
  Test<Estrin, 6, 0, 0>();
  Test<Estrin, 7, 0, 0>();
  Test<Estrin, 8, 0, 0>();
  Test<Estrin, 9, 0, 0>();
  Test<Estrin, 10, 0, 0>();
  Test<Estrin, 11, 0, 0>();
  Test<Estrin, 12, 0, 0>();
  Test<Estrin, 13, 0, 0>();
  Test<Estrin, 14, 0, 0>();
  Test<Estrin, 15, 0, 0>();
  Test<Estrin, 16, 0, 2>();
  Test<Estrin, 17, 0, 4>();
}

TEST_F(PolynomialEvaluatorTest, Horner) {
  Test<Horner, 0, 0, 0>();
  Test<Horner, 1, 0, 0>();
  Test<Horner, 2, 0, 0>();
  Test<Horner, 3, 0, 0>();
  Test<Horner, 4, 0, 0>();
  Test<Horner, 5, 0, 0>();
  Test<Horner, 6, 0, 0>();
  Test<Horner, 7, 0, 0>();
  Test<Horner, 8, 0, 0>();
  Test<Horner, 9, 0, 0>();
  Test<Horner, 10, 0, 0>();
  Test<Horner, 11, 0, 0>();
  Test<Horner, 12, 0, 0>();
  Test<Horner, 13, 0, 0>();
  Test<Horner, 14, 0, 0>();
  Test<Horner, 15, 0, 1>();
  Test<Horner, 16, 0, 3>();
  Test<Horner, 17, 0, 1>();
}

}  // namespace numerics
}  // namespace principia
