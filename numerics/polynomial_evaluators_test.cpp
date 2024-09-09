#include "numerics/polynomial_evaluators.hpp"

#include <cstddef>
#include <utility>

#include "gtest/gtest.h"
#include "numerics/combinatorics.hpp"

namespace principia {
namespace numerics {

using namespace principia::numerics::_combinatorics;
using namespace principia::numerics::_polynomial_evaluators;

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
           int degree>
  void Test() {
    using E = Evaluator<double, double, degree>;
    auto const binomial_coefficients =
        MakeBinomialTuple<typename E::Coefficients, degree>(
            std::make_index_sequence<degree + 1>());
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
  Test<Estrin, 1>();
  Test<Estrin, 2>();
  Test<Estrin, 3>();
  Test<Estrin, 4>();
  Test<Estrin, 5>();
  Test<Estrin, 6>();
  Test<Estrin, 7>();
  Test<Estrin, 8>();
  Test<Estrin, 9>();
  Test<Estrin, 10>();
  Test<Estrin, 11>();
  Test<Estrin, 12>();
  Test<Estrin, 13>();
  Test<Estrin, 14>();
}

TEST_F(PolynomialEvaluatorTest, Horner) {
  Test<Horner, 1>();
  Test<Horner, 2>();
  Test<Horner, 3>();
  Test<Horner, 4>();
  Test<Horner, 5>();
  Test<Horner, 6>();
  Test<Horner, 7>();
  Test<Horner, 8>();
  Test<Horner, 9>();
  Test<Horner, 10>();
  Test<Horner, 11>();
  Test<Horner, 12>();
  Test<Horner, 13>();
  Test<Horner, 14>();
}

}  // namespace numerics
}  // namespace principia
