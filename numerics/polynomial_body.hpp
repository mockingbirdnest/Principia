#pragma once

#include "numerics/polynomial.hpp"

#include <memory>
#include <utility>

#include "absl/strings/str_cat.h"
#include "absl/strings/str_join.h"
#include "numerics/polynomial_in_monomial_basis.hpp"

namespace principia {
namespace numerics {
namespace _polynomial {
namespace internal {

using namespace principia::numerics::_polynomial_in_monomial_basis;

#define PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(value)                       \
  case value:                                                               \
    return make_not_null_unique<                                            \
        PolynomialInMonomialBasis<Value, Argument, value>>(                 \
        PolynomialInMonomialBasis<Value, Argument, value>::ReadFromMessage( \
            message))

#define PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE_WITH_EVALUATOR(value,         \
                                                              evaluator_tag) \
  case value:                                                                \
    return make_not_null_unique<                                             \
        PolynomialInMonomialBasis<Value, Argument, value>>(                  \
        PolynomialInMonomialBasis<Value, Argument, value>::ReadFromMessage(  \
            message, evaluator_tag))

template<typename Value_, typename Argument_>
not_null<std::unique_ptr<Polynomial<Value_, Argument_>>>
Polynomial<Value_, Argument_>::ReadFromMessage(
    serialization::Polynomial const& message) {
  // 24 is the largest exponent that we can serialize for Quantity.
  switch (message.degree()) {
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(1);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(2);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(3);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(4);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(5);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(6);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(7);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(8);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(9);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(10);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(11);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(12);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(13);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(14);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(15);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(16);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(17);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(18);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(19);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(20);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(21);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(22);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(23);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(24);
    default:
      LOG(FATAL) << "Unexpected degree " << message.degree();
      break;
  }
}

template<typename Value_, typename Argument_>
template<template<typename, typename, int> typename Evaluator_>
not_null<std::unique_ptr<Polynomial<Value_, Argument_>>>
Polynomial<Value_, Argument_>::ReadFromMessage(
    serialization::Polynomial const& message,
    with_evaluator_t<Evaluator_> evaluator_tag) {
  // 24 is the largest exponent that we can serialize for Quantity.
  switch (message.degree()) {
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE_WITH_EVALUATOR(1, evaluator_tag);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE_WITH_EVALUATOR(2, evaluator_tag);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE_WITH_EVALUATOR(3, evaluator_tag);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE_WITH_EVALUATOR(4, evaluator_tag);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE_WITH_EVALUATOR(5, evaluator_tag);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE_WITH_EVALUATOR(6, evaluator_tag);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE_WITH_EVALUATOR(7, evaluator_tag);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE_WITH_EVALUATOR(8, evaluator_tag);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE_WITH_EVALUATOR(9, evaluator_tag);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE_WITH_EVALUATOR(10, evaluator_tag);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE_WITH_EVALUATOR(11, evaluator_tag);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE_WITH_EVALUATOR(12, evaluator_tag);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE_WITH_EVALUATOR(13, evaluator_tag);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE_WITH_EVALUATOR(14, evaluator_tag);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE_WITH_EVALUATOR(15, evaluator_tag);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE_WITH_EVALUATOR(16, evaluator_tag);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE_WITH_EVALUATOR(17, evaluator_tag);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE_WITH_EVALUATOR(18, evaluator_tag);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE_WITH_EVALUATOR(19, evaluator_tag);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE_WITH_EVALUATOR(20, evaluator_tag);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE_WITH_EVALUATOR(21, evaluator_tag);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE_WITH_EVALUATOR(22, evaluator_tag);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE_WITH_EVALUATOR(23, evaluator_tag);
    PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE_WITH_EVALUATOR(24, evaluator_tag);
    default:
      LOG(FATAL) << "Unexpected degree " << message.degree();
      break;
  }
}

#undef PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE_WITH_EVALUATOR
#undef PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE

}  // namespace internal
}  // namespace _polynomial
}  // namespace numerics
}  // namespace principia
