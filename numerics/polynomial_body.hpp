#pragma once

#include "numerics/polynomial.hpp"

#include <memory>
#include <utility>

#include "numerics/polynomial_evaluators.hpp"
#include "numerics/polynomial_in_monomial_basis.hpp"

namespace principia {
namespace numerics {
namespace _polynomial {
namespace internal {

using namespace principia::numerics::_polynomial_evaluators;
using namespace principia::numerics::_polynomial_in_monomial_basis;

template<affine Value_, affine Argument_>
  requires homogeneous_affine_module<Value_, Difference<Argument_>>
Value_ PRINCIPIA_VECTORCALL
Polynomial<Value_, Argument_>::operator()(Argument argument) const {
  return VirtualEvaluate(argument);
}

template<affine Value_, affine Argument_>
  requires homogeneous_affine_module<Value_, Difference<Argument_>>
Derivative<Value_, Argument_> PRINCIPIA_VECTORCALL
Polynomial<Value_, Argument_>::EvaluateDerivative(Argument argument) const {
  return VirtualEvaluateDerivative(argument);
}

template<affine Value_, affine Argument_>
  requires homogeneous_affine_module<Value_, Difference<Argument_>>
void PRINCIPIA_VECTORCALL Polynomial<Value_, Argument_>::EvaluateWithDerivative(
    Argument argument,
    Value& value,
    Derivative<Value, Argument>& derivative) const {
  VirtualEvaluateWithDerivative(argument, value, derivative);
}

template<affine Value_, affine Argument_>
  requires homogeneous_affine_module<Value_, Difference<Argument_>>
not_null<std::unique_ptr<Polynomial<Value_, Argument_>>>
Polynomial<Value_, Argument_>::ReadFromMessage(
    serialization::Polynomial const& message) {
  CHECK(
      message.HasExtension(serialization::PolynomialInMonomialBasis::extension))
      << "Polynomial serialization is only supported for "
         "PolynomialInMonomialBasis";
  auto const& extension =
      message.GetExtension(serialization::PolynomialInMonomialBasis::extension);
  CHECK(extension.has_evaluator())
      << "No evaluator specified for pre-Καραθεοδωρή deserialization "
      << extension.DebugString();
  switch (extension.evaluator().kind()) {
    case serialization::PolynomialInMonomialBasis::Evaluator::ESTRIN:
      return ReadFromMessage<Estrin>(message);
    case serialization::PolynomialInMonomialBasis::Evaluator::HORNER:
      return ReadFromMessage<Horner>(message);
    case serialization::PolynomialInMonomialBasis::Evaluator::
        ESTRIN_WITHOUT_FMA:
      return ReadFromMessage<EstrinWithoutFMA>(message);
    case serialization::PolynomialInMonomialBasis::Evaluator::
        HORNER_WITHOUT_FMA:
      return ReadFromMessage<HornerWithoutFMA>(message);
    default:
      LOG(FATAL) << "Unexpected evaluator " << message.DebugString();
  }
}

#define PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE(value)                  \
  case value:                                                          \
    return make_not_null_unique<                                       \
        PolynomialInMonomialBasis<Value, Argument, value, Evaluator>>( \
        PolynomialInMonomialBasis<Value, Argument, value, Evaluator>:: \
            ReadFromMessage(message))

template<affine Value_, affine Argument_>
  requires homogeneous_affine_module<Value_, Difference<Argument_>>
template<template<typename, typename, int> typename Evaluator>
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

#undef PRINCIPIA_POLYNOMIAL_DEGREE_VALUE_CASE

}  // namespace internal
}  // namespace _polynomial
}  // namespace numerics
}  // namespace principia
