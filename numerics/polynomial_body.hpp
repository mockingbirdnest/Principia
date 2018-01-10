#pragma once

#include "numerics/polynomial.hpp"

#include "base/macros.hpp"
#include "numerics/polynomial_evaluators.hpp"

namespace principia {
namespace numerics {
namespace internal_polynomial {

template<typename Value, typename Argument, int degree>
PolynomialInMonomialBasis<Value, Argument, degree>::PolynomialInMonomialBasis(
    Coefficients const& coefficients)
    : coefficients_(coefficients) {}

template<typename Value, typename Argument, int degree>
Value PolynomialInMonomialBasis<Value, Argument, degree>::Evaluate(
    Argument const& argument) const {
#if 0
  return HornerEvaluator<Value, Argument, degree>::Evaluate(
      coefficients_, argument);
#else
  return EstrinEvaluator<Value, Argument, degree>::Evaluate(
      coefficients_,argument);
#endif
}

template<typename Value, typename Argument, int degree>
Derivative<Value, Argument>
PolynomialInMonomialBasis<Value, Argument, degree>::EvaluateDerivative(
    Argument const& argument) const {
  return HornerEvaluator<Value, Argument, degree>::EvaluateDerivative(
      coefficients_, argument);
}

}  // namespace internal_polynomial
}  // namespace numerics
}  // namespace principia
