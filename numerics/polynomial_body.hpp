#pragma once

#include "numerics/polynomial.hpp"

#include "base/macros.hpp"
#include "numerics/polynomial_evaluators.hpp"

namespace principia {
namespace numerics {
namespace internal_polynomial {

template<typename Value, typename Argument, int degree,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Value, Argument, degree, Evaluator>::
PolynomialInMonomialBasis(Coefficients const& coefficients)
    : coefficients_(coefficients) {}

template<typename Value, typename Argument, int degree,
         template<typename, typename, int> class Evaluator>
Value PolynomialInMonomialBasis<Value, Argument, degree, Evaluator>::
Evaluate(Argument const& argument) const {
  return Evaluator<Value, Argument, degree>::Evaluate(coefficients_, argument);
}

template<typename Value, typename Argument, int degree,
         template<typename, typename, int> class Evaluator>
Derivative<Value, Argument>
PolynomialInMonomialBasis<Value, Argument, degree, Evaluator>::
EvaluateDerivative(Argument const& argument) const {
  return HornerEvaluator<Value, Argument, degree>::EvaluateDerivative(
      coefficients_, argument);
}

}  // namespace internal_polynomial
}  // namespace numerics
}  // namespace principia
