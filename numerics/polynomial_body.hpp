#pragma once

#include "numerics/polynomial.hpp"

namespace principia {
namespace numerics {
namespace internal_polynomial {

template<typename Value, typename Argument>
Argument const& Polynomial<Value, Argument>::argument_min() const {
  return argument_min_;
}

template<typename Value, typename Argument>
Argument const& Polynomial<Value, Argument>::argument_max() const {
  return argument_max_;
}

template<typename Value, typename Argument>
Polynomial<Value, Argument>::Polynomial(Argument const& argument_min,
                                        Argument const& argument_max)
    : argument_min_(argument_min),
      argument_max_(argument_max) {}

template<typename Value, typename Argument, int degree>
PolynomialInMonomialBasis<Value, Argument, degree>::PolynomialInMonomialBasis(
    Coefficients const& coeffients,
    Argument const& argument_min,
    Argument const& argument_max)
    : Polynomial(argument_min, argument_max),
      coefficients_(coefficients) {}

template<typename Value, typename Argument, int degree>
Value PolynomialInMonomialBasis<Value, Argument, degree>::Evaluate(
    Argument const& argument) const {
}

template<typename Value, typename Argument, int degree>
Derivative<Value, Argument>
PolynomialInMonomialBasis<Value, Argument, degree>::EvaluateDerivative(
    Argument const& argument) const {
}

}  // namespace internal_polynomial
}  // namespace numerics
}  // namespace principia
