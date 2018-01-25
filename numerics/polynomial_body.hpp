#pragma once

#include "numerics/polynomial.hpp"

#include "base/macros.hpp"

namespace principia {
namespace numerics {
namespace internal_polynomial {

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator>::
PolynomialInMonomialBasis(Coefficients const& coefficients)
    : coefficients_(coefficients) {}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
Value PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator>::
Evaluate(Argument const& argument) const {
  return Evaluator<Value, Argument, degree_>::Evaluate(coefficients_, argument);
}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
Derivative<Value, Argument>
PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator>::
EvaluateDerivative(Argument const& argument) const {
  return Evaluator<Value, Argument, degree_>::EvaluateDerivative(
      coefficients_, argument);
}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
constexpr int
PolynomialInMonomialBasis<Value, Argument, degree_, Evaluator>::degree() const {
  return degree_;
}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
PolynomialInMonomialBasis<Value, Point<Argument>, degree_, Evaluator>::
PolynomialInMonomialBasis(Coefficients const& coefficients,
                          Point<Argument> const& origin)
    : coefficients_(coefficients),
      origin_(origin) {}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
Value PolynomialInMonomialBasis<Value, Point<Argument>, degree_, Evaluator>::
Evaluate(Point<Argument> const& argument) const {
  return Evaluator<Value, Argument, degree_>::Evaluate(
      coefficients_, argument - origin_);
}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
Derivative<Value, Argument>
PolynomialInMonomialBasis<Value, Point<Argument>, degree_, Evaluator>::
EvaluateDerivative(Point<Argument> const& argument) const {
  return Evaluator<Value, Argument, degree_>::EvaluateDerivative(
      coefficients_, argument - origin_);
}

template<typename Value, typename Argument, int degree_,
         template<typename, typename, int> class Evaluator>
constexpr int
PolynomialInMonomialBasis<Value, Point<Argument>, degree_, Evaluator>::
degree() const {
  return degree_;
}

}  // namespace internal_polynomial
}  // namespace numerics
}  // namespace principia
