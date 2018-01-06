#pragma once

#include "numerics/polynomial.hpp"

#include "base/macros.hpp"

namespace principia {
namespace numerics {
namespace internal_polynomial {

namespace {

// The structs below are only needed because the language won't let us have (1)
// a partial specialization of a function (2) an explicit specialization of a
// member function template in an unspecialized class template.  Sigh.
// We use FORCE_INLINE because we have to write this recursively, but we really
// want linear code.

template<typename Value, typename Argument, int degree, int k>
struct HornerEvaluator {
  using Coefficients =
      typename PolynomialInMonomialBasis<Value, Argument, degree>::Coefficients;

  static FORCE_INLINE NthDerivative<Value, Argument, k>
  Evaluate(Coefficients const& coefficients,
           Difference<Argument> const& argument);
  static FORCE_INLINE NthDerivative<Value, Argument, k>
  EvaluateDerivative(Coefficients const& coefficients,
                     Difference<Argument> const& argument);
};

template<typename Value, typename Argument, int degree>
struct HornerEvaluator<Value, Argument, degree, degree> {
  using Coefficients =
      typename PolynomialInMonomialBasis<Value, Argument, degree>::Coefficients;

  static FORCE_INLINE NthDerivative<Value, Argument, degree>
  Evaluate(Coefficients const& coefficients,
           Difference<Argument> const& argument);
  static FORCE_INLINE NthDerivative<Value, Argument, degree>
  EvaluateDerivative(Coefficients const& coefficients,
                     Difference<Argument> const& argument);
};

template<typename Value, typename Argument, int degree, int k>
NthDerivative<Value, Argument, k>
HornerEvaluator<Value, Argument, degree, k>::Evaluate(
    Coefficients const& coefficients,
    Difference<Argument> const& argument) {
  return std::get<k>(coefficients) +
         argument * 
         HornerEvaluator<Value, Argument, degree, k + 1>::Evaluate(
             coefficients, argument);
}

template<typename Value, typename Argument, int degree, int k>
NthDerivative<Value, Argument, k>
HornerEvaluator<Value, Argument, degree, k>::EvaluateDerivative(
    Coefficients const& coefficients,
    Difference<Argument> const& argument) {
  return std::get<k>(coefficients) * k +
         argument *
         HornerEvaluator<Value, Argument, degree, k + 1>::EvaluateDerivative(
             coefficients, argument);
}

template<typename Value, typename Argument, int degree>
NthDerivative<Value, Argument, degree>
HornerEvaluator<Value, Argument, degree, degree>::Evaluate(
    Coefficients const& coefficients,
    Difference<Argument> const& argument) {
  return std::get<degree>(coefficients);
}

template<typename Value, typename Argument, int degree>
NthDerivative<Value, Argument, degree>
HornerEvaluator<Value, Argument, degree, degree>::EvaluateDerivative(
    Coefficients const& coefficients,
    Difference<Argument> const& argument) {
  return std::get<degree>(coefficients) * degree;
}

}  // namespace

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
    Coefficients const& coefficients,
    Argument const& argument_min,
    Argument const& argument_max)
    : Polynomial<Value, Argument>(argument_min, argument_max),
      coefficients_(coefficients) {}

template<typename Value, typename Argument, int degree>
Value PolynomialInMonomialBasis<Value, Argument, degree>::Evaluate(
    Argument const& argument) const {
  return HornerEvaluator<Value, Argument, degree, 0>::Evaluate(
      coefficients_, argument - argument_min());
}

template<typename Value, typename Argument, int degree>
Derivative<Value, Argument>
PolynomialInMonomialBasis<Value, Argument, degree>::EvaluateDerivative(
    Argument const& argument) const {
  return HornerEvaluator<Value, Argument, degree, 1>::EvaluateDerivative(
      coefficients_, argument - argument_min());
}

}  // namespace internal_polynomial
}  // namespace numerics
}  // namespace principia
