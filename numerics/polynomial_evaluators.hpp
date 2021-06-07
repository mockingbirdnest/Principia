#pragma once

#include "base/macros.hpp"
#include "numerics/polynomial.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_polynomial_evaluators {

using quantities::Derivative;
using quantities::Square;

// The structs below are only needed because the language won't let us have (1)
// a partial specialization of a function (2) an explicit specialization of a
// member function template in an unspecialized class template.  Sigh.
// We use FORCE_INLINE because we have to write this recursively, but we really
// want linear code.

template<typename Value, typename Argument, int degree, bool allow_fma = true>
struct EstrinEvaluator {
  // The fully qualified name below designates the template, not the current
  // instance.
  using Coefficients = typename PolynomialInMonomialBasis<
      Value,
      Argument,
      degree,
      internal_polynomial_evaluators::EstrinEvaluator>::Coefficients;

  template<bool fma>
  FORCE_INLINE(static) Value Evaluate(Coefficients const& coefficients,
                                      Argument const& argument);

  template<bool fma>
  FORCE_INLINE(static) Derivative<Value, Argument>
  EvaluateDerivative(Coefficients const& coefficients,
                     Argument const& argument);
};

template<typename Value, typename Argument, int degree>
using EstrinEvaluatorWithoutFMA =
    EstrinEvaluator<Value, Argument, degree, /*allow_fma=*/false>;

template<typename Value, typename Argument, int degree, bool allow_fma = true>
struct HornerEvaluator {
  // The fully qualified name below designates the template, not the current
  // instance.
  using Coefficients = typename PolynomialInMonomialBasis<
      Value,
      Argument,
      degree,
      internal_polynomial_evaluators::HornerEvaluator>::Coefficients;

  template<bool fma>
  FORCE_INLINE(static) Value Evaluate(Coefficients const& coefficients,
                                      Argument const& argument);

  template<bool fma>
  FORCE_INLINE(static) Derivative<Value, Argument>
  EvaluateDerivative(Coefficients const& coefficients,
                     Argument const& argument);
};

template<typename Value, typename Argument, int degree>
using HornerEvaluatorWithoutFMA =
    HornerEvaluator<Value, Argument, degree, /*allow_fma=*/false>;

}  // namespace internal_polynomial_evaluators

using internal_polynomial_evaluators::EstrinEvaluator;
using internal_polynomial_evaluators::EstrinEvaluatorWithoutFMA;
using internal_polynomial_evaluators::HornerEvaluator;
using internal_polynomial_evaluators::HornerEvaluatorWithoutFMA;

}  // namespace numerics
}  // namespace principia

#include "numerics/polynomial_evaluators_body.hpp"
