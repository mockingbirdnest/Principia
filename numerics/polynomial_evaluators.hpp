#pragma once

#include "base/macros.hpp"
#include "numerics/polynomial.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_polynomial_evaluators {

using quantities::Derivative;
using quantities::Square;

template<typename Value, typename Argument, int degree, bool allow_fma>
struct EstrinEvaluator;
template<typename Value, typename Argument, int degree, bool allow_fma>
struct HornerEvaluator;

namespace exported {
template<typename Value, typename Argument, int degree>
using EstrinEvaluator = internal_polynomial_evaluators::
    EstrinEvaluator<Value, Argument, degree, /*allow_fma=*/true>;
template<typename Value, typename Argument, int degree>
using EstrinEvaluatorWithoutFMA = internal_polynomial_evaluators::
    EstrinEvaluator<Value, Argument, degree, /*allow_fma=*/false>;
template<typename Value, typename Argument, int degree>
using HornerEvaluator = internal_polynomial_evaluators::
    HornerEvaluator<Value, Argument, degree, /*allow_fma=*/true>;
template<typename Value, typename Argument, int degree>
using HornerEvaluatorWithoutFMA = internal_polynomial_evaluators::
    HornerEvaluator<Value, Argument, degree, /*allow_fma=*/false>;
}  // namespace exported

// We use FORCE_INLINE because we have to write this recursively, but we really
// want linear code.

template<typename Value, typename Argument, int degree, bool allow_fma>
struct EstrinEvaluator {
  // The fully qualified name below designates the template, not the current
  // instance.
  using Coefficients = typename PolynomialInMonomialBasis<
      Value,
      Argument,
      degree,
      exported::EstrinEvaluator>::Coefficients;

  FORCE_INLINE(static) Value Evaluate(Coefficients const& coefficients,
                                      Argument const& argument);
  FORCE_INLINE(static) Derivative<Value, Argument>
  EvaluateDerivative(Coefficients const& coefficients,
                     Argument const& argument);
};

template<typename Value, typename Argument, int degree, bool allow_fma>
struct HornerEvaluator {
  // The fully qualified name below designates the template, not the current
  // instance.
  using Coefficients = typename PolynomialInMonomialBasis<
      Value,
      Argument,
      degree,
      exported::HornerEvaluator>::Coefficients;

  FORCE_INLINE(static) Value Evaluate(Coefficients const& coefficients,
                                      Argument const& argument);
  FORCE_INLINE(static) Derivative<Value, Argument>
  EvaluateDerivative(Coefficients const& coefficients,
                     Argument const& argument);
};

}  // namespace internal_polynomial_evaluators

using internal_polynomial_evaluators::exported::EstrinEvaluator;
using internal_polynomial_evaluators::exported::EstrinEvaluatorWithoutFMA;
using internal_polynomial_evaluators::exported::HornerEvaluator;
using internal_polynomial_evaluators::exported::HornerEvaluatorWithoutFMA;

}  // namespace numerics
}  // namespace principia

#include "numerics/polynomial_evaluators_body.hpp"
