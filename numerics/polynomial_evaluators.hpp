#pragma once

#include "base/macros.hpp"
#include "numerics/polynomial.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace _polynomial_evaluators {
namespace internal {

using namespace principia::quantities::_named_quantities;

template<typename Value, typename Argument, int degree, bool allow_fma>
struct EstrinEvaluator;
template<typename Value, typename Argument, int degree, bool allow_fma>
struct HornerEvaluator;
}  // namespace internal

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

namespace _polynomial_evaluators {
namespace internal {
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
      numerics::EstrinEvaluator>::Coefficients;

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
      numerics::HornerEvaluator>::Coefficients;

  FORCE_INLINE(static) Value Evaluate(Coefficients const& coefficients,
                                      Argument const& argument);
  FORCE_INLINE(static) Derivative<Value, Argument>
  EvaluateDerivative(Coefficients const& coefficients,
                     Argument const& argument);
};

}  // namespace internal
}  // namespace _polynomial_evaluators
}  // namespace _polynomial_evaluators
}  // namespace numerics
}  // namespace principia

namespace principia::numerics {
using namespace principia::numerics::_polynomial_evaluators;
}  // namespace principia::numerics

#include "numerics/polynomial_evaluators_body.hpp"
