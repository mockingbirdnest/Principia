#pragma once

#include "quantities/named_quantities.hpp"
#include "quantities/tuples.hpp"

namespace principia {
namespace numerics {
namespace _polynomial_evaluators {
namespace internal {

using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_tuples;

template<typename Value, typename Argument, int degree>
struct Evaluator;

template<typename Value, typename Argument, int degree, bool allow_fma>
struct EstrinEvaluator;
template<typename Value, typename Argument, int degree, bool allow_fma>
struct HornerEvaluator;

template<template<typename, typename, int> typename Evaluator_>
struct with_evaluator_t {};

template<template<typename, typename, int> typename Evaluator_>
static constexpr with_evaluator_t<Evaluator_> with_evaluator;

}  // namespace internal

template<typename Value, typename Argument, int degree>
using EstrinEvaluator = internal::
    EstrinEvaluator<Value, Argument, degree, /*allow_fma=*/true>;
template<typename Value, typename Argument, int degree>
using EstrinEvaluatorWithoutFMA = internal::
    EstrinEvaluator<Value, Argument, degree, /*allow_fma=*/false>;
template<typename Value, typename Argument, int degree>
using HornerEvaluator = internal::
    HornerEvaluator<Value, Argument, degree, /*allow_fma=*/true>;
template<typename Value, typename Argument, int degree>
using HornerEvaluatorWithoutFMA = internal::
    HornerEvaluator<Value, Argument, degree, /*allow_fma=*/false>;

namespace internal {

//TODO(phl)Cleanup
template<typename Value, typename Argument, int degree>
struct Evaluator {
  // This definition is replicated from |PolynomialInMonomialBasis| to avoid
  // circular dependencies.
  using Coefficients = Derivatives<Value, Argument, degree + 1>;

  virtual Value Evaluate(Coefficients const& coefficients,
                         Argument const& argument) const = 0;
  virtual Derivative<Value, Argument> EvaluateDerivative(
      Coefficients const& coefficients,
      Argument const& argument) const = 0;
};


// We use FORCE_INLINE because we have to write this recursively, but we really
// want linear code.

template<typename Value, typename Argument, int degree, bool allow_fma>
struct EstrinEvaluator : public Evaluator<Value, Argument, degree> {
  using typename Evaluator<Value, Argument, degree>::Coefficients;

  static Evaluator<Value, Argument, degree> const* Singleton();

  FORCE_INLINE(inline) Value Evaluate(Coefficients const& coefficients,
                                      Argument const& argument) const override;
  FORCE_INLINE(inline) Derivative<Value, Argument>
  EvaluateDerivative(Coefficients const& coefficients,
                     Argument const& argument) const override;
};

template<typename Value, typename Argument, int degree, bool allow_fma>
struct HornerEvaluator : public Evaluator<Value, Argument, degree> {
  using typename Evaluator<Value, Argument, degree>::Coefficients;

  static Evaluator<Value, Argument, degree> const* Singleton();

  FORCE_INLINE(inline) Value Evaluate(Coefficients const& coefficients,
                                      Argument const& argument) const override;
  FORCE_INLINE(inline) Derivative<Value, Argument>
  EvaluateDerivative(Coefficients const& coefficients,
                     Argument const& argument) const override;
};

}  // namespace internal
}  // namespace _polynomial_evaluators
}  // namespace numerics
}  // namespace principia

#include "numerics/polynomial_evaluators_body.hpp"
