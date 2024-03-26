#pragma once

#include "quantities/named_quantities.hpp"
#include "quantities/tuples.hpp"

namespace principia {
namespace numerics {
namespace _polynomial_evaluators {
namespace internal {

using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_tuples;

template<template<typename, typename, int> typename Evaluator_>
struct with_evaluator_t {};

template<template<typename, typename, int> typename Evaluator_>
static constexpr with_evaluator_t<Evaluator_> with_evaluator;

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

template<typename Value, typename Argument, int degree, bool allow_fma>
struct EstrinEvaluator;
template<typename Value, typename Argument, int degree, bool allow_fma>
struct HornerEvaluator;

}  // namespace internal

template<typename Value, typename Argument, int degree>
using Estrin = internal::
    EstrinEvaluator<Value, Argument, degree, /*allow_fma=*/true>;
template<typename Value, typename Argument, int degree>
using EstrinWithoutFMA = internal::
    EstrinEvaluator<Value, Argument, degree, /*allow_fma=*/false>;
template<typename Value, typename Argument, int degree>
using Horner = internal::
    HornerEvaluator<Value, Argument, degree, /*allow_fma=*/true>;
template<typename Value, typename Argument, int degree>
using HornerWithoutFMA = internal::
    HornerEvaluator<Value, Argument, degree, /*allow_fma=*/false>;

using internal::Evaluator;
using internal::with_evaluator;
using internal::with_evaluator_t;

}  // namespace _polynomial_evaluators
}  // namespace numerics
}  // namespace principia

#include "numerics/polynomial_evaluators_body.hpp"
