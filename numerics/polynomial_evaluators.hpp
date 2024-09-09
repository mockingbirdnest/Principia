#pragma once

#include "base/not_null.hpp"
#include "numerics/fma.hpp"
#include "quantities/concepts.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/tuples.hpp"
#include "serialization/numerics.pb.h"

namespace principia {
namespace numerics {
namespace _polynomial_evaluators {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::numerics::_fma;
using namespace principia::quantities::_concepts;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_tuples;

template<template<typename, typename, int> typename Evaluator_>
struct with_evaluator_t {};

template<template<typename, typename, int> typename Evaluator_>
static constexpr with_evaluator_t<Evaluator_> with_evaluator;

template<typename Value, typename Argument, int degree>
  requires additive_group<Argument>
struct Evaluator {
  // This definition is replicated from `PolynomialInMonomialBasis` to avoid
  // circular dependencies.
  using Coefficients = Derivatives<Value, Argument, degree + 1>;

  FORCE_INLINE(static)
  Value Evaluate(
      Coefficients const& coefficients,
      Argument const& argument,
      not_null<Evaluator const*> evaluator);
  FORCE_INLINE(static)
  Derivative<Value, Argument> EvaluateDerivative(
      Coefficients const& coefficients,
      Argument const& argument,
      not_null<Evaluator const*> evaluator);

  static void WriteToMessage(
      not_null<serialization::PolynomialInMonomialBasis::Evaluator*> message,
      not_null<Evaluator const*> evaluator);
  static not_null<Evaluator const*> ReadFromMessage(
      serialization::PolynomialInMonomialBasis::Evaluator const& message);
};

template<typename Value, typename Argument, int degree, FMAPolicy fma_policy>
class EstrinEvaluator;
template<typename Value, typename Argument, int degree, FMAPolicy fma_policy>
class HornerEvaluator;

}  // namespace internal

template<typename Value, typename Argument, int degree>
using Estrin = internal::
    EstrinEvaluator<Value, Argument, degree, internal::FMAPolicy::Auto>;
template<typename Value, typename Argument, int degree>
using EstrinWithoutFMA = internal::
    EstrinEvaluator<Value, Argument, degree, internal::FMAPolicy::Disallow>;
template<typename Value, typename Argument, int degree>
using Horner = internal::
    HornerEvaluator<Value, Argument, degree, internal::FMAPolicy::Auto>;
template<typename Value, typename Argument, int degree>
using HornerWithoutFMA = internal::
    HornerEvaluator<Value, Argument, degree, internal::FMAPolicy::Disallow>;

using internal::EstrinEvaluator;
using internal::Evaluator;
using internal::HornerEvaluator;
using internal::with_evaluator;
using internal::with_evaluator_t;

}  // namespace _polynomial_evaluators
}  // namespace numerics
}  // namespace principia

#include "numerics/polynomial_evaluators_body.hpp"
