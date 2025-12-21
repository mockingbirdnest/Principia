#pragma once

#include "base/not_null.hpp"
#include "geometry/concepts.hpp"
#include "numerics/fma.hpp"
#include "quantities/arithmetic.hpp"
#include "quantities/tuples.hpp"
#include "serialization/numerics.pb.h"

namespace principia {
namespace numerics {

class PolynomialEvaluatorTest;

namespace _polynomial_evaluators {
namespace internal {

using namespace principia::base::_not_null;
using namespace principia::geometry::_concepts;
using namespace principia::numerics::_fma;
using namespace principia::quantities::_arithmetic;
using namespace principia::quantities::_tuples;

// This definition is replicated from `PolynomialInMonomialBasis` to avoid
// circular dependencies.  Don't export it from here.
template<typename Value, typename Argument, int degree>
using Coefficients = Derivatives<Value, Argument, degree + 1>;

template<typename Value, typename Argument, int degree,
         FMAPolicy fma_policy, FMAPresence fma_presence>
class EstrinEvaluator {
  using Coefficients = internal::Coefficients<Value, Argument, degree>;

 public:
  FORCE_INLINE(static) Value
  Evaluate(Coefficients const& coefficients,
           Argument const& argument);
  FORCE_INLINE(static) Derivative<Value, Argument>
  EvaluateDerivative(Coefficients const& coefficients,
                     Argument const& argument);
  FORCE_INLINE(static) void
  EvaluateWithDerivative(Coefficients const& coefficients,
                         Argument const& argument,
                         Value& value,
                         Derivative<Value, Argument>& derivative);

  static void WriteToMessage(
      not_null<serialization::PolynomialInMonomialBasis::Evaluator*> message);

 private:
  EstrinEvaluator() = default;

  friend class PolynomialEvaluatorTest;
};

template<typename Value, typename Argument, int degree,
         FMAPolicy fma_policy, FMAPresence fma_presence>
class HornerEvaluator {
  using Coefficients = internal::Coefficients<Value, Argument, degree>;

 public:
  FORCE_INLINE(static) Value
  Evaluate(Coefficients const& coefficients,
           Argument const& argument);
  FORCE_INLINE(static) Derivative<Value, Argument>
  EvaluateDerivative(Coefficients const& coefficients,
                     Argument const& argument);
  FORCE_INLINE(static) void
  EvaluateWithDerivative(Coefficients const& coefficients,
                         Argument const& argument,
                         Value& value,
                         Derivative<Value, Argument>& derivative);

  static void WriteToMessage(
      not_null<serialization::PolynomialInMonomialBasis::Evaluator*> message);

 private:
  HornerEvaluator() = default;

  friend class PolynomialEvaluatorTest;
};

}  // namespace internal

template<typename Value, typename Argument, int degree>
using Estrin = internal::EstrinEvaluator<Value, Argument, degree,
                                         internal::FMAPolicy::Auto,
                                         internal::FMAPresence::Unknown>;
template<typename Value, typename Argument, int degree>
using EstrinWithoutFMA =
    internal::EstrinEvaluator<Value, Argument, degree,
                              internal::FMAPolicy::Disallow,
                              internal::FMAPresence::Unknown>;
template<typename Value, typename Argument, int degree>
using Horner = internal::HornerEvaluator<Value, Argument, degree,
                                         internal::FMAPolicy::Auto,
                                         internal::FMAPresence::Unknown>;
template<typename Value, typename Argument, int degree>
using HornerWithoutFMA =
    internal::HornerEvaluator<Value, Argument, degree,
                              internal::FMAPolicy::Disallow,
                              internal::FMAPresence::Unknown>;

using internal::EstrinEvaluator;
using internal::HornerEvaluator;

}  // namespace _polynomial_evaluators
}  // namespace numerics
}  // namespace principia

#include "numerics/polynomial_evaluators_body.hpp"
