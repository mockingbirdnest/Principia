#pragma once

#include "base/macros.hpp"
#include "numerics/polynomial.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {
namespace internal_polynomial_evaluator {

using quantities::Derivative;
using quantities::Square;

// The structs below are only needed because the language won't let us have (1)
// a partial specialization of a function (2) an explicit specialization of a
// member function template in an unspecialized class template.  Sigh.
// We use FORCE_INLINE because we have to write this recursively, but we really
// want linear code.

// Generator for repeated squaring:
//   SquareGenerator<Length, 0>::Type is Length
//   SquareGenerator<Length, 1>::Type is Exponentiation<Length, 2>
//   SquareGenerator<Length, 2>::Type is Exponentiation<Length, 4>
//   SquareGenerator<Length, n>::Type is Exponentiation<Length, 2^n>
// etc.
template<typename Argument, int n>
struct SquareGenerator {
  using Type = Square<typename SquareGenerator<Argument, n - 1>::Type>;
  static Type Evaluate(Argument const& argument);
};
template<typename Argument>
struct SquareGenerator<Argument, 0> {
  using Type = Argument;
  static Type Evaluate(Argument const& argument);
};

template<typename Argument, typename>
struct SquaresGenerator;
template<typename Argument, int... orders>
struct SquaresGenerator<Argument, std::integer_sequence<int, orders...>> {
  using Type = std::tuple<typename SquareGenerator<Argument, orders>::Type...>;
  static Type Evaluate(Argument const& argument);
};

template<typename Value, typename Argument, int degree>
struct EstrinEvaluator {
  using Coefficients =
      typename PolynomialInMonomialBasis<Value, Argument, degree>::Coefficients;

  static FORCE_INLINE Value Evaluate(Coefficients const& coefficients,
                                     Argument const& argument);
};

template<typename Value, typename Argument, int degree>
struct HornerEvaluator {
  using Coefficients =
      typename PolynomialInMonomialBasis<Value, Argument, degree>::Coefficients;

  static FORCE_INLINE Value Evaluate(Coefficients const& coefficients,
                                     Argument const& argument);
  static FORCE_INLINE Derivative<Value, Argument>
  EvaluateDerivative(Coefficients const& coefficients,
                     Argument const& argument);
};

}  // namespace internal_polynomial_evaluator

using internal_polynomial_evaluator::EstrinEvaluator;
using internal_polynomial_evaluator::HornerEvaluator;

}  // namespace numerics
}  // namespace principia

#include "numerics/polynomial_evaluators_body.hpp"
