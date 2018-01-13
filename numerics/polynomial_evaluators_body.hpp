#pragma once

#include "numerics/polynomial_evaluators.hpp"

#include <tuple>

namespace principia {
namespace numerics {
namespace internal_polynomial_evaluators {

namespace {

// Greatest power of 2 less than or equal to n.  8 -> 8, 7 -> 4.
constexpr int FloorOfPowerOf2(int const n) {
  return n == 0 ? 0 : n == 1 ? 1 : FloorOfPowerOf2(n >> 1) << 1;
}

// Ceiling log2 of n.  8 -> 3, 7 -> 2.
constexpr int CeilingLog2(int const n) {
  return n == 1 ? 0 : CeilingLog2(n >> 1) + 1;
}

}  // namespace

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

template<typename Argument, int n>
auto SquareGenerator<Argument, n>::Evaluate(Argument const& argument) -> Type {
  auto const argument_n_minus_1 =
      SquareGenerator<Argument, n - 1>::Evaluate(argument);
  return argument_n_minus_1 * argument_n_minus_1;
}

template<typename Argument>
auto SquareGenerator<Argument, 0>::Evaluate(Argument const& argument) -> Type {
  return argument;
}

template<typename Argument, int... orders>
auto SquaresGenerator<Argument, std::integer_sequence<int, orders...>>::
    Evaluate(Argument const& argument) -> Type {
  return std::make_tuple(
      SquareGenerator<Argument, orders>::Evaluate(argument)...);
}

// Internal helper for Estrin evaluation.  |degree| is the degree of the overall
// polynomial, |low| and |high| defines the subpolynomial that we currently
// evaluate, i.e., the one with a constant term coefficient
// |std::get<low>(coefficients)| and a highest degree coefficient
// |std::get<high>(coefficients)|.
template<typename Value, typename Argument, int degree, int low, int high>
struct InternalEstrinEvaluator {
  using ArgumentSquaresGenerator =
      SquaresGenerator<Argument,
                       std::make_integer_sequence<int,
                                                  CeilingLog2(degree) + 1>>;
  using ArgumentSquares = typename ArgumentSquaresGenerator::Type;
  using Coefficients =
      typename PolynomialInMonomialBasis<Value, Argument, degree,
                                         EstrinEvaluator>::Coefficients;

  static FORCE_INLINE NthDerivative<Value, Argument, low> Evaluate(
      Coefficients const& coefficients,
      Argument const& argument,
      ArgumentSquares const& argument_squares);
};

template<typename Value, typename Argument, int degree, int low>
struct InternalEstrinEvaluator<Value, Argument, degree, low, low> {
  using ArgumentSquaresGenerator =
      SquaresGenerator<Argument,
                       std::make_integer_sequence<int,
                                                  CeilingLog2(degree) + 1>>;
  using ArgumentSquares = typename ArgumentSquaresGenerator::Type;
  using Coefficients =
      typename PolynomialInMonomialBasis<Value, Argument, degree,
                                         EstrinEvaluator>::Coefficients;

  static FORCE_INLINE NthDerivative<Value, Argument, low> Evaluate(
      Coefficients const& coefficients,
      Argument const& argument,
      ArgumentSquares const& argument_squares);
};

template<typename Value, typename Argument, int degree, int low, int high>
NthDerivative<Value, Argument, low>
InternalEstrinEvaluator<Value, Argument, degree, low, high>::Evaluate(
    Coefficients const& coefficients,
    Argument const& argument,
    ArgumentSquares const& argument_squares) {
  constexpr int n = CeilingLog2(high - low);
  constexpr int m = FloorOfPowerOf2(high - low);  // m = 2^n
  return InternalEstrinEvaluator<Value, Argument, degree, low, low + m - 1>::
             Evaluate(coefficients, argument, argument_squares) +
         std::get<n>(argument_squares) *
             InternalEstrinEvaluator<Value, Argument, degree, low + m, high>::
                 Evaluate(coefficients, argument, argument_squares);
}

template<typename Value, typename Argument, int degree, int low>
NthDerivative<Value, Argument, low>
InternalEstrinEvaluator<Value, Argument, degree, low, low>::Evaluate(
    Coefficients const& coefficients,
    Argument const& argument,
    ArgumentSquares const& argument_squares) {
  return std::get<low>(coefficients);
}

template<typename Value, typename Argument, int degree>
Value EstrinEvaluator<Value, Argument, degree>::Evaluate(
    Coefficients const& coefficients,
    Argument const& argument) {
  using InternalEvaluator = InternalEstrinEvaluator<Value,
                                                    Argument,
                                                    degree,
                                                    /*low=*/0,
                                                    /*high=*/degree>;
  return InternalEvaluator::Evaluate(
      coefficients,
      argument,
      InternalEvaluator::ArgumentSquaresGenerator::Evaluate(argument));
}

// Internal helper for Horner evaluation.  |degree| is the degree of the overall
// polynomial, |low| defines the subpolynomial that we currently evaluate, i.e.,
// the one with a constant term coefficient |std::get<low>(coefficients)|.
template<typename Value, typename Argument, int degree, int low>
struct InternalHornerEvaluator {
  using Coefficients =
      typename PolynomialInMonomialBasis<Value, Argument, degree,
                                         HornerEvaluator>::Coefficients;

  static FORCE_INLINE NthDerivative<Value, Argument, low>
  Evaluate(Coefficients const& coefficients,
           Argument const& argument);
  static FORCE_INLINE NthDerivative<Value, Argument, low>
  EvaluateDerivative(Coefficients const& coefficients,
                     Argument const& argument);
};

template<typename Value, typename Argument, int degree>
struct InternalHornerEvaluator<Value, Argument, degree, degree> {
  using Coefficients =
      typename PolynomialInMonomialBasis<Value, Argument, degree,
                                         HornerEvaluator>::Coefficients;

  static FORCE_INLINE NthDerivative<Value, Argument, degree>
  Evaluate(Coefficients const& coefficients,
           Argument const& argument);
  static FORCE_INLINE NthDerivative<Value, Argument, degree>
  EvaluateDerivative(Coefficients const& coefficients,
                     Argument const& argument);
};

template<typename Value, typename Argument, int degree, int low>
NthDerivative<Value, Argument, low>
InternalHornerEvaluator<Value, Argument, degree, low>::Evaluate(
    Coefficients const& coefficients,
    Argument const& argument) {
  return std::get<low>(coefficients) +
         argument *
         InternalHornerEvaluator<Value, Argument, degree, low + 1>::
             Evaluate(coefficients, argument);
}

template<typename Value, typename Argument, int degree, int low>
NthDerivative<Value, Argument, low>
InternalHornerEvaluator<Value, Argument, degree, low>::EvaluateDerivative(
    Coefficients const& coefficients,
    Argument const& argument) {
  return std::get<low>(coefficients) * low +
         argument *
         InternalHornerEvaluator<Value, Argument, degree, low + 1>::
             EvaluateDerivative(coefficients, argument);
}

template<typename Value, typename Argument, int degree>
NthDerivative<Value, Argument, degree>
InternalHornerEvaluator<Value, Argument, degree, degree>::Evaluate(
    Coefficients const& coefficients,
    Argument const& argument) {
  return std::get<degree>(coefficients);
}

template<typename Value, typename Argument, int degree>
NthDerivative<Value, Argument, degree>
InternalHornerEvaluator<Value, Argument, degree, degree>::EvaluateDerivative(
    Coefficients const& coefficients,
    Argument const& argument) {
  return std::get<degree>(coefficients) * degree;
}

template<typename Value, typename Argument, int degree>
Value HornerEvaluator<Value, Argument, degree>::Evaluate(
    Coefficients const& coefficients,
    Argument const& argument) {
  return InternalHornerEvaluator<Value, Argument, degree, /*low=*/0>::Evaluate(
      coefficients, argument);
}

template<typename Value, typename Argument, int degree>
Derivative<Value, Argument>
HornerEvaluator<Value, Argument, degree>::EvaluateDerivative(
    Coefficients const& coefficients,
    Argument const& argument) {
  return InternalHornerEvaluator<Value, Argument, degree, /*low=*/1>::
      EvaluateDerivative(coefficients, argument);
}

}  // namespace internal_polynomial_evaluators
}  // namespace numerics
}  // namespace principia
