#pragma once

#include "numerics/polynomial.hpp"

#include "base/macros.hpp"

namespace principia {
namespace numerics {
namespace internal_polynomial {

using quantities::Exponentiation;
using quantities::Square;

namespace {

// Greatest power of 2 less than or equal to n.  8 -> 8, 7 -> 4.
constexpr int flp2(int const n) {
  return n == 0 ? 0 : n == 1 ? 1 : flp2(n >> 1) << 1;
}

// Ceiling log2 of n.  8 -> 3, 7 -> 2.
constexpr int log2(int const n) {
  return n == 1 ? 0 : flp2(n >> 1) + 1;
}

// Generator for repeated squaring:
//   SquareGenerator<Length, 0>::Type is Length
//   SquareGenerator<Length, 1>::Type is Exponentiation<Length, 2>
//   SquareGenerator<Length, 2>::Type is Exponentiation<Length, 4>
//   SquareGenerator<Length, n>::Type is Exponentiation<Length, 2^n>
// etc.
template<typename Argument, int n>
struct SquareGenerator {
  using Type = Square<typename SquareGenerator<Argument, n - 1>::Type>;

  Type Evaluate(Argument const& argument) {
    auto const argument_n_minus_1 =
        SquareGenerator<Argument, n - 1>::Evaluate(argument);
    return argument_n_minus_1 * argument_n_minus_1;
  }
};

template<typename Argument>
struct SquareGenerator<Argument, 0> {
  using Type = Argument;

  Type Evaluate(Argument const& argument) {
    return argument;
  }
};

template<typename Argument, typename>
struct SquaresGenerator;
template<typename Argument, int... orders>
struct SquaresGenerator<Argument, std::integer_sequence<int, orders...>> {
  using Type = std::tuple<typename SquareGenerator<Argument, orders>::Type...>;

  Type Evaluate(Argument const& argument) {
    return std::make_tuple<SquareGenerator<Argument, orders>::
               Evaluate(argument)...>;
  }
};

template<typename Value,
         typename Argument,
         int degree,
         int low = 0,
         int high = degree>
struct EstrinEvaluator {
  using ArgumentSquares = typename SquaresGenerator<
      Argument,
      std::make_integer_sequence<int, log2(degree)>>::Type;
  using Coefficients =
      typename PolynomialInMonomialBasis<Value, Argument, degree>::Coefficients;

  static NthDerivative<Value, Argument, low>
  Evaluate(Coefficients const& coefficients,
           Argument const& argument,
           ArgumentSquares& argument_squares);
};

template<typename Value, typename Argument, int degree, int low>
struct EstrinEvaluator<Value, Argument, degree, low, low> {
  using ArgumentSquares = typename SquaresGenerator<
      Argument,
      std::make_integer_sequence<int, log2(degree)>>::Type;
  using Coefficients =
      typename PolynomialInMonomialBasis<Value, Argument, degree>::Coefficients;

  static NthDerivative<Value, Argument, low> Evaluate(
      Coefficients const& coefficients,
      Argument const& argument,
      ArgumentSquares& argument_squares);
};

template<typename Value, typename Argument, int degree, int low, int high>
NthDerivative<Value, Argument, low>
EstrinEvaluator<Value, Argument, degree, low, high>::Evaluate(
    Coefficients const& coefficients,
    Argument const& argument,
    ArgumentSquares& argument_squares) {
  constexpr int n = flp2(high - low);
  return EstrinEvaluator<Value, Argument, degree, low, low + n - 1>::Evaluate(
             coefficients, argument, argument_squares) +
         std::get<n>(argument_squares) *
             EstrinEvaluator<Value, Argument, degree, low + n, high>::Evaluate(
                 coefficients, argument, argument_squares);
}

template<typename Value, typename Argument, int degree, int low>
NthDerivative<Value, Argument, low>
EstrinEvaluator<Value, Argument, degree, low, low>::Evaluate(
    Coefficients const& coefficients,
    Argument const& argument,
    ArgumentSquares& argument_squares) {
  return std::get<low>(coefficients);
}

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
           Argument const& argument);
  static FORCE_INLINE NthDerivative<Value, Argument, k>
  EvaluateDerivative(Coefficients const& coefficients,
                     Argument const& argument);
};

template<typename Value, typename Argument, int degree>
struct HornerEvaluator<Value, Argument, degree, degree> {
  using Coefficients =
      typename PolynomialInMonomialBasis<Value, Argument, degree>::Coefficients;

  static FORCE_INLINE NthDerivative<Value, Argument, degree>
  Evaluate(Coefficients const& coefficients,
           Argument const& argument);
  static FORCE_INLINE NthDerivative<Value, Argument, degree>
  EvaluateDerivative(Coefficients const& coefficients,
                     Argument const& argument);
};

template<typename Value, typename Argument, int degree, int k>
NthDerivative<Value, Argument, k>
HornerEvaluator<Value, Argument, degree, k>::Evaluate(
    Coefficients const& coefficients,
    Argument const& argument) {
  return std::get<k>(coefficients) +
         argument *
         HornerEvaluator<Value, Argument, degree, k + 1>::Evaluate(
             coefficients, argument);
}

template<typename Value, typename Argument, int degree, int k>
NthDerivative<Value, Argument, k>
HornerEvaluator<Value, Argument, degree, k>::EvaluateDerivative(
    Coefficients const& coefficients,
    Argument const& argument) {
  return std::get<k>(coefficients) * k +
         argument *
         HornerEvaluator<Value, Argument, degree, k + 1>::EvaluateDerivative(
             coefficients, argument);
}

template<typename Value, typename Argument, int degree>
NthDerivative<Value, Argument, degree>
HornerEvaluator<Value, Argument, degree, degree>::Evaluate(
    Coefficients const& coefficients,
    Argument const& argument) {
  return std::get<degree>(coefficients);
}

template<typename Value, typename Argument, int degree>
NthDerivative<Value, Argument, degree>
HornerEvaluator<Value, Argument, degree, degree>::EvaluateDerivative(
    Coefficients const& coefficients,
    Argument const& argument) {
  return std::get<degree>(coefficients) * degree;
}

}  // namespace

template<typename Value, typename Argument, int degree>
PolynomialInMonomialBasis<Value, Argument, degree>::PolynomialInMonomialBasis(
    Coefficients const& coefficients)
    : coefficients_(coefficients) {}

template<typename Value, typename Argument, int degree>
Value PolynomialInMonomialBasis<Value, Argument, degree>::Evaluate(
    Argument const& argument) const {
  return HornerEvaluator<Value, Argument, degree, 0>::Evaluate(
      coefficients_, argument);
}

template<typename Value, typename Argument, int degree>
Derivative<Value, Argument>
PolynomialInMonomialBasis<Value, Argument, degree>::EvaluateDerivative(
    Argument const& argument) const {
  return HornerEvaluator<Value, Argument, degree, 1>::EvaluateDerivative(
      coefficients_, argument);
}

}  // namespace internal_polynomial
}  // namespace numerics
}  // namespace principia
