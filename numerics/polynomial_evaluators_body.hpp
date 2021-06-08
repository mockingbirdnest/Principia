#pragma once

#include "numerics/polynomial_evaluators.hpp"

#include <cstddef>
#include <tuple>

#include "base/bits.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace internal_polynomial_evaluators {

using base::FloorLog2;
using base::PowerOf2Le;
using quantities::FusedMultiplyAdd;

// Generator for repeated squaring:
//   SquareGenerator<Length, 0>::Type is Exponentiation<Length, 2>
//   SquareGenerator<Length, 1>::Type is Exponentiation<Length, 4>
//   SquareGenerator<Length, n>::Type is Exponentiation<Length, 2^(n + 1)>
// etc.
template<typename Argument, int n>
struct SquareGenerator {
  using Type = Square<typename SquareGenerator<Argument, n - 1>::Type>;
  static Type Evaluate(Argument const& argument);
};
template<typename Argument>
struct SquareGenerator<Argument, 0> {
  using Type = Square<Argument>;
  static Type Evaluate(Argument const& argument);
};

template<typename Argument, typename>
struct SquaresGenerator;
template<typename Argument, std::size_t... orders>
struct SquaresGenerator<Argument, std::index_sequence<orders...>> {
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
  return argument * argument;
}

template<typename Argument, std::size_t... orders>
auto SquaresGenerator<Argument, std::index_sequence<orders...>>::
    Evaluate(Argument const& argument) -> Type {
  return std::make_tuple(
      SquareGenerator<Argument, orders>::Evaluate(argument)...);
}

// Internal helper for Estrin evaluation.  |degree| is the degree of the overall
// polynomial, |low| and |subdegree| defines the subpolynomial that we currently
// evaluate, i.e., the one with a constant term coefficient
// |std::get<low>(coefficients)| and degree |subdegree|.
template<typename Value, typename Argument,
         int degree, bool fma, int low, int subdegree>
struct InternalEstrinEvaluator {
  using ArgumentSquaresGenerator =
      SquaresGenerator<Argument, std::make_index_sequence<FloorLog2(degree)>>;
  using ArgumentSquares = typename ArgumentSquaresGenerator::Type;
  using Coefficients =
      typename PolynomialInMonomialBasis<
          Value, Argument, degree, numerics::EstrinEvaluator>::Coefficients;

  FORCE_INLINE(static) Derivative<Value, Argument, low> Evaluate(
      Coefficients const& coefficients,
      Argument const& argument,
      ArgumentSquares const& argument_squares);
  FORCE_INLINE(static) Derivative<Value, Argument, low> EvaluateDerivative(
      Coefficients const& coefficients,
      Argument const& argument,
      ArgumentSquares const& argument_squares);
};

template<typename Value, typename Argument, int degree, bool fma, int low>
struct InternalEstrinEvaluator<Value, Argument, degree, fma, low, 1> {
  using ArgumentSquaresGenerator =
      SquaresGenerator<Argument, std::make_index_sequence<FloorLog2(degree)>>;
  using ArgumentSquares = typename ArgumentSquaresGenerator::Type;
  using Coefficients =
      typename PolynomialInMonomialBasis<
          Value, Argument, degree, numerics::EstrinEvaluator>::Coefficients;

  FORCE_INLINE(static) Derivative<Value, Argument, low> Evaluate(
      Coefficients const& coefficients,
      Argument const& argument,
      ArgumentSquares const& argument_squares);
  FORCE_INLINE(static) Derivative<Value, Argument, low> EvaluateDerivative(
      Coefficients const& coefficients,
      Argument const& argument,
      ArgumentSquares const& argument_squares);
};

template<typename Value, typename Argument, int degree, bool fma, int low>
struct InternalEstrinEvaluator<Value, Argument, degree, fma, low, 0> {
  using ArgumentSquaresGenerator =
      SquaresGenerator<Argument, std::make_index_sequence<FloorLog2(degree)>>;
  using ArgumentSquares = typename ArgumentSquaresGenerator::Type;
  using Coefficients =
      typename PolynomialInMonomialBasis<
          Value, Argument, degree, numerics::EstrinEvaluator>::Coefficients;

  FORCE_INLINE(static) Derivative<Value, Argument, low> Evaluate(
      Coefficients const& coefficients,
      Argument const& argument,
      ArgumentSquares const& argument_squares);
  FORCE_INLINE(static) Derivative<Value, Argument, low> EvaluateDerivative(
      Coefficients const& coefficients,
      Argument const& argument,
      ArgumentSquares const& argument_squares);
};

template<typename Value, typename Argument,
         int degree, bool fma, int low, int subdegree>
FORCE_INLINE(inline) Derivative<Value, Argument, low>
InternalEstrinEvaluator<Value, Argument, degree, fma, low, subdegree>::Evaluate(
    Coefficients const& coefficients,
    Argument const& argument,
    ArgumentSquares const& argument_squares) {
  static_assert(subdegree >= 2,
                "Unexpected subdegree in InternalEstrinEvaluator::Evaluate");
  // |n| is used to select |argument^(2^(n + 1))| = |argument^m|.
  constexpr int n = FloorLog2(subdegree) - 1;
  // |m| is |2^(n + 1)|.
  constexpr int m = PowerOf2Le(subdegree);
  auto const& xᵐ = std::get<n>(argument_squares);
  auto const a =
      InternalEstrinEvaluator<Value, Argument,
                              degree, fma, low + m, subdegree - m>::
          Evaluate(coefficients, argument, argument_squares);
  auto const b =
      InternalEstrinEvaluator<Value, Argument, degree, fma, low, m - 1>::
          Evaluate(coefficients, argument, argument_squares);
  if constexpr (fma) {
    return FusedMultiplyAdd(a, xᵐ, b);
  } else {
    return  a * xᵐ + b;
  }
}

template<typename Value, typename Argument,
         int degree, bool fma, int low, int subdegree>
FORCE_INLINE(inline) Derivative<Value, Argument, low>
InternalEstrinEvaluator<Value, Argument, degree, fma, low, subdegree>::
EvaluateDerivative(Coefficients const& coefficients,
                   Argument const& argument,
                   ArgumentSquares const& argument_squares) {
  static_assert(subdegree >= 2,
                "Unexpected subdegree in InternalEstrinEvaluator::"
                "EvaluateDerivative");
  // |n| is used to select |argument^(2^(n + 1))| = |argument^m|.
  constexpr int n = FloorLog2(subdegree) - 1;
  // |m| is |2^(n + 1)|.
  constexpr int m = PowerOf2Le(subdegree);
  auto const& xᵐ = std::get<n>(argument_squares);
  auto const a =
      InternalEstrinEvaluator<Value, Argument,
                              degree, fma, low + m, subdegree - m>::
          EvaluateDerivative(coefficients, argument, argument_squares);
  auto const b =
      InternalEstrinEvaluator<Value, Argument, degree, fma, low, m - 1>::
          EvaluateDerivative(coefficients, argument, argument_squares);
  if constexpr (fma) {
    return FusedMultiplyAdd(a, xᵐ, b);
  } else {
    return a * xᵐ + b;
  }
}

template<typename Value, typename Argument, int degree, bool fma, int low>
FORCE_INLINE(inline) Derivative<Value, Argument, low>
InternalEstrinEvaluator<Value, Argument, degree, fma, low, 1>::Evaluate(
    Coefficients const& coefficients,
    Argument const& argument,
    ArgumentSquares const& argument_squares) {
  auto const& x = argument;
  auto const& a = std::get<low + 1>(coefficients);
  auto const& b = std::get<low>(coefficients);
  if constexpr (fma) {
    return FusedMultiplyAdd(a, x, b);
  } else {
    return a * x + b;
  }
}

template<typename Value, typename Argument, int degree, bool fma, int low>
FORCE_INLINE(inline) Derivative<Value, Argument, low>
InternalEstrinEvaluator<Value, Argument, degree, fma, low, 0>::Evaluate(
    Coefficients const& coefficients,
    Argument const& argument,
    ArgumentSquares const& argument_squares) {
  return std::get<low>(coefficients);
}

template<typename Value, typename Argument, int degree, bool fma, int low>
FORCE_INLINE(inline) Derivative<Value, Argument, low>
InternalEstrinEvaluator<Value, Argument, degree, fma, low, 1>::
    EvaluateDerivative(Coefficients const& coefficients,
                       Argument const& argument,
                       ArgumentSquares const& argument_squares) {
  auto const& x = argument;
  auto const& a = (low + 1) * std::get<low + 1>(coefficients);
  auto const& b = low * std::get<low>(coefficients);
  if constexpr (fma) {
    return FusedMultiplyAdd(a, x, b);
  } else {
    return a * x + b;
  }
}

template<typename Value, typename Argument, int degree, bool fma, int low>
FORCE_INLINE(inline) Derivative<Value, Argument, low>
InternalEstrinEvaluator<Value, Argument, degree, fma, low, 0>::
    EvaluateDerivative(Coefficients const& coefficients,
                       Argument const& argument,
                       ArgumentSquares const& argument_squares) {
  return low * std::get<low>(coefficients);
}

template<typename Value, typename Argument, int degree, bool allow_fma>
FORCE_INLINE(inline) Value
EstrinEvaluator<Value, Argument, degree, allow_fma>::Evaluate(
    Coefficients const& coefficients,
    Argument const& argument) {
  if (allow_fma && UseHardwareFMA) {
    using InternalEvaluator = InternalEstrinEvaluator<Value,
                                                      Argument,
                                                      degree,
                                                      /*fma=*/true,
                                                      /*low=*/0,
                                                      /*subdegree=*/degree>;
    return InternalEvaluator::Evaluate(
        coefficients,
        argument,
        InternalEvaluator::ArgumentSquaresGenerator::Evaluate(argument));
  } else {
    using InternalEvaluator = InternalEstrinEvaluator<Value,
                                                      Argument,
                                                      degree,
                                                      /*fma=*/false,
                                                      /*low=*/0,
                                                      /*subdegree=*/degree>;
    return InternalEvaluator::Evaluate(
        coefficients,
        argument,
        InternalEvaluator::ArgumentSquaresGenerator::Evaluate(argument));
  }
}

template<typename Value, typename Argument, int degree, bool allow_fma>
FORCE_INLINE(inline) Derivative<Value, Argument>
EstrinEvaluator<Value, Argument, degree, allow_fma>::EvaluateDerivative(
    Coefficients const& coefficients,
    Argument const& argument) {
  if constexpr (degree == 0) {
    return Derivative<Value, Argument>{};
  } else if (allow_fma && UseHardwareFMA) {
    using InternalEvaluator =
        InternalEstrinEvaluator<Value,
                                Argument,
                                degree,
                                true,
                                /*low=*/1,
                                /*subdegree=*/degree - 1>;
    return InternalEvaluator::EvaluateDerivative(
        coefficients,
        argument,
        InternalEvaluator::ArgumentSquaresGenerator::Evaluate(argument));
  } else {
    using InternalEvaluator =
        InternalEstrinEvaluator<Value,
                                Argument,
                                degree,
                                false,
                                /*low=*/1,
                                /*subdegree=*/degree - 1>;
    return InternalEvaluator::EvaluateDerivative(
        coefficients,
        argument,
        InternalEvaluator::ArgumentSquaresGenerator::Evaluate(argument));
  }
}

// Internal helper for Horner evaluation.  |degree| is the degree of the overall
// polynomial, |low| defines the subpolynomial that we currently evaluate, i.e.,
// the one with a constant term coefficient |std::get<low>(coefficients)|.
template<typename Value, typename Argument, int degree, bool fma, int low>
struct InternalHornerEvaluator {
  using Coefficients =
      typename PolynomialInMonomialBasis<
          Value, Argument, degree, numerics::HornerEvaluator>::Coefficients;

  FORCE_INLINE(static) Derivative<Value, Argument, low>
  Evaluate(Coefficients const& coefficients,
           Argument const& argument);
  FORCE_INLINE(static) Derivative<Value, Argument, low>
  EvaluateDerivative(Coefficients const& coefficients,
                     Argument const& argument);
};

template<typename Value, typename Argument, int degree, bool fma>
struct InternalHornerEvaluator<Value, Argument, degree, fma, degree> {
  using Coefficients =
      typename PolynomialInMonomialBasis<
          Value, Argument, degree, numerics::HornerEvaluator>::Coefficients;

  FORCE_INLINE(static) Derivative<Value, Argument, degree>
  Evaluate(Coefficients const& coefficients,
           Argument const& argument);
  FORCE_INLINE(static) Derivative<Value, Argument, degree>
  EvaluateDerivative(Coefficients const& coefficients,
                     Argument const& argument);
};

template<typename Value, typename Argument, int degree, bool fma, int low>
FORCE_INLINE(inline) Derivative<Value, Argument, low>
InternalHornerEvaluator<Value, Argument, degree, fma, low>::Evaluate(
    Coefficients const& coefficients,
    Argument const& argument) {
  auto const& x = argument;
  auto const a =
      InternalHornerEvaluator<Value, Argument, degree, fma, low + 1>::Evaluate(
          coefficients, argument);
  auto const& b = std::get<low>(coefficients);
  if constexpr (fma) {
    return FusedMultiplyAdd(a, x, b);
  } else {
    return a * x + b;
  }
}

template<typename Value, typename Argument, int degree, bool fma, int low>
FORCE_INLINE(inline) Derivative<Value, Argument, low>
InternalHornerEvaluator<Value, Argument, degree, fma, low>::EvaluateDerivative(
    Coefficients const& coefficients,
    Argument const& argument) {
  auto const& x = argument;
  auto const a =
      InternalHornerEvaluator<Value, Argument, degree, fma, low + 1>::
          EvaluateDerivative(coefficients, argument);
  auto const b = std::get<low>(coefficients) * low;
  if constexpr (fma) {
    return FusedMultiplyAdd(a, x, b);
  } else {
    return a * x + b;
  }
}

template<typename Value, typename Argument, int degree, bool fma>
FORCE_INLINE(inline) Derivative<Value, Argument, degree>
InternalHornerEvaluator<Value, Argument, degree, fma, degree>::Evaluate(
    Coefficients const& coefficients,
    Argument const& argument) {
  return std::get<degree>(coefficients);
}

template<typename Value, typename Argument, int degree, bool fma>
FORCE_INLINE(inline) Derivative<Value, Argument, degree>
InternalHornerEvaluator<Value, Argument, degree, fma, degree>::
    EvaluateDerivative(Coefficients const& coefficients,
                        Argument const& argument) {
  return std::get<degree>(coefficients) * degree;
}

template<typename Value, typename Argument, int degree, bool allow_fma>
FORCE_INLINE(inline) Value
HornerEvaluator<Value, Argument, degree, allow_fma>::Evaluate(
    Coefficients const& coefficients,
    Argument const& argument) {
  if (allow_fma && UseHardwareFMA) {
    return InternalHornerEvaluator<
        Value, Argument, degree, /*fma=*/true, /*low=*/0>::Evaluate(
            coefficients, argument);
  } else {
    return InternalHornerEvaluator<
        Value, Argument, degree, /*fma=*/false, /*low=*/0>::Evaluate(
            coefficients, argument);
  }
}

template<typename Value, typename Argument, int degree, bool allow_fma>
FORCE_INLINE(inline) Derivative<Value, Argument>
HornerEvaluator<Value, Argument, degree, allow_fma>::EvaluateDerivative(
    Coefficients const& coefficients,
    Argument const& argument) {
  if constexpr (degree == 0) {
    return Derivative<Value, Argument>{};
  } else if (allow_fma && UseHardwareFMA) {
    return InternalHornerEvaluator<
        Value, Argument, degree,
        /*fma=*/true, /*low=*/1>::EvaluateDerivative(
            coefficients, argument);
  } else {
    return InternalHornerEvaluator<
        Value, Argument, degree,
        /*fma=*/false, /*low=*/1>::EvaluateDerivative(
            coefficients, argument);
  }
}

}  // namespace internal_polynomial_evaluators
}  // namespace numerics
}  // namespace principia
