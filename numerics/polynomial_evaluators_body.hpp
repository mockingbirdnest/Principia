#pragma once

#include "numerics/polynomial_evaluators.hpp"

#include <cstddef>
#include <tuple>

#include "base/bits.hpp"
#include "numerics/fma.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace _polynomial_evaluators {
namespace internal {

using namespace principia::base::_bits;
using namespace principia::numerics::_fma;
using namespace principia::quantities::_elementary_functions;

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
      typename Evaluator<Value, Argument, degree>::Coefficients;

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
      typename Evaluator<Value, Argument, degree>::Coefficients;

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
      typename Evaluator<Value, Argument, degree>::Coefficients;

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
    using quantities::_elementary_functions::FusedMultiplyAdd;
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
    using quantities::_elementary_functions::FusedMultiplyAdd;
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
    using quantities::_elementary_functions::FusedMultiplyAdd;
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
    using quantities::_elementary_functions::FusedMultiplyAdd;
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

// Internal helper for Horner evaluation.  |degree| is the degree of the overall
// polynomial, |low| defines the subpolynomial that we currently evaluate, i.e.,
// the one with a constant term coefficient |std::get<low>(coefficients)|.
template<typename Value, typename Argument, int degree, bool fma, int low>
struct InternalHornerEvaluator {
  using Coefficients =
      typename Evaluator<Value, Argument, degree>::Coefficients;

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
      typename Evaluator<Value, Argument, degree>::Coefficients;

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
    using quantities::_elementary_functions::FusedMultiplyAdd;
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
    using quantities::_elementary_functions::FusedMultiplyAdd;
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


template<typename Value, typename Argument, int degree>
  requires additive_group<Argument>
Value Evaluator<Value, Argument, degree>::Evaluate(
    Coefficients const& coefficients,
    Argument const& argument,
    not_null<Evaluator const*> const evaluator) {
  // We cannot make this function virtual, as that prevents inlining and the
  // performance evaluation is severe (polynomial evaluation takes nanoseconds).
  // Instead we dispatch by hand.  We use the singleton as an identity for the
  // types since the constructors are hidden and there is only ever one object
  // of each evaluator type.  There are no calls in this code.
  if (evaluator == Estrin<Value, Argument, degree>::Singleton()) {
    return Estrin<Value, Argument, degree>::Evaluate(
        coefficients, argument);
  } else if (evaluator == Horner<Value, Argument, degree>::Singleton()) {
    return Horner<Value, Argument, degree>::Evaluate(
        coefficients, argument);
  } else if (evaluator ==
             EstrinWithoutFMA<Value, Argument, degree>::Singleton()) {
    return EstrinWithoutFMA<Value, Argument, degree>::Evaluate(
        coefficients, argument);
  } else {
    /*evaluator == HornerWithoutFMA<Value, Argument, degree>::Singleton())*/
    return HornerWithoutFMA<Value, Argument, degree>::Evaluate(
        coefficients, argument);
  }
}

template<typename Value, typename Argument, int degree>
  requires additive_group<Argument>
Derivative<Value, Argument>
Evaluator<Value, Argument, degree>::EvaluateDerivative(
    Coefficients const& coefficients,
    Argument const& argument,
    not_null<Evaluator const*> const evaluator) {
  if (evaluator == Estrin<Value, Argument, degree>::Singleton()) {
    return Estrin<Value, Argument, degree>::EvaluateDerivative(
        coefficients, argument);
  } else if (evaluator == Horner<Value, Argument, degree>::Singleton()) {
    return Horner<Value, Argument, degree>::EvaluateDerivative(
        coefficients, argument);
  } else if (evaluator ==
             EstrinWithoutFMA<Value, Argument, degree>::Singleton()) {
    return EstrinWithoutFMA<Value, Argument, degree>::EvaluateDerivative(
        coefficients, argument);
  } else {
    /*evaluator == HornerWithoutFMA<Value, Argument, degree>::Singleton())*/
    return HornerWithoutFMA<Value, Argument, degree>::EvaluateDerivative(
        coefficients, argument);
  }
}

template<typename Value, typename Argument, int degree>
  requires additive_group<Argument>
void Evaluator<Value, Argument, degree>::WriteToMessage(
    not_null<serialization::PolynomialInMonomialBasis::Evaluator*> message,
    not_null<Evaluator const*> evaluator) {
  serialization::PolynomialInMonomialBasis::Evaluator::Kind kind;
  if (evaluator == Estrin<Value, Argument, degree>::Singleton()) {
    kind = serialization::PolynomialInMonomialBasis::Evaluator::ESTRIN;
  } else if (evaluator == Horner<Value, Argument, degree>::Singleton()) {
    kind = serialization::PolynomialInMonomialBasis::Evaluator::HORNER;
  } else if (evaluator ==
             EstrinWithoutFMA<Value, Argument, degree>::Singleton()) {
    kind =
        serialization::PolynomialInMonomialBasis::Evaluator::ESTRIN_WITHOUT_FMA;
  } else {
    /*evaluator == HornerWithoutFMA<Value, Argument, degree>::Singleton())*/
    kind =
        serialization::PolynomialInMonomialBasis::Evaluator::HORNER_WITHOUT_FMA;
  }
  message->set_kind(kind);
}

template<typename Value, typename Argument, int degree>
  requires additive_group<Argument>
not_null<Evaluator<Value, Argument, degree> const*>
Evaluator<Value, Argument, degree>::ReadFromMessage(
    serialization::PolynomialInMonomialBasis::Evaluator const& message) {
  switch (message.kind()) {
    case serialization::PolynomialInMonomialBasis::Evaluator::ESTRIN:
      return Estrin<Value, Argument, degree>::Singleton();
    case serialization::PolynomialInMonomialBasis::Evaluator::HORNER:
      return Horner<Value, Argument, degree>::Singleton();
    case serialization::PolynomialInMonomialBasis::Evaluator::
        ESTRIN_WITHOUT_FMA:
      return EstrinWithoutFMA<Value, Argument, degree>::Singleton();
    case serialization::PolynomialInMonomialBasis::Evaluator::
        HORNER_WITHOUT_FMA:
      return HornerWithoutFMA<Value, Argument, degree>::Singleton();
  }
  LOG(FATAL) << "Unexpected evaluator " << message.DebugString();
}

template<typename Value, typename Argument, int degree, bool allow_fma>
class EstrinEvaluator : public Evaluator<Value, Argument, degree> {
 public:
  using Coefficients =
      typename Evaluator<Value, Argument, degree>::Coefficients;

  static constexpr Evaluator<Value, Argument, degree> const* Singleton();

  FORCE_INLINE(static) Value
  Evaluate(Coefficients const& coefficients,
           Argument const& argument);
  FORCE_INLINE(static) Derivative<Value, Argument>
  EvaluateDerivative(Coefficients const& coefficients,
                     Argument const& argument);

 private:
  EstrinEvaluator() = default;
};

template<typename Value, typename Argument, int degree, bool allow_fma>
constexpr Evaluator<Value, Argument, degree> const*
EstrinEvaluator<Value, Argument, degree, allow_fma>::Singleton() {
  static constexpr EstrinEvaluator singleton;
  return &singleton;
}

template<typename Value, typename Argument, int degree, bool allow_fma>
Value EstrinEvaluator<Value, Argument, degree, allow_fma>::Evaluate(
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
Derivative<Value, Argument>
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

template<typename Value, typename Argument, int degree, bool allow_fma>
class HornerEvaluator : public Evaluator<Value, Argument, degree> {
 public:
  using Coefficients =
      typename Evaluator<Value, Argument, degree>::Coefficients;

  static constexpr Evaluator<Value, Argument, degree> const* Singleton();

  FORCE_INLINE(static) Value
  Evaluate(Coefficients const& coefficients,
           Argument const& argument);
  FORCE_INLINE(static) Derivative<Value, Argument>
  EvaluateDerivative(Coefficients const& coefficients,
                     Argument const& argument);

 private:
  HornerEvaluator() = default;
};

template<typename Value, typename Argument, int degree, bool allow_fma>
constexpr Evaluator<Value, Argument, degree> const*
HornerEvaluator<Value, Argument, degree, allow_fma>::Singleton() {
  static constexpr HornerEvaluator singleton;
  return &singleton;
}

template<typename Value, typename Argument, int degree, bool allow_fma>
Value HornerEvaluator<Value, Argument, degree, allow_fma>::Evaluate(
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
Derivative<Value, Argument>
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

}  // namespace internal
}  // namespace _polynomial_evaluators
}  // namespace numerics
}  // namespace principia
