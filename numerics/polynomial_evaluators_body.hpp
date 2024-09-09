#pragma once

#include "numerics/polynomial_evaluators.hpp"

#include <cstddef>
#include <tuple>

#include "base/bits.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace _polynomial_evaluators {
namespace internal {

using namespace principia::base::_bits;
using namespace principia::quantities::_elementary_functions;


// Internal helper for Estrin evaluation.  `degree` is the degree of the overall
// polynomial, `low` and `subdegree` defines the subpolynomial that we currently
// evaluate, i.e., the one with a constant term coefficient
// `std::get<low>(coefficients)` and degree `subdegree`.
template<typename Value, typename Argument,
         int degree, bool fma, int low, int subdegree>
struct InternalEstrinEvaluator {
  using Coefficients =
      typename Evaluator<Value, Argument, degree>::Coefficients;

  FORCE_INLINE(static) Derivative<Value, Argument, low> Evaluate(
      Coefficients const& coefficients,
      Argument const& argument);
  FORCE_INLINE(static) Derivative<Value, Argument, low> EvaluateDerivative(
      Coefficients const& coefficients,
      Argument const& argument);
};

template<typename Value, typename Argument, int degree, bool fma, int low>
struct InternalEstrinEvaluator<Value, Argument, degree, fma, low, 1> {
  using Coefficients =
      typename Evaluator<Value, Argument, degree>::Coefficients;

  FORCE_INLINE(static) Derivative<Value, Argument, low> Evaluate(
      Coefficients const& coefficients,
      Argument const& argument);
  FORCE_INLINE(static) Derivative<Value, Argument, low> EvaluateDerivative(
      Coefficients const& coefficients,
      Argument const& argument);
};

template<typename Value, typename Argument, int degree, bool fma, int low>
struct InternalEstrinEvaluator<Value, Argument, degree, fma, low, 0> {
  using Coefficients =
      typename Evaluator<Value, Argument, degree>::Coefficients;

  FORCE_INLINE(static) Derivative<Value, Argument, low> Evaluate(
      Coefficients const& coefficients,
      Argument const& argument);
  FORCE_INLINE(static) Derivative<Value, Argument, low> EvaluateDerivative(
      Coefficients const& coefficients,
      Argument const& argument);
};

template<typename Value, typename Argument,
         int degree, bool fma, int low, int subdegree>
FORCE_INLINE(inline) Derivative<Value, Argument, low>
InternalEstrinEvaluator<Value, Argument, degree, fma, low, subdegree>::Evaluate(
    Coefficients const& coefficients,
    Argument const& argument) {
  static_assert(subdegree >= 2,
                "Unexpected subdegree in InternalEstrinEvaluator::Evaluate");
  constexpr int m = PowerOf2Le(subdegree);

  auto const a =
      InternalEstrinEvaluator<Value, Argument,
                              degree, fma, low + m, subdegree - m>::
          Evaluate(coefficients, argument);
  auto const b =
      InternalEstrinEvaluator<Value, Argument, degree, fma, low, m - 1>::
          Evaluate(coefficients, argument);

  // Because this function is heavily inline (including the recursive calls) we
  // rely on common subexpression elimination to detect that we repeatedly
  // compute x², x⁴, etc. and to only compute each power once.  Inspecting the
  // generated code, we verified that the right thing actually happens.
  // This relies on the fact that `Pow` is computed using the Russian peasant
  // algorithm.
  auto const xᵐ = Pow<m>(argument);
  if constexpr (fma) {
    using quantities::_elementary_functions::FusedMultiplyAdd;
    return FusedMultiplyAdd(a, xᵐ, b);
  } else {
    return a * xᵐ + b;
  }
}

template<typename Value, typename Argument,
         int degree, bool fma, int low, int subdegree>
FORCE_INLINE(inline) Derivative<Value, Argument, low>
InternalEstrinEvaluator<Value, Argument, degree, fma, low, subdegree>::
EvaluateDerivative(Coefficients const& coefficients,
                   Argument const& argument) {
  static_assert(subdegree >= 2,
                "Unexpected subdegree in InternalEstrinEvaluator::"
                "EvaluateDerivative");
  constexpr int m = PowerOf2Le(subdegree);

  auto const a =
      InternalEstrinEvaluator<Value, Argument,
                              degree, fma, low + m, subdegree - m>::
          EvaluateDerivative(coefficients, argument);
  auto const b =
      InternalEstrinEvaluator<Value, Argument, degree, fma, low, m - 1>::
          EvaluateDerivative(coefficients, argument);

  auto const xᵐ = Pow<m>(argument);
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
    Argument const& argument) {
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
    Argument const& argument) {
  return std::get<low>(coefficients);
}

template<typename Value, typename Argument, int degree, bool fma, int low>
FORCE_INLINE(inline) Derivative<Value, Argument, low>
InternalEstrinEvaluator<Value, Argument, degree, fma, low, 1>::
EvaluateDerivative(Coefficients const& coefficients,
                    Argument const& argument) {
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
                    Argument const& argument) {
  return low * std::get<low>(coefficients);
}

// Internal helper for Horner evaluation.  `degree` is the degree of the overall
// polynomial, `low` defines the subpolynomial that we currently evaluate, i.e.,
// the one with a constant term coefficient `std::get<low>(coefficients)`.
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
  if constexpr (degree <= 2) {
    // Horner and Estrin are identical.
    if (evaluator == Estrin<Value, Argument, degree>::Singleton() ||
        evaluator == Horner<Value, Argument, degree>::Singleton()) {
      return Horner<Value, Argument, degree>::Evaluate(coefficients, argument);
    } else {
      /*evaluator == EstrinWithoutFMA<Value, Argument, degree>::Singleton() ||
        evaluator == HornerWithoutFMA<Value, Argument, degree>::Singleton())*/
      return HornerWithoutFMA<Value, Argument, degree>::Evaluate(
          coefficients, argument);
    }
  } else {
    if (evaluator == Estrin<Value, Argument, degree>::Singleton()) {
      return Estrin<Value, Argument, degree>::Evaluate(coefficients, argument);
    } else if (evaluator == Horner<Value, Argument, degree>::Singleton()) {
      return Horner<Value, Argument, degree>::Evaluate(coefficients, argument);
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
}

template<typename Value, typename Argument, int degree>
  requires additive_group<Argument>
Derivative<Value, Argument>
Evaluator<Value, Argument, degree>::EvaluateDerivative(
    Coefficients const& coefficients,
    Argument const& argument,
    not_null<Evaluator const*> const evaluator) {
  // For some reason using a simpler control flow for `degree <= 3`, like we do
  // above, degrades the performance of `Evaluate`.  So let's not do that.
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

template<typename Value, typename Argument, int degree, FMAPolicy fma_policy>
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

template<typename Value, typename Argument, int degree, FMAPolicy fma_policy>
constexpr Evaluator<Value, Argument, degree> const*
EstrinEvaluator<Value, Argument, degree, fma_policy>::Singleton() {
  static constexpr EstrinEvaluator singleton;
  return &singleton;
}

template<typename Value, typename Argument, int degree, FMAPolicy fma_policy>
Value EstrinEvaluator<Value, Argument, degree, fma_policy>::Evaluate(
    Coefficients const& coefficients,
    Argument const& argument) {
  if ((fma_policy == FMAPolicy::Force && CanEmitFMAInstructions) ||
      (fma_policy == FMAPolicy::Auto && UseHardwareFMA)) {
    using InternalEvaluator = InternalEstrinEvaluator<Value,
                                                      Argument,
                                                      degree,
                                                      /*fma=*/true,
                                                      /*low=*/0,
                                                      /*subdegree=*/degree>;
    return InternalEvaluator::Evaluate(coefficients, argument);
  } else {
    using InternalEvaluator = InternalEstrinEvaluator<Value,
                                                      Argument,
                                                      degree,
                                                      /*fma=*/false,
                                                      /*low=*/0,
                                                      /*subdegree=*/degree>;
    return InternalEvaluator::Evaluate(coefficients, argument);
  }
}

template<typename Value, typename Argument, int degree, FMAPolicy fma_policy>
Derivative<Value, Argument>
EstrinEvaluator<Value, Argument, degree, fma_policy>::EvaluateDerivative(
    Coefficients const& coefficients,
    Argument const& argument) {
  if constexpr (degree == 0) {
    return Derivative<Value, Argument>{};
  } else if ((fma_policy == FMAPolicy::Force && CanEmitFMAInstructions) ||
             (fma_policy == FMAPolicy::Auto && UseHardwareFMA)) {
    using InternalEvaluator =
        InternalEstrinEvaluator<Value,
                                Argument,
                                degree,
                                true,
                                /*low=*/1,
                                /*subdegree=*/degree - 1>;
    return InternalEvaluator::EvaluateDerivative(coefficients, argument);
  } else {
    using InternalEvaluator =
        InternalEstrinEvaluator<Value,
                                Argument,
                                degree,
                                false,
                                /*low=*/1,
                                /*subdegree=*/degree - 1>;
    return InternalEvaluator::EvaluateDerivative(coefficients, argument);
  }
}

template<typename Value, typename Argument, int degree, FMAPolicy fma_policy>
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

template<typename Value, typename Argument, int degree, FMAPolicy fma_policy>
constexpr Evaluator<Value, Argument, degree> const*
HornerEvaluator<Value, Argument, degree, fma_policy>::Singleton() {
  static constexpr HornerEvaluator singleton;
  return &singleton;
}

template<typename Value, typename Argument, int degree, FMAPolicy fma_policy>
Value HornerEvaluator<Value, Argument, degree, fma_policy>::Evaluate(
    Coefficients const& coefficients,
    Argument const& argument) {
  if ((fma_policy == FMAPolicy::Force && CanEmitFMAInstructions) ||
      (fma_policy == FMAPolicy::Auto && UseHardwareFMA)) {
    return InternalHornerEvaluator<
        Value, Argument, degree, /*fma=*/true, /*low=*/0>::Evaluate(
            coefficients, argument);
  } else {
    return InternalHornerEvaluator<
        Value, Argument, degree, /*fma=*/false, /*low=*/0>::Evaluate(
            coefficients, argument);
  }
}

template<typename Value, typename Argument, int degree, FMAPolicy fma_policy>
Derivative<Value, Argument>
HornerEvaluator<Value, Argument, degree, fma_policy>::EvaluateDerivative(
    Coefficients const& coefficients,
    Argument const& argument) {
  if constexpr (degree == 0) {
    return Derivative<Value, Argument>{};
  } else if ((fma_policy == FMAPolicy::Force && CanEmitFMAInstructions) ||
             (fma_policy == FMAPolicy::Auto && UseHardwareFMA)) {
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
