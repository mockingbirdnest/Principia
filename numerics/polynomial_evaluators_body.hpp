#pragma once

#include "numerics/polynomial_evaluators.hpp"

#include <cstddef>
#include <tuple>

#include "base/bits.hpp"
#include "numerics/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace _polynomial_evaluators {
namespace internal {

using namespace principia::base::_bits;
using namespace principia::numerics::_elementary_functions;

// Internal helper for Estrin evaluation.  `degree` is the degree of the overall
// polynomial, `low` and `subdegree` defines the subpolynomial that we currently
// evaluate, i.e., the one with a constant term coefficient
// `std::get<low>(coefficients)` and degree `subdegree`.
template<typename Value, typename Argument,
         int degree, bool fma, int low, int subdegree>
struct InternalEstrinEvaluator {
  using Coefficients = internal::Coefficients<Value, Argument, degree>;

  FORCE_INLINE(static) Derivative<Value, Argument, low> Evaluate(
      Coefficients const& coefficients,
      Argument const& argument);
  FORCE_INLINE(static) Derivative<Value, Argument, low> EvaluateDerivative(
      Coefficients const& coefficients,
      Argument const& argument);
};

template<typename Value, typename Argument, int degree, bool fma, int low>
struct InternalEstrinEvaluator<Value, Argument, degree, fma, low, 1> {
  using Coefficients = internal::Coefficients<Value, Argument, degree>;

  FORCE_INLINE(static) Derivative<Value, Argument, low> Evaluate(
      Coefficients const& coefficients,
      Argument const& argument);
  FORCE_INLINE(static) Derivative<Value, Argument, low> EvaluateDerivative(
      Coefficients const& coefficients,
      Argument const& argument);
};

template<typename Value, typename Argument, int degree, bool fma, int low>
struct InternalEstrinEvaluator<Value, Argument, degree, fma, low, 0> {
  using Coefficients = internal::Coefficients<Value, Argument, degree>;

  FORCE_INLINE(static) Derivative<Value, Argument, low> Evaluate(
      Coefficients const& coefficients,
      Argument const& argument);
  FORCE_INLINE(static) Derivative<Value, Argument, low> EvaluateDerivative(
      Coefficients const& coefficients,
      Argument const& argument);
};

template<typename Value, typename Argument, int degree, bool fma, int low>
struct InternalEstrinEvaluator<Value, Argument, degree, fma, low, -1> {
  using Coefficients = internal::Coefficients<Value, Argument, degree>;

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
InternalEstrinEvaluator<Value, Argument, degree, fma, low, subdegree>::
Evaluate(Coefficients const& coefficients,
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

  // Because this function is heavily inlined (including the recursive calls) we
  // rely on common subexpression elimination to detect that we repeatedly
  // compute x², x⁴, etc. and to only compute each power once.  Inspecting the
  // generated code, we verified that the right thing actually happens.
  // This relies on the fact that `Pow` is computed using the Russian peasant
  // algorithm.
  auto const xᵐ = Pow<m>(argument);
  if constexpr (fma) {
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
    return FusedMultiplyAdd(a, xᵐ, b);
  } else {
    return a * xᵐ + b;
  }
}

template<typename Value, typename Argument, int degree, bool fma, int low>
FORCE_INLINE(inline) Derivative<Value, Argument, low>
InternalEstrinEvaluator<Value, Argument, degree, fma, low, 1>::
Evaluate(Coefficients const& coefficients,
         Argument const& argument) {
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
InternalEstrinEvaluator<Value, Argument, degree, fma, low, 0>::
Evaluate(Coefficients const& coefficients,
         Argument const& argument) {
  return std::get<low>(coefficients);
}

template<typename Value, typename Argument, int degree, bool fma, int low>
FORCE_INLINE(inline) Derivative<Value, Argument, low>
InternalEstrinEvaluator<Value, Argument, degree, fma, low, -1>::
Evaluate(Coefficients const& coefficients,
         Argument const& argument) {
  return Derivative<Value, Argument, low>{};
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

template<typename Value, typename Argument, int degree, bool fma, int low>
FORCE_INLINE(inline) Derivative<Value, Argument, low>
InternalEstrinEvaluator<Value, Argument, degree, fma, low, -1>::
EvaluateDerivative(Coefficients const& coefficients,
                   Argument const& argument) {
  return Derivative<Value, Argument, low>{};
}

// Internal helper for Horner evaluation.  `degree` is the degree of the overall
// polynomial, `low` defines the subpolynomial that we currently evaluate, i.e.,
// the one with a constant term coefficient `std::get<low>(coefficients)`.
template<typename Value, typename Argument, int degree, bool fma, int low>
struct InternalHornerEvaluator {
  using Coefficients = internal::Coefficients<Value, Argument, degree>;

  FORCE_INLINE(static) Derivative<Value, Argument, low>
  Evaluate(Coefficients const& coefficients,
           Argument const& argument);
  FORCE_INLINE(static) Derivative<Value, Argument, low>
  EvaluateDerivative(Coefficients const& coefficients,
                     Argument const& argument);
};

template<typename Value, typename Argument, int degree, bool fma>
struct InternalHornerEvaluator<Value, Argument, degree, fma, degree> {
  using Coefficients = internal::Coefficients<Value, Argument, degree>;

  FORCE_INLINE(static) Derivative<Value, Argument, degree>
  Evaluate(Coefficients const& coefficients,
           Argument const& argument);
  FORCE_INLINE(static) Derivative<Value, Argument, degree>
  EvaluateDerivative(Coefficients const& coefficients,
                     Argument const& argument);
};

template<typename Value, typename Argument, int degree, bool fma, int low>
FORCE_INLINE(inline) Derivative<Value, Argument, low>
InternalHornerEvaluator<Value, Argument, degree, fma, low>::
Evaluate(Coefficients const& coefficients,
         Argument const& argument) {
  if constexpr (low > degree) {
    return Derivative<Value, Argument, low>{};
  } else {
    auto const& x = argument;
    auto const a =
        InternalHornerEvaluator<Value, Argument, degree, fma, low + 1>::
            Evaluate(coefficients, argument);
    auto const& b = std::get<low>(coefficients);
    if constexpr (fma) {
      return FusedMultiplyAdd(a, x, b);
    } else {
      return a * x + b;
    }
  }
}

template<typename Value, typename Argument, int degree, bool fma, int low>
FORCE_INLINE(inline) Derivative<Value, Argument, low>
InternalHornerEvaluator<Value, Argument, degree, fma, low>::
EvaluateDerivative(Coefficients const& coefficients,
                   Argument const& argument) {
  if constexpr (low > degree) {
    return Derivative<Value, Argument, low>{};
  } else {
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
}

template<typename Value, typename Argument, int degree, bool fma>
FORCE_INLINE(inline) Derivative<Value, Argument, degree>
InternalHornerEvaluator<Value, Argument, degree, fma, degree>::
Evaluate(Coefficients const& coefficients,
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

template<typename Value, typename Argument, int degree,
         FMAPolicy fma_policy, FMAPresence fma_presence>
Value EstrinEvaluator<Value, Argument, degree, fma_policy, fma_presence>::
Evaluate(Coefficients const& coefficients,
         Argument const& argument) {
  if constexpr (fma_policy == FMAPolicy::Auto) {
    if ((fma_presence == FMAPresence::Present ||
         (fma_presence == FMAPresence::Unknown && CanUseHardwareFMA))) {
      using InternalEvaluator = InternalEstrinEvaluator<Value,
                                                        Argument,
                                                        degree,
                                                        /*fma=*/true,
                                                        /*low=*/0,
                                                        /*subdegree=*/degree>;
      return InternalEvaluator::Evaluate(coefficients, argument);
    }
  }
  using InternalEvaluator = InternalEstrinEvaluator<Value,
                                                    Argument,
                                                    degree,
                                                    /*fma=*/false,
                                                    /*low=*/0,
                                                    /*subdegree=*/degree>;
  return InternalEvaluator::Evaluate(coefficients, argument);
}

template<typename Value, typename Argument, int degree,
         FMAPolicy fma_policy, FMAPresence fma_presence>
Derivative<Value, Argument>
EstrinEvaluator<Value, Argument, degree, fma_policy, fma_presence>::
EvaluateDerivative(Coefficients const& coefficients,
                   Argument const& argument) {
  if constexpr (fma_policy == FMAPolicy::Auto) {
    if (fma_presence == FMAPresence::Present ||
        (fma_presence == FMAPresence::Unknown && CanUseHardwareFMA)) {
        using InternalEvaluator =
            InternalEstrinEvaluator<Value,
                                    Argument,
                                    degree,
                                    /*fma=*/true,
                                    /*low=*/1,
                                    /*subdegree=*/degree - 1>;
        return InternalEvaluator::EvaluateDerivative(coefficients, argument);
      }
  }
  using InternalEvaluator = InternalEstrinEvaluator<Value,
                                                    Argument,
                                                    degree,
                                                    /*fma=*/false,
                                                    /*low=*/1,
                                                    /*subdegree=*/degree - 1>;
  return InternalEvaluator::EvaluateDerivative(coefficients, argument);
}

template<typename Value, typename Argument, int degree,
         FMAPolicy fma_policy, FMAPresence fma_presence>
void EstrinEvaluator<Value, Argument, degree, fma_policy, fma_presence>::
EvaluateWithDerivative(Coefficients const& coefficients,
                       Argument const& argument,
                       Value& value,
                       Derivative<Value, Argument>& derivative) {
  if constexpr (fma_policy == FMAPolicy::Auto) {
    if (fma_presence == FMAPresence::Present ||
        (fma_presence == FMAPresence::Unknown && CanUseHardwareFMA)) {
        value = InternalEstrinEvaluator<
            Value,
            Argument,
            degree,
            /*fma=*/true,
            /*low=*/0,
            /*subdegree=*/degree>::Evaluate(coefficients, argument);
        derivative =
            InternalEstrinEvaluator<Value,
                                    Argument,
                                    degree,
                                    /*fma=*/true,
                                    /*low=*/1,
                                    /*subdegree=*/degree -
                                        1>::EvaluateDerivative(coefficients,
                                                               argument);
        return;
      }
  }
  value = InternalEstrinEvaluator<Value,
                                  Argument,
                                  degree,
                                  /*fma=*/false,
                                  /*low=*/0,
                                  /*subdegree=*/degree>::Evaluate(coefficients,
                                                                  argument);
  derivative = InternalEstrinEvaluator<Value,
                                       Argument,
                                       degree,
                                       /*fma=*/false,
                                       /*low=*/1,
                                       /*subdegree=*/degree -
                                           1>::EvaluateDerivative(coefficients,
                                                                  argument);
}

template<typename Value, typename Argument, int degree,
         FMAPolicy fma_policy, FMAPresence fma_presence>
void EstrinEvaluator<Value, Argument, degree, fma_policy, fma_presence>::
WriteToMessage(
    not_null<serialization::PolynomialInMonomialBasis::Evaluator*> message) {
  switch (fma_policy) {
    case FMAPolicy::Auto:
      message->set_kind(
          serialization::PolynomialInMonomialBasis::Evaluator::ESTRIN);
      break;
    case FMAPolicy::Disallow:
      message->set_kind(serialization::PolynomialInMonomialBasis::Evaluator::
                            ESTRIN_WITHOUT_FMA);
      break;
  }
}


template<typename Value, typename Argument, int degree,
         FMAPolicy fma_policy, FMAPresence fma_presence>
Value HornerEvaluator<Value, Argument, degree, fma_policy, fma_presence>::
Evaluate(Coefficients const& coefficients,
         Argument const& argument) {
  if constexpr (fma_policy == FMAPolicy::Auto) {
    if (fma_presence == FMAPresence::Present ||
        (fma_presence == FMAPresence::Unknown && CanUseHardwareFMA)) {
      return InternalHornerEvaluator<Value,
                                     Argument,
                                     degree,
                                     /*fma=*/true,
                                     /*low=*/0>::Evaluate(coefficients,
                                                          argument);
    }
  }
  return InternalHornerEvaluator<Value,
                                 Argument,
                                 degree,
                                 /*fma=*/false,
                                 /*low=*/0>::Evaluate(coefficients, argument);
}

template<typename Value, typename Argument, int degree,
         FMAPolicy fma_policy, FMAPresence fma_presence>
Derivative<Value, Argument>
HornerEvaluator<Value, Argument, degree, fma_policy, fma_presence>::
    EvaluateDerivative(Coefficients const& coefficients,
                       Argument const& argument) {
  if constexpr (fma_policy == FMAPolicy::Auto) {
    if (fma_presence == FMAPresence::Present ||
        (fma_presence == FMAPresence::Unknown && CanUseHardwareFMA)) {
      return InternalHornerEvaluator<
          Value,
          Argument,
          degree,
          /*fma=*/true,
          /*low=*/1>::EvaluateDerivative(coefficients, argument);
    }
  }
  return InternalHornerEvaluator<Value,
                                 Argument,
                                 degree,
                                 /*fma=*/false,
                                 /*low=*/1>::EvaluateDerivative(coefficients,
                                                                argument);
}

template<typename Value, typename Argument, int degree,
         FMAPolicy fma_policy, FMAPresence fma_presence>
void HornerEvaluator<Value, Argument, degree, fma_policy, fma_presence>::
EvaluateWithDerivative(Coefficients const& coefficients,
                       Argument const& argument,
                       Value& value,
                           Derivative<Value, Argument>& derivative) {
  if constexpr (fma_policy == FMAPolicy::Auto) {
    if (fma_presence == FMAPresence::Present ||
        (fma_presence == FMAPresence::Unknown && CanUseHardwareFMA)) {
      value =
          InternalHornerEvaluator<Value,
                                  Argument,
                                  degree,
                                  /*fma=*/true,
                                  /*low=*/0>::Evaluate(coefficients, argument);
      derivative =
          InternalHornerEvaluator<Value,
                                  Argument,
                                  degree,
                                  /*fma=*/true,
                                  /*low=*/1>::EvaluateDerivative(coefficients,
                                                                 argument);
    }
  }
  value = InternalHornerEvaluator<Value,
                                  Argument,
                                  degree,
                                  /*fma=*/false,
                                  /*low=*/0>::Evaluate(coefficients, argument);
  derivative =
      InternalHornerEvaluator<Value,
                              Argument,
                              degree,
                              /*fma=*/false,
                              /*low=*/1>::EvaluateDerivative(coefficients,
                                                             argument);
}

template<typename Value, typename Argument, int degree,
         FMAPolicy fma_policy, FMAPresence fma_presence>
void HornerEvaluator<Value, Argument, degree, fma_policy, fma_presence>::
WriteToMessage(
    not_null<serialization::PolynomialInMonomialBasis::Evaluator*> message) {
  switch (fma_policy) {
    case FMAPolicy::Auto:
      message->set_kind(
          serialization::PolynomialInMonomialBasis::Evaluator::HORNER);
      break;
    case FMAPolicy::Disallow:
      message->set_kind(serialization::PolynomialInMonomialBasis::Evaluator::
                            HORNER_WITHOUT_FMA);
      break;
  }
}

}  // namespace internal
}  // namespace _polynomial_evaluators
}  // namespace numerics
}  // namespace principia
