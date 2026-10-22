#pragma once

#include "numerics/polynomial_evaluators.hpp"

#include <cstddef>
#include <tuple>

#include "base/bits.hpp"
#include "numerics/elementary_functions.hpp"
#include "numerics/wide.hpp"

namespace principia {
namespace numerics {
namespace _polynomial_evaluators {
namespace internal {

using namespace principia::base::_bits;
using namespace principia::numerics::_elementary_functions;
using namespace principia::numerics::_wide;

template<typename Value, typename Argument>
concept wide_enough_is_wide =
    requires(Derivative<Value, Argument> a, Difference<Argument> x) {
      { a * Wide(x) };
    };

template<typename Value, typename Argument>
struct WideEnoughGenerator {
  using type = Argument;
};

template<typename Value, typename Argument>
  requires wide_enough_is_wide<Value, Argument>
struct WideEnoughGenerator<Value, Argument> {
  using type = Wide<Argument>;
};

template<typename Value, typename Argument>
using WideEnough = typename WideEnoughGenerator<Value, Argument>::type;

template<typename Value, typename Argument, int low, int lowʹ>
struct ValueAndDerivative {
  Derivative<Value, Argument, low> value;
  //TODO(phl)what is the derivative exactly?
  Derivative<Value, Argument, lowʹ> derivative;
};

// This will look weird.  This function saves a meaningful number of cycles
// (roughly 10 for degree 10, 22 for degree 17) when evaluating the value and
// derivative of `RelativeDegreesOfFreedom`-valued polynomials (see the
// nanobenchmarks).  It saves a multiplication by 1.0 and a memory access to
// load 2.0, but that's not where the performance gain comes from.  Apparently,
// there is a lot less spilling of registers to memory when using this function
// (xmm6-xmm11 instead of xmm6-xmm15), and also no temporaries touching memory.
// Note that passing `x` by reference kills the performance gain.
// I am leaving a `static_assert` because we'll need to recheck (using the
// nanobenchmarks) that this function is still beneficial after any change of
// compiler.
#if 0 && PRINCIPIA_COMPILER_MSVC
template<int m, typename Value>
  requires(m > 0)
Value Mul(Value const x) {
  static_assert(_MSC_FULL_VER == 195'236'725 || _MSC_FULL_VER == 195'236'815,
                "Check if the Mul function is still beneficial");
  if constexpr (m == 1) {
    return x;
  } else if constexpr (m == 2) {
    return x + x;
  } else {
    return m * x;
  }
}
#define MUL(m, x) Mul<m>(x)
#else
#define MUL(m, x) ((m) * (x))
#endif

// Internal helper for Estrin evaluation.  `degree` is the degree of the overall
// polynomial, `low` and `degree` defines the subpolynomial that we currently
// evaluate, i.e., the one with a constant term coefficient
// `get<low>(coefficients)` and degree `degree`.
template<typename Value, typename Argument, typename Coefficients,
         bool fma, int low, int degree>
struct InternalEstrinEvaluator {
  using WideEnoughArgument = WideEnough<Value, Argument>;

  FORCE_INLINE static Derivative<Value, Argument, low> Evaluate(
      Coefficients const& coefficients,
      WideEnoughArgument const& argument);
  FORCE_INLINE static Derivative<Value, Argument, low> EvaluateDerivative(
      Coefficients const& coefficients,
      WideEnoughArgument const& argument);
};

template<typename Value, typename Argument, typename Coefficients,
         bool fma, int low, int degree, int lowʹ, int degreeʹ>
struct InternalEstrinEvaluator2 {
  using WideEnoughArgument = WideEnough<Value, Argument>;

  FORCE_INLINE static ValueAndDerivative<Value, Argument, low, lowʹ>
  EvaluateWithDerivative(Coefficients const& coefficients,
                         WideEnoughArgument const& argument);
};

template<typename Value, typename Argument, typename Coefficients,
         bool fma, int low>
struct InternalEstrinEvaluator<Value, Argument, Coefficients, fma, low, 1> {
  using WideEnoughArgument = WideEnough<Value, Argument>;

  FORCE_INLINE static Derivative<Value, Argument, low> Evaluate(
      Coefficients const& coefficients,
      WideEnoughArgument const& argument);
  FORCE_INLINE static Derivative<Value, Argument, low> EvaluateDerivative(
      Coefficients const& coefficients,
      WideEnoughArgument const& argument);
};

template<typename Value, typename Argument, typename Coefficients,
         bool fma, int low, int lowʹ>
struct InternalEstrinEvaluator2<Value, Argument, Coefficients,
                                fma, low, 1, lowʹ, 1> {
  using WideEnoughArgument = WideEnough<Value, Argument>;

  FORCE_INLINE static ValueAndDerivative<Value, Argument, low, lowʹ>
  EvaluateWithDerivative(Coefficients const& coefficients,
                         WideEnoughArgument const& argument);
};

template<typename Value, typename Argument, typename Coefficients,
         bool fma, int low, int lowʹ>
struct InternalEstrinEvaluator2<Value, Argument, Coefficients,
                                fma, low, 1, lowʹ, 0> {
  using WideEnoughArgument = WideEnough<Value, Argument>;

  FORCE_INLINE static ValueAndDerivative<Value, Argument, low, lowʹ>
  EvaluateWithDerivative(Coefficients const& coefficients,
                         WideEnoughArgument const& argument);
};

template<typename Value, typename Argument, typename Coefficients,
         bool fma, int low>
struct InternalEstrinEvaluator<Value, Argument, Coefficients, fma, low, 0> {
  using WideEnoughArgument = WideEnough<Value, Argument>;

  FORCE_INLINE static Derivative<Value, Argument, low> Evaluate(
      Coefficients const& coefficients,
      WideEnoughArgument const& argument);
  FORCE_INLINE static Derivative<Value, Argument, low> EvaluateDerivative(
      Coefficients const& coefficients,
      WideEnoughArgument const& argument);
};

template<typename Value, typename Argument, typename Coefficients,
         bool fma, int low, int lowʹ>
struct InternalEstrinEvaluator2<Value, Argument, Coefficients,
                                fma, low, 0, lowʹ, 0> {
  using WideEnoughArgument = WideEnough<Value, Argument>;

  FORCE_INLINE static ValueAndDerivative<Value, Argument, low, lowʹ>
  EvaluateWithDerivative(Coefficients const& coefficients,
                         WideEnoughArgument const& argument);
};

template<typename Value, typename Argument, typename Coefficients,
         bool fma, int low>
struct InternalEstrinEvaluator<Value, Argument, Coefficients, fma, low, -1> {
  using WideEnoughArgument = WideEnough<Value, Argument>;

  FORCE_INLINE static Derivative<Value, Argument, low> Evaluate(
      Coefficients const& coefficients,
      WideEnoughArgument const& argument);
  FORCE_INLINE static Derivative<Value, Argument, low> EvaluateDerivative(
      Coefficients const& coefficients,
      WideEnoughArgument const& argument);
};

template<typename Value, typename Argument, typename Coefficients,
         bool fma, int low, int degree>
FORCE_INLINE Derivative<Value, Argument, low>
InternalEstrinEvaluator<Value, Argument, Coefficients, fma, low, degree>::
Evaluate(Coefficients const& coefficients,
         WideEnoughArgument const& argument) {
  static_assert(degree >= 2,
                "Unexpected degree in InternalEstrinEvaluator::Evaluate");
  //TODO(phl)Comment m.
  constexpr int m₀ = PowerOf2Lt(degree);
  constexpr int m = degree - m₀ < m₀ - 1 ? m₀ / 2 : m₀;

  auto const a =
      InternalEstrinEvaluator<Value, Argument, Coefficients,
                              fma, low + m, degree - m>::
          Evaluate(coefficients, argument);
  auto const b =
      InternalEstrinEvaluator<Value, Argument, Coefficients,
                              fma, low, m - 1>::
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

template<typename Value, typename Argument, typename Coefficients,
         bool fma, int low, int degree>
FORCE_INLINE Derivative<Value, Argument, low>
InternalEstrinEvaluator<Value, Argument, Coefficients, fma, low, degree>::
EvaluateDerivative(Coefficients const& coefficients,
                   WideEnoughArgument const& argument) {
  static_assert(degree >= 2,
                "Unexpected degree in InternalEstrinEvaluator::"
                "EvaluateDerivative");
  constexpr int m₀ = PowerOf2Lt(degree);
  constexpr int m = degree - m₀ < m₀ - 1 ? m₀ / 2 : m₀;

  auto const a =
      InternalEstrinEvaluator<Value, Argument, Coefficients,
                              fma, low + m, degree - m>::
          EvaluateDerivative(coefficients, argument);
  auto const b =
      InternalEstrinEvaluator<Value, Argument, Coefficients,
                              fma, low, m - 1>::
          EvaluateDerivative(coefficients, argument);

  auto const xᵐ = Pow<m>(argument);
  if constexpr (fma) {
    return FusedMultiplyAdd(a, xᵐ, b);
  } else {
    return a * xᵐ + b;
  }
}

template<typename Value, typename Argument, typename Coefficients,
         bool fma, int low, int degree, int lowʹ, int degreeʹ>
FORCE_INLINE ValueAndDerivative<Value, Argument, low, lowʹ>
InternalEstrinEvaluator2<Value, Argument, Coefficients,
                         fma, low, degree, lowʹ, degreeʹ>::
EvaluateWithDerivative(Coefficients const& coefficients,
                       WideEnoughArgument const& argument) {
  static_assert(degree >= 2,
                "Unexpected degree in InternalEstrinEvaluator::"
                "EvaluateWithDerivative");
  static_assert(degree >= degreeʹ);
  static_assert(low <= lowʹ);
  constexpr int m₀ = PowerOf2Lt(degree);
  constexpr int m = degree - m₀ < m₀ - 1 ? m₀ / 2 : m₀;

  auto const a =
      InternalEstrinEvaluator2<Value, Argument, Coefficients,
                              fma, low + m, degree - m, lowʹ + m, degreeʹ -m>::
          EvaluateWithDerivative(coefficients, argument);
  auto const b =
      InternalEstrinEvaluator2<Value, Argument, Coefficients,
                              fma, low, m - 1, lowʹ, m - 1>::
          EvaluateWithDerivative(coefficients, argument);

  auto const& x = argument;
  auto const xᵐ = Pow<m>(x);
  if constexpr (fma) {
    return {.value = FusedMultiplyAdd(a.value, xᵐ, b.value),
            .derivative = FusedMultiplyAdd(a.derivative, xᵐ, b.derivative)};
  } else {
    return {.value = a.value * xᵐ + b.value,
            .derivative = a.derivative * xᵐ + b.derivative};
  }
}

template<typename Value, typename Argument, typename Coefficients,
         bool fma, int low>
FORCE_INLINE Derivative<Value, Argument, low>
InternalEstrinEvaluator<Value, Argument, Coefficients, fma, low, 1>::
Evaluate(Coefficients const& coefficients,
         WideEnoughArgument const& argument) {
  auto const& x = argument;
  auto const& a = get<low + 1>(coefficients);
  auto const& b = get<low>(coefficients);
  if constexpr (fma) {
    return FusedMultiplyAdd(a, x, b);
  } else {
    return a * x + b;
  }
}

template<typename Value, typename Argument, typename Coefficients,
         bool fma, int low>
FORCE_INLINE Derivative<Value, Argument, low>
InternalEstrinEvaluator<Value, Argument, Coefficients, fma, low, 1>::
EvaluateDerivative(Coefficients const& coefficients,
                   WideEnoughArgument const& argument) {
  auto const& x = argument;
  auto const a = MUL(low + 1, get<low + 1>(coefficients));
  auto const b = MUL(low, get<low>(coefficients));
  if constexpr (fma) {
    return FusedMultiplyAdd(a, x, b);
  } else {
    return a * x + b;
  }
}

template<typename Value, typename Argument, typename Coefficients,
         bool fma, int low, int lowʹ>
FORCE_INLINE ValueAndDerivative<Value, Argument, low, lowʹ>
InternalEstrinEvaluator2<Value, Argument, Coefficients, fma, low, 1, lowʹ, 1>::
EvaluateWithDerivative(Coefficients const& coefficients,
                       WideEnoughArgument const& argument) {
  auto const& x = argument;
  auto const& a = get<low + 1>(coefficients);
  auto const& b = get<low>(coefficients);
  auto const aʹ = MUL(lowʹ + 1, get<lowʹ + 1>(coefficients));
  auto const bʹ = MUL(lowʹ, get<lowʹ>(coefficients));
  if constexpr (fma) {
    return {.value = FusedMultiplyAdd(a, x, b),
            .derivative = FusedMultiplyAdd(aʹ, x, bʹ)};
  } else {
    return {.value = a * x + b, .derivative = aʹ * x + bʹ};
  }
}

template<typename Value, typename Argument, typename Coefficients,
         bool fma, int low, int lowʹ>
FORCE_INLINE ValueAndDerivative<Value, Argument, low, lowʹ>
InternalEstrinEvaluator2<Value, Argument, Coefficients, fma, low, 1, lowʹ, 0>::
EvaluateWithDerivative(Coefficients const& coefficients,
                        WideEnoughArgument const& argument) {
  auto const& x = argument;
  auto const& a = get<low + 1>(coefficients);
  auto const& b = get<low>(coefficients);
  auto const bʹ = MUL(lowʹ, get<lowʹ>(coefficients));
  if constexpr (fma) {
    return {.value = FusedMultiplyAdd(a, x, b), .derivative = bʹ};
  } else {
    return {.value = a * x + b, .derivative = bʹ};
  }
}

template<typename Value, typename Argument, typename Coefficients,
         bool fma, int low>
FORCE_INLINE Derivative<Value, Argument, low>
InternalEstrinEvaluator<Value, Argument, Coefficients, fma, low, 0>::
Evaluate(Coefficients const& coefficients,
         WideEnoughArgument const& /*argument*/) {
  return get<low>(coefficients);
}

template<typename Value, typename Argument, typename Coefficients,
         bool fma, int low>
FORCE_INLINE Derivative<Value, Argument, low>
InternalEstrinEvaluator<Value, Argument, Coefficients, fma, low, 0>::
EvaluateDerivative(Coefficients const& coefficients,
                   WideEnoughArgument const& /*argument*/) {
  return MUL(low, get<low>(coefficients));
}

template<typename Value, typename Argument, typename Coefficients,
         bool fma, int low, int lowʹ>
FORCE_INLINE ValueAndDerivative<Value, Argument, low, lowʹ>
InternalEstrinEvaluator2<Value, Argument, Coefficients, fma, low, 0, lowʹ, 0>::
EvaluateWithDerivative(Coefficients const& coefficients,
                        WideEnoughArgument const& /*argument*/) {
  auto const& b = get<low>(coefficients);
  auto const bʹ = MUL(lowʹ, get<lowʹ>(coefficients));
  return {.value = b, .derivative = bʹ};
}

template<typename Value, typename Argument, typename Coefficients,
         bool fma, int low>
FORCE_INLINE Derivative<Value, Argument, low>
InternalEstrinEvaluator<Value, Argument, Coefficients, fma, low, -1>::
Evaluate(Coefficients const& /*coefficients*/,
         WideEnoughArgument const& /*argument*/) {
  return Derivative<Value, Argument, low>{};
}

template<typename Value, typename Argument, typename Coefficients,
         bool fma, int low>
FORCE_INLINE Derivative<Value, Argument, low>
InternalEstrinEvaluator<Value, Argument, Coefficients, fma, low, -1>::
EvaluateDerivative(Coefficients const& /*coefficients*/,
                   WideEnoughArgument const& /*argument*/) {
  return Derivative<Value, Argument, low>{};
}

// Internal helper for Horner evaluation.  `degree` is the degree of the overall
// polynomial, `low` defines the subpolynomial that we currently evaluate, i.e.,
// the one with a constant term coefficient `get<low>(coefficients)`.
template<typename Value, typename Argument, int degree, bool fma, int low>
struct InternalHornerEvaluator {
  using Coefficients = internal::Coefficients<Value, Argument, degree>;
  using WideEnoughArgument = WideEnough<Value, Argument>;

  FORCE_INLINE static Derivative<Value, Argument, low>
  Evaluate(Coefficients const& coefficients,
           WideEnoughArgument const& argument);
  FORCE_INLINE static Derivative<Value, Argument, low>
  EvaluateDerivative(Coefficients const& coefficients,
                     WideEnoughArgument const& argument);
};

template<typename Value, typename Argument, int degree, bool fma>
struct InternalHornerEvaluator<Value, Argument, degree, fma, degree> {
  using Coefficients = internal::Coefficients<Value, Argument, degree>;
  using WideEnoughArgument = WideEnough<Value, Argument>;

  FORCE_INLINE static Derivative<Value, Argument, degree>
  Evaluate(Coefficients const& coefficients,
           WideEnoughArgument const& argument);
  FORCE_INLINE static Derivative<Value, Argument, degree>
  EvaluateDerivative(Coefficients const& coefficients,
                     WideEnoughArgument const& argument);
};

template<typename Value, typename Argument, int degree, bool fma, int low>
FORCE_INLINE Derivative<Value, Argument, low>
InternalHornerEvaluator<Value, Argument, degree, fma, low>::
Evaluate(Coefficients const& coefficients,
         WideEnoughArgument const& argument) {
  if constexpr (low > degree) {
    return Derivative<Value, Argument, low>{};
  } else {
    auto const& x = argument;
    auto const a =
        InternalHornerEvaluator<Value, Argument, degree, fma, low + 1>::
            Evaluate(coefficients, argument);
    auto const& b = get<low>(coefficients);
    if constexpr (fma) {
      return FusedMultiplyAdd(a, x, b);
    } else {
      return a * x + b;
    }
  }
}

template<typename Value, typename Argument, int degree, bool fma, int low>
FORCE_INLINE Derivative<Value, Argument, low>
InternalHornerEvaluator<Value, Argument, degree, fma, low>::
EvaluateDerivative(Coefficients const& coefficients,
                   WideEnoughArgument const& argument) {
  if constexpr (low > degree) {
    return Derivative<Value, Argument, low>{};
  } else {
    auto const& x = argument;
    auto const a =
        InternalHornerEvaluator<Value, Argument, degree, fma, low + 1>::
            EvaluateDerivative(coefficients, argument);
    auto const b = MUL(low, get<low>(coefficients));
    if constexpr (fma) {
      return FusedMultiplyAdd(a, x, b);
    } else {
      return a * x + b;
    }
  }
}

template<typename Value, typename Argument, int degree, bool fma>
FORCE_INLINE Derivative<Value, Argument, degree>
InternalHornerEvaluator<Value, Argument, degree, fma, degree>::
Evaluate(Coefficients const& coefficients,
         WideEnoughArgument const& /*argument*/) {
  return get<degree>(coefficients);
}

template<typename Value, typename Argument, int degree, bool fma>
FORCE_INLINE Derivative<Value, Argument, degree>
InternalHornerEvaluator<Value, Argument, degree, fma, degree>::
EvaluateDerivative(Coefficients const& coefficients,
                   WideEnoughArgument const& /*argument*/) {
  return get<degree>(coefficients) * degree;
}

template<typename Value, typename Argument, int degree,
         FMAPolicy fma_policy, FMAPresence fma_presence>
Value EstrinEvaluator<Value, Argument, degree, fma_policy, fma_presence>::
Evaluate(Coefficients const& coefficients,
         Argument const& argument) {
  auto const wide_enough_argument = WideEnough<Value, Argument>(argument);
  if constexpr (fma_policy == FMAPolicy::Auto) {
    if ((fma_presence == FMAPresence::Present ||
         (fma_presence == FMAPresence::Unknown && CanUseHardwareFMA))) {
      using InternalEvaluator = InternalEstrinEvaluator<Value,
                                                        Argument,
                                                        Coefficients,
                                                        /*fma=*/true,
                                                        /*low=*/0,
                                                        /*degree=*/degree>;
      return InternalEvaluator::Evaluate(coefficients, wide_enough_argument);
    }
  }
  using InternalEvaluator = InternalEstrinEvaluator<Value,
                                                    Argument,
                                                    Coefficients,
                                                    /*fma=*/false,
                                                    /*low=*/0,
                                                    /*degree=*/degree>;
  return InternalEvaluator::Evaluate(coefficients, wide_enough_argument);
}

template<typename Value, typename Argument, int degree,
         FMAPolicy fma_policy, FMAPresence fma_presence>
Derivative<Value, Argument>
EstrinEvaluator<Value, Argument, degree, fma_policy, fma_presence>::
EvaluateDerivative(Coefficients const& coefficients,
                   Argument const& argument) {
  auto const wide_enough_argument = WideEnough<Value, Argument>(argument);
  if constexpr (fma_policy == FMAPolicy::Auto) {
    if (fma_presence == FMAPresence::Present ||
        (fma_presence == FMAPresence::Unknown && CanUseHardwareFMA)) {
      using InternalEvaluator =
          InternalEstrinEvaluator<Value,
                                  Argument,
                                  Coefficients,
                                  /*fma=*/true,
                                  /*low=*/1,
                                  /*degree=*/degree - 1>;
      return InternalEvaluator::EvaluateDerivative(coefficients,
                                                   wide_enough_argument);
    }
  }
  using InternalEvaluator = InternalEstrinEvaluator<Value,
                                                    Argument,
                                                    Coefficients,
                                                    /*fma=*/false,
                                                    /*low=*/1,
                                                    /*degree=*/degree - 1>;
  return InternalEvaluator::EvaluateDerivative(coefficients,
                                               wide_enough_argument);
}

template<typename Value, typename Argument, int degree,
         FMAPolicy fma_policy, FMAPresence fma_presence>
void EstrinEvaluator<Value, Argument, degree, fma_policy, fma_presence>::
EvaluateWithDerivative(Coefficients const& coefficients,
                       Argument const& argument,
                       Value& value,
                       Derivative<Value, Argument>& derivative) {
  if constexpr (degree == 0) {
    value = get<0>(coefficients);
    derivative = Derivative<Value, Argument>{};
    return;
  } else {
    auto const wide_enough_argument = WideEnough<Value, Argument>(argument);
    if constexpr (fma_policy == FMAPolicy::Auto) {
      if (fma_presence == FMAPresence::Present ||
          (fma_presence == FMAPresence::Unknown && CanUseHardwareFMA)) {
        auto const [v, d] = InternalEstrinEvaluator2<Value,
                                                     Argument,
                                                     Coefficients,
                                                     /*fma=*/true,
                                                     /*low=*/0,
                                                     /*degree=*/degree,
                                                     /*lowʹ=*/1,
                                                     /*degreeʹ=*/degree - 1>::
            EvaluateWithDerivative(coefficients, wide_enough_argument);
        value = v;
        derivative = d;
        return;
      }
    }
    auto const [v, d] = InternalEstrinEvaluator2<
        Value,
        Argument,
        Coefficients,
        /*fma=*/false,
        /*low=*/0,
        /*degree=*/degree,
        /*lowʹ=*/1,
        /*degreeʹ=*/degree - 1>::EvaluateWithDerivative(coefficients,
                                                        wide_enough_argument);
    value = v;
    derivative = d;
  }
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
  auto const wide_enough_argument = WideEnough<Value, Argument>(argument);
  if constexpr (fma_policy == FMAPolicy::Auto) {
    if (fma_presence == FMAPresence::Present ||
        (fma_presence == FMAPresence::Unknown && CanUseHardwareFMA)) {
      return InternalHornerEvaluator<Value,
                                     Argument,
                                     degree,
                                     /*fma=*/true,
                                     /*low=*/0>::Evaluate(coefficients,
                                                          wide_enough_argument);
    }
  }
  return InternalHornerEvaluator<Value,
                                 Argument,
                                 degree,
                                 /*fma=*/false,
                                 /*low=*/0>::Evaluate(coefficients,
                                                      wide_enough_argument);
}

template<typename Value, typename Argument, int degree,
         FMAPolicy fma_policy, FMAPresence fma_presence>
Derivative<Value, Argument>
HornerEvaluator<Value, Argument, degree, fma_policy, fma_presence>::
    EvaluateDerivative(Coefficients const& coefficients,
                       Argument const& argument) {
  auto const wide_enough_argument = WideEnough<Value, Argument>(argument);
  if constexpr (fma_policy == FMAPolicy::Auto) {
    if (fma_presence == FMAPresence::Present ||
        (fma_presence == FMAPresence::Unknown && CanUseHardwareFMA)) {
      return InternalHornerEvaluator<
          Value,
          Argument,
          degree,
          /*fma=*/true,
          /*low=*/1>::EvaluateDerivative(coefficients, wide_enough_argument);
    }
  }
  return InternalHornerEvaluator<
      Value,
      Argument,
      degree,
      /*fma=*/false,
      /*low=*/1>::EvaluateDerivative(coefficients, wide_enough_argument);
}

template<typename Value, typename Argument, int degree,
         FMAPolicy fma_policy, FMAPresence fma_presence>
void HornerEvaluator<Value, Argument, degree, fma_policy, fma_presence>::
EvaluateWithDerivative(Coefficients const& coefficients,
                       Argument const& argument,
                       Value& value,
                           Derivative<Value, Argument>& derivative) {
  auto const wide_enough_argument = WideEnough<Value, Argument>(argument);
  if constexpr (fma_policy == FMAPolicy::Auto) {
    if (fma_presence == FMAPresence::Present ||
        (fma_presence == FMAPresence::Unknown && CanUseHardwareFMA)) {
      value =
          InternalHornerEvaluator<Value,
                                  Argument,
                                  degree,
                                  /*fma=*/true,
                                  /*low=*/0>::Evaluate(coefficients,
                                                       wide_enough_argument);
      derivative = InternalHornerEvaluator<
          Value,
          Argument,
          degree,
          /*fma=*/true,
          /*low=*/1>::EvaluateDerivative(coefficients, wide_enough_argument);
    }
  }
  value = InternalHornerEvaluator<Value,
                                  Argument,
                                  degree,
                                  /*fma=*/false,
                                  /*low=*/0>::Evaluate(coefficients,
                                                       wide_enough_argument);
  derivative = InternalHornerEvaluator<
      Value,
      Argument,
      degree,
      /*fma=*/false,
      /*low=*/1>::EvaluateDerivative(coefficients, wide_enough_argument);
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
