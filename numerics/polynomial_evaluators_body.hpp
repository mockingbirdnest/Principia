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
  Derivative<Value, Argument, lowʹ> derivative;
};

// This function saves a multiplication by 1.0 and a multiplication by 2.0.
// The compiler is apparently unable to optimize these because the code is
// produced by intrinsics in `R3Element::operator*`.
#if PRINCIPIA_COMPILER_MSVC
template<int m, typename Value>
  requires(m > 0)
Value Mul(Value const x) {
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

// Similar to the above helper, but evaluates both the value and the derivative
// of the polynomial.  `lowʹ` and `degreeʹ` define the subpolynomial that we
// currently evaluate for the derivative.
template<typename Value, typename Argument, typename Coefficients,
         bool fma, int low, int degree, int lowʹ, int degreeʹ>
struct InternalEstrinEvaluatorWithDerivative {
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
struct InternalEstrinEvaluatorWithDerivative<Value, Argument, Coefficients,
                                             fma, low, 1, lowʹ, 1> {
  using WideEnoughArgument = WideEnough<Value, Argument>;

  FORCE_INLINE static ValueAndDerivative<Value, Argument, low, lowʹ>
  EvaluateWithDerivative(Coefficients const& coefficients,
                         WideEnoughArgument const& argument);
};

template<typename Value, typename Argument, typename Coefficients,
         bool fma, int low, int lowʹ>
struct InternalEstrinEvaluatorWithDerivative<Value, Argument, Coefficients,
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
struct InternalEstrinEvaluatorWithDerivative<Value, Argument, Coefficients,
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
  // There is a wide flexibility in the choice of `m`, which is not documented
  // in the literature.  [Mul16], p. 47 gives the example of degree 7 and says
  // that the algorithm "can be easily extended to ... any degree".  [Est60]
  // takes `2ᵐ < degree` which is a reasonable choice, but he doesn't justify
  // it.  https://en.wikipedia.org/wiki/Estrin%27s_scheme takes `2ᵐ ≤ degree`
  // which is a lame choice when `degree` is a power of 2 because the first
  // recursion has `a` of degree 0 and doesn't usefully split the polynomial.
  // Some choices of `m` are unambiguously bad.  For instance, `m = 2` doesn't
  // require evaluating powers above `x²`, but it sequentializes the FMAs
  // above the innermost ones.  Here and in `EvaluateDerivative` we use Estrin's
  // choice.
  constexpr int m = PowerOf2Lt(degree);

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
  constexpr int m = PowerOf2Lt(degree);

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
InternalEstrinEvaluatorWithDerivative<Value, Argument, Coefficients,
                                      fma, low, degree, lowʹ, degreeʹ>::
EvaluateWithDerivative(Coefficients const& coefficients,
                       WideEnoughArgument const& argument) {
  static_assert(degree >= 2,
                "Unexpected degree in InternalEstrinEvaluatorWithDerivative::"
                "EvaluateWithDerivative");
  static_assert(degree == degreeʹ || degree == degreeʹ + 1);
  static_assert(low == lowʹ || low == lowʹ - 1);
  // Here we choose `m` as small as possible (to avoid computing large powers of
  // `x`), but large enough to split the polynomial into two pieces of roughly
  // equal degree, to preserve the independence of the FMAs.  The evaluation
  // tree may end up a bit deeper with this choice, and the leaves look more
  // like Horner, but the benchmarks demonstrate that this is preferable to
  // Estrin's choice (because it relieves register pressure).
  constexpr int m = ClosestPowerOf2(degree) / 2;
  // In order to comply with the assertions above, `mʹ` must be either `m` or
  // `m - 1`.  Choosing it to be `degreeʹ < degree && m > 1 ? m - 1 : m`
  // would result in a better "alignment" of the computations (`lowʹ == low`
  // most of the time, so the end of the recursion only accesses two
  // coefficients instead of three).  Unfortunately, this would require the
  // extra computation of `Pow<mʹ>(x)`, which makes things more expensive.
  constexpr int mʹ = m;

  auto const a = InternalEstrinEvaluatorWithDerivative<
                     Value, Argument, Coefficients,
                     fma, low + m, degree - m, lowʹ + mʹ, degreeʹ - mʹ>::
                 EvaluateWithDerivative(coefficients, argument);
  auto const b = InternalEstrinEvaluatorWithDerivative<
                     Value, Argument, Coefficients,
                     fma, low, m - 1, lowʹ, mʹ - 1>::
                 EvaluateWithDerivative(coefficients, argument);

  auto const& x = argument;
  auto const xᵐ = Pow<m>(x);
  auto const xᵐʹ = xᵐ;
  if constexpr (fma) {
    return {.value = FusedMultiplyAdd(a.value, xᵐ, b.value),
            .derivative = FusedMultiplyAdd(a.derivative, xᵐʹ, b.derivative)};
  } else {
    return {.value = a.value * xᵐ + b.value,
            .derivative = a.derivative * xᵐʹ + b.derivative};
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
InternalEstrinEvaluatorWithDerivative<Value, Argument, Coefficients,
                                      fma, low, 1, lowʹ, 1>::
EvaluateWithDerivative(Coefficients const& coefficients,
                       WideEnoughArgument const& argument) {
  // No point in trying to factor out the inner FMA in the computation of the
  // derivative below because the degrees are never aligned.  If this ever
  // changed, this assertion would fire.
  static_assert(low == lowʹ - 1);
  // Don't use references in these declarations, it degrades the performance of
  // the generated code.
  auto const x = argument;
  auto const a = get<low + 1>(coefficients);
  auto const b = get<low>(coefficients);
  auto const aʹ = get<lowʹ + 1>(coefficients);
  auto const bʹ = get<lowʹ>(coefficients);
  if constexpr (fma) {
    // The strange formula for the derivative replaces two multiplications and
    // one FMA by two FMAs.  It also avoids having to reify the value of `lowʹ`
    // in the generated code.  Note that this does not appear to be beneficial
    // when we only compute the derivative.
    return {.value = FusedMultiplyAdd(a, x, b),
            .derivative = FusedMultiplySubtract(
                FusedMultiplyAdd(aʹ, x, bʹ), lowʹ + 1, bʹ)};
  } else {
    return {.value = a * x + b, .derivative = (aʹ * x + bʹ) * (lowʹ + 1) - bʹ};
  }
}

template<typename Value, typename Argument, typename Coefficients,
         bool fma, int low, int lowʹ>
FORCE_INLINE ValueAndDerivative<Value, Argument, low, lowʹ>
InternalEstrinEvaluatorWithDerivative<Value, Argument, Coefficients,
                                      fma, low, 1, lowʹ, 0>::
EvaluateWithDerivative(Coefficients const& coefficients,
                       WideEnoughArgument const& argument) {
  auto const x = argument;
  auto const a = get<low + 1>(coefficients);
  auto const b = get<low>(coefficients);
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
InternalEstrinEvaluatorWithDerivative<Value, Argument, Coefficients,
                                      fma, low, 0, lowʹ, 0>::
EvaluateWithDerivative(Coefficients const& coefficients,
                       WideEnoughArgument const& /*argument*/) {
  auto const b = get<low>(coefficients);
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
  } else {
    auto const wide_enough_argument = WideEnough<Value, Argument>(argument);
    if constexpr (fma_policy == FMAPolicy::Auto) {
      if (fma_presence == FMAPresence::Present ||
          (fma_presence == FMAPresence::Unknown && CanUseHardwareFMA)) {
        auto const [v, d] =
            InternalEstrinEvaluatorWithDerivative<Value,
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
    auto const [v, d] = InternalEstrinEvaluatorWithDerivative<
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
