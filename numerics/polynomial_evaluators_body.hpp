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
// You'd think that the compiler would optimize this but it doesn't.
// I am leaving a `static_assert` because we'll need to recheck (using the
// nanobenchmarks) that this function is still beneficial after any change of
// compiler.
#if PRINCIPIA_COMPILER_MSVC
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
  // There is a wide flexibility in the choice of `m`, which is not documented
  // in the literature.  [Mul16], p. 47 gives the example of degree 7 and says
  // that the algorithm "can be easily extended to ... any degree".  [Est60]
  // takes `2ᵐ < degree` which is a reasonable choice, but doesn't justify it.
  // https://en.wikipedia.org/wiki/Estrin%27s_scheme takes `2ᵐ ≤ degree` which
  // is a lame choice when `degree` is a power of 2 because the first recursion
  // has `a` of degree 0 and doesn't usefully split the polynomial.
  // Some choices of `m` are unambiguously bad.  For instance, `m = 2` doesn't
  // requires to evaluate powers above `x²`, but it sequentializes the FMAs
  // above the innermost ones.  We want to choose `m` as small as possible (to
  // avoid computing large powers of `x`), but large enough to split the
  // polynomial into two pieces of roughly equal degree, to preserve the
  // independence of the FMAs.  The choice below achieves this.
  constexpr int m = ClosestPowerOf2(degree) / 2;

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
  constexpr int m = ClosestPowerOf2(degree) / 2;

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
  constexpr int m = ClosestPowerOf2(degree) / 2;
  // In order to comply with the assertions above, `mʹ` must be either `m` or
  // `m - 1`.  Choosing it to be `degreeʹ < degree && m > 1 ? m - 1 : m`
  // would result in a better "alignment" of the computations (`lowʹ == low`
  // most of the time, so the end of the recursion only accesses two
  // coefficients instead of three).  Unfortunately, this would require the
  // extra computation of `Pow<mʹ>(x)`.  Since the compiler doesn't take
  // advantage of the improved "alignment", we'll go with the simplest choice.
  constexpr int mʹ = m;

  auto const a = InternalEstrinEvaluator2<
                     Value, Argument, Coefficients,
                     fma, low + m, degree - m, lowʹ + mʹ, degreeʹ - mʹ>::
                 EvaluateWithDerivative(coefficients, argument);
  auto const b = InternalEstrinEvaluator2<
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
#if 0
  } else if constexpr (degree == 17) {
    auto const wide_enough_argument = WideEnough<Value, Argument>(argument);
    auto const& x = wide_enough_argument;

    auto const aaaaa = get<17>(coefficients);
    auto const aaaab = get<16>(coefficients);
    auto const aaaabʹ = MUL(17, get<17>(coefficients));
    auto const aaaa = FusedMultiplyAdd(aaaaa, x, aaaab);
    auto const aaaaʹ = aaaabʹ;

    auto const aaaba = get<15>(coefficients);
    auto const aaabb = get<14>(coefficients);
    auto const aaabaʹ = MUL(16, get<16>(coefficients));
    auto const aaabbʹ = MUL(15, get<15>(coefficients));
    auto const aaab = FusedMultiplyAdd(aaaba, x, aaabb);
    auto const aaabʹ = FusedMultiplyAdd(aaabaʹ, x, aaabbʹ);

    auto const xaaa = Pow<2>(x);
    auto const aaa = FusedMultiplyAdd(aaaa, xaaa, aaab);
    auto const aaaʹ = FusedMultiplyAdd(aaaaʹ, xaaa, aaabʹ);

    auto const aaba = get<13>(coefficients);
    auto const aabb = get<12>(coefficients);
    auto const aabaʹ = MUL(14, get<14>(coefficients));
    auto const aabbʹ = MUL(13, get<13>(coefficients));
    auto const aab = FusedMultiplyAdd(aaba, x, aabb);
    auto const aabʹ = FusedMultiplyAdd(aabaʹ, x, aabbʹ);

    auto const xaa = Pow<2>(x);
    auto const aa = FusedMultiplyAdd(aaa, xaa, aab);
    auto const aaʹ = FusedMultiplyAdd(aaaʹ, xaa, aabʹ);

    auto const abaa = get<11>(coefficients);
    auto const abab = get<10>(coefficients);
    auto const abaaʹ = MUL(12, get<12>(coefficients));
    auto const ababʹ = MUL(11, get<11>(coefficients));
    auto const aba = FusedMultiplyAdd(abaa, x, abab);
    auto const abaʹ = FusedMultiplyAdd(abaaʹ, x, ababʹ);

    auto const abba = get<9>(coefficients);
    auto const abbb = get<8>(coefficients);
    auto const abbaʹ = MUL(10, get<10>(coefficients));
    auto const abbbʹ = MUL(9, get<9>(coefficients));
    auto const abb = FusedMultiplyAdd(abba, x, abbb);
    auto const abbʹ = FusedMultiplyAdd(abbaʹ, x, abbbʹ);

    auto const xab = Pow<2>(x);
    auto const ab = FusedMultiplyAdd(aba, xab, abb);
    auto const abʹ = FusedMultiplyAdd(abaʹ, xab, abbʹ);

    auto const xa = Pow<4>(x);
    auto const a = FusedMultiplyAdd(aa, xa, ab);
    auto const aʹ = FusedMultiplyAdd(aaʹ, xa, abʹ);

    auto const baaa = get<7>(coefficients);
    auto const baab = get<6>(coefficients);
    auto const baaaʹ = MUL(8, get<8>(coefficients));
    auto const baabʹ = MUL(7, get<7>(coefficients));
    auto const baa = FusedMultiplyAdd(baaa, x, baab);
    auto const baaʹ = FusedMultiplyAdd(baaaʹ, x, baabʹ);

    auto const baba = get<5>(coefficients);
    auto const babb = get<4>(coefficients);
    auto const babaʹ = MUL(6, get<6>(coefficients));
    auto const babbʹ = MUL(5, get<5>(coefficients));
    auto const bab = FusedMultiplyAdd(baba, x, babb);
    auto const babʹ = FusedMultiplyAdd(babaʹ, x, babbʹ);

    auto const xba = Pow<2>(x);
    auto const ba = FusedMultiplyAdd(baa, xba, bab);
    auto const baʹ = FusedMultiplyAdd(baaʹ, xba, babʹ);

    auto const bbaa = get<3>(coefficients);
    auto const bbab = get<2>(coefficients);
    auto const bbaaʹ = MUL(4, get<4>(coefficients));
    auto const bbabʹ = MUL(3, get<3>(coefficients));
    auto const bba = FusedMultiplyAdd(bbaa, x, bbab);
    auto const bbaʹ = FusedMultiplyAdd(bbaaʹ, x, bbabʹ);

    auto const bbba = get<1>(coefficients);
    auto const bbbb = get<0>(coefficients);
    auto const bbbaʹ = MUL(2, get<2>(coefficients));
    auto const bbbbʹ = MUL(1, get<1>(coefficients));
    auto const bbb = FusedMultiplyAdd(bbba, x, bbbb);
    auto const bbbʹ = FusedMultiplyAdd(bbbaʹ, x, bbbbʹ);

    auto const xbb = Pow<2>(x);
    auto const bb = FusedMultiplyAdd(bba, xbb, bbb);
    auto const bbʹ = FusedMultiplyAdd(bbaʹ, xbb, bbbʹ);

    auto const xb = Pow<4>(x);
    auto const b = FusedMultiplyAdd(ba, xb, bb);
    auto const bʹ = FusedMultiplyAdd(baʹ, xb, bbʹ);

    auto const x8 = Pow<8>(x);
    value = FusedMultiplyAdd(a, x8, b);
    derivative = FusedMultiplyAdd(aʹ, x8, bʹ);
#endif
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
