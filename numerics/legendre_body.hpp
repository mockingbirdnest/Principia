
#pragma once

#include "numerics/legendre.hpp"

#include "base/macros.hpp"

namespace principia {
namespace numerics {
namespace internal_legendre {

// Apparently, FORCE_INLINE has to be on the definition for it to work on
// namespace-level functions.
template<int degree_, template<typename, typename, int> class Evaluator>
FORCE_INLINE(constexpr)
PolynomialInMonomialBasis<double, double, degree_, Evaluator>
LegendrePolynomial() {
  using Pn = PolynomialInMonomialBasis<double, double, degree_, Evaluator>;
  if constexpr (degree_ == 0) {
    return Pn(std::make_tuple(1));
  } else if constexpr (degree_ == 1) {
    return Pn({0, 1});
  } else {
    constexpr int n = degree_;
    // NOTE(phl): Would like to make this declaration constexpr, but the type is
    // not a literal type.
    PolynomialInMonomialBasis<double, double, 1, Evaluator> const
        multiplier({0, 2 * n - 1});
    return (multiplier * LegendrePolynomial<degree_ - 1, Evaluator>() -
            (n - 1) * LegendrePolynomial<degree_ - 2, Evaluator>()) / n;
  }
}

}  // namespace internal_legendre
}  // namespace numerics
}  // namespace principia
