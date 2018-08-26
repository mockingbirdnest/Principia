
#pragma once

#include "base/macros.hpp"
#include "numerics/polynomial.hpp"

namespace principia {
namespace numerics {
namespace internal_legendre {

template<int degree_, template<typename, typename, int> class Evaluator>
FORCE_INLINE(constexpr)
PolynomialInMonomialBasis<double, double, degree_, Evaluator>
LegendrePolynomial();

}  // namespace internal_legendre
}  // namespace numerics
}  // namespace principia

#include "numerics/legendre_body.hpp"
