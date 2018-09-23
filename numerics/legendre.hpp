
#pragma once

#include "numerics/polynomial.hpp"

namespace principia {
namespace numerics {
namespace internal_legendre {

template<int degree_, template<typename, typename, int> class Evaluator>
constexpr PolynomialInMonomialBasis<double, double, degree_, Evaluator>
LegendrePolynomial();

// Multiplying a normalized Cnm or Snm coefficient by this factor yields an
// unnormalized coefficient.  Dividing an unnormalized Cnm or Snm coefficient by
// this factor yields a normalized coefficient.
double LegendreNormalizationFactor(int degree, int order);

}  // namespace internal_legendre

using internal_legendre::LegendreNormalizationFactor;
using internal_legendre::LegendrePolynomial;

}  // namespace numerics
}  // namespace principia

#include "numerics/legendre_body.hpp"
