#pragma once

#include "numerics/polynomial_in_monomial_basis.hpp"

namespace principia {
namespace numerics {
namespace _legendre {
namespace internal {

using namespace principia::numerics::_polynomial_in_monomial_basis;

template<int degree_>
constexpr PolynomialInMonomialBasis<double, double, degree_>
LegendrePolynomial();

}  // namespace internal

using internal::LegendrePolynomial;

}  // namespace _legendre
}  // namespace numerics
}  // namespace principia

#include "numerics/legendre_body.hpp"
