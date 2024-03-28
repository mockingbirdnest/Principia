#pragma once

#include "numerics/legendre.hpp"

#include <tuple>

#include "numerics/combinatorics.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace _legendre {
namespace internal {

using namespace principia::numerics::_combinatorics;
using namespace principia::quantities::_elementary_functions;

template<int degree, typename>
struct LegendreCoefficientsGenerator;

template<int degree, std::size_t... indices>
struct LegendreCoefficientsGenerator<degree, std::index_sequence<indices...>> {
  // This computation follows
  // https://en.wikipedia.org/wiki/Legendre_polynomials, fourth formula in the
  // "Explicit representations" section.  The formula has been rewritten to
  // eliminate references to the Γ function.
  static constexpr auto coefficients = std::make_tuple(
      (degree - indices) % 2 == 0
          ? ((degree - indices) % 4 == 0 ? 1 : -1) *
                DoubleFactorial(degree + indices - 1) /
                static_cast<double>(Factorial(indices) *
                                    DoubleFactorial(degree - indices))
          : 0 ...);
};

// Apparently, FORCE_INLINE has to be on the definition for it to work on
// namespace-level functions.
template<int degree_>
FORCE_INLINE(constexpr)
PolynomialInMonomialBasis<double, double, degree_>
LegendrePolynomial() {
  return PolynomialInMonomialBasis<double, double, degree_>(
      LegendreCoefficientsGenerator<
          degree_,
          std::make_index_sequence<degree_ + 1>>::coefficients);
}

}  // namespace internal
}  // namespace _legendre
}  // namespace numerics
}  // namespace principia
