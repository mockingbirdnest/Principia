
#pragma once

#include "numerics/legendre.hpp"

#include <tuple>

#include "base/macros.hpp"
#include "numerics/combinatorics.hpp"
#include "quantities/elementary_functions.hpp"

namespace principia {
namespace numerics {
namespace internal_legendre {

using quantities::Sqrt;

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
template<int degree_, template<typename, typename, int> class Evaluator>
FORCE_INLINE(constexpr)
PolynomialInMonomialBasis<double, double, degree_, Evaluator>
LegendrePolynomial() {
  auto c = LegendreCoefficientsGenerator<
          degree_,
          std::make_index_sequence<degree_ + 1>>::coefficients;
  return PolynomialInMonomialBasis<double, double, degree_, Evaluator>(
      LegendreCoefficientsGenerator<
          degree_,
          std::make_index_sequence<degree_ + 1>>::coefficients);
}

}  // namespace internal_legendre
}  // namespace numerics
}  // namespace principia
