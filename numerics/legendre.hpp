
#pragma once

#include "numerics/polynomial.hpp"

namespace principia {
namespace numerics {
namespace internal_legendre {

template<int degree_, template<typename, typename, int> class Evaluator>
class LegendrePolynomial : public 
  PolynomialInMonomialBasis<double, double, degree_, Evaluator> {
 public:
  LegendrePolynomial();
};

}  // namespace internal_legendre
}  // namespace numerics
}  // namespace principia

#include "numerics/legendre_body.hpp"
