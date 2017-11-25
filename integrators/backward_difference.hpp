#pragma once

#include "numerics/fixed_arrays.hpp"

namespace principia {
namespace integrators {
namespace internal_backward_difference {

using numerics::FixedVector;

// Definition of a backward difference formula.  |order_| is the order of the
// approximation.  The formula requires |order_ + 1| points.  numerators[0]
// corresponds to f(x₀), numerators[1] to f(x₋₁), etc.
template<int order_>
struct BackwardDifference final {
  static constexpr int order = order_;
  FixedVector<double, order + 1> numerators;
  double denominator;
};

}  // namespace internal_backward_difference

using internal_backward_difference::BackwardDifference;

template<int order>
BackwardDifference<order> const& FirstDerivativeBackwardDifference();

}  // namespace integrators
}  // namespace principia

#include "integrators/backward_difference_body.hpp"
