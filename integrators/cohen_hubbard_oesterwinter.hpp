
#pragma once

#include "numerics/fixed_arrays.hpp"

namespace principia {
namespace integrators {
namespace internal_cohen_hubbard_oesterwinter {

using numerics::FixedVector;

// Definition of a Cohen-Hubbard-Oesterwinter formula.  |order_| is the order of
// the approximation, that is, the error on the derivative is h^(order + 1).
// The formula requires |order_| values of the acceleration.  numerators[0]
// corresponds to f"(x₀), numerators[1] to f"(x₋₁), etc.
// TODO(phl): This struct shows up in many places.  Factor it out.
template<int order_>
struct CohenHubbardOesterwinter final {
  static constexpr int order = order_;
  FixedVector<double, order> numerators;
  double denominator;
};

}  // namespace internal_cohen_hubbard_oesterwinter

using internal_cohen_hubbard_oesterwinter::CohenHubbardOesterwinter;

template<int order>
CohenHubbardOesterwinter<order> const& CohenHubbardOesterwinterOrder();

}  // namespace integrators
}  // namespace principia

#include "integrators/cohen_hubbard_oesterwinter_body.hpp"
