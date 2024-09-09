#pragma once

#include "numerics/fixed_arrays.hpp"

namespace principia {
namespace integrators {
namespace _adams_moulton_integrator {
namespace internal {

using namespace principia::numerics::_fixed_arrays;

// Definition of an Adams-Moulton method.  Not technically an "integrator" as
// it is not a subclass of `integrators::Integrator`.
template<int order_>
struct AdamsMoulton final {
  static constexpr int order = order_;
  FixedVector<double, order> numerators;
  double denominator;
};

}  // namespace internal

using internal::AdamsMoulton;

template<int order>
AdamsMoulton<order> const& AdamsMoultonOrder();

}  // namespace _adams_moulton_integrator
}  // namespace integrators
}  // namespace principia

#include "integrators/adams_moulton_integrator_body.hpp"
