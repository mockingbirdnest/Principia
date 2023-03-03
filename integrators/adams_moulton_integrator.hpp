#pragma once

#include "numerics/fixed_arrays.hpp"

namespace principia {
namespace integrators {
namespace internal_adams_moulton_integrator {

using namespace principia::numerics::_fixed_arrays;

// Definition of an Adams-Moulton method.  Not technically an "integrator" as
// it is not a subclass of |integrators::Integrator|.
template<int order_>
struct AdamsMoulton final {
  static constexpr int order = order_;
  FixedVector<double, order> numerators;
  double denominator;
};

}  // namespace internal_adams_moulton_integrator

using internal_adams_moulton_integrator::AdamsMoulton;

template<int order>
AdamsMoulton<order> const& AdamsMoultonOrder();

}  // namespace integrators
}  // namespace principia

#include "integrators/adams_moulton_integrator_body.hpp"
