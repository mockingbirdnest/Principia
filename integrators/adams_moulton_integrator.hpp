#pragma once

#include "numerics/fixed_arrays.hpp"

namespace principia {

using numerics::FixedVector;

namespace integrators {

// Definition of an Adams-Moulton method.  Not technically an "integrator" as
// it is not a subclass of |integrators::Integrator|.
template<int order_>
struct AdamsMoulton {
  static constexpr int order = order_;
  FixedVector<double, order> numerators;
  double denominator;
};

template<int order>
AdamsMoulton<order> const& AdamsMoultonOrder();

}  // namespace integrators
}  // namespace principia

#include "integrators/adams_moulton_integrator_body.hpp"
