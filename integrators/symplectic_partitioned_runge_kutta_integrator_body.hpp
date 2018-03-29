
#pragma once

#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"

#include "base/mod.hpp"

namespace principia {
namespace integrators {
namespace internal_symplectic_runge_kutta_nyström_integrator {

using base::mod;

template<typename Method, typename Position>
SymplecticPartitionedRungeKuttaIntegrator<Method, Position>::
    SymplecticPartitionedRungeKuttaIntegrator() {
  // TODO(phl): This might be turned into a static_assert.
  if (Method::first_same_as_last) {
    CHECK_EQ(0.0, Method::a[Method::stages - 1]);
  }
  if (Method::time_reversible) {
    CHECK(Method::first_same_as_last);
    for (int i = 0; i < Method::stages - 1; ++i) {
      CHECK_EQ(Method::a[i], Method::a[Method::stages - 2 - i]);
    }
    for (int i = 0; i < Method::stages; ++i) {
      CHECK_EQ(Method::b[i], Method::b[stages - 1 - i]);
    }
  }
}

}  // namespace internal_symplectic_runge_kutta_nyström_integrator


}  // namespace integrators
}  // namespace principia
