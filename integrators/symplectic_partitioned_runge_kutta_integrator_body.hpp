
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
  if (first_same_as_last) {
    CHECK_EQ(0.0, a_[stages_ - 1]);
  }
  if (time_reversible) {
    CHECK(first_same_as_last);
    for (int i = 0; i < stages_ - 1; ++i) {
      CHECK_EQ(a_[i], a_[stages_ - 2 - i]);
    }
    for (int i = 0; i < stages_; ++i) {
      CHECK_EQ(b_[i], b_[stages_ - 1 - i]);
    }
  }
}

}  // namespace internal_symplectic_runge_kutta_nyström_integrator


}  // namespace integrators
}  // namespace principia
