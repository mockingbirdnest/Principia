#pragma once

#include "integrators/explicit_embedded_runge_kutta_nyström_integrator.hpp"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <vector>

#include "glog/logging.h"
#include "quantities/quantities.hpp"

namespace principia {

using quantities::Difference;
using quantities::Quotient;

namespace integrators {

template<typename Position>
void ExplicitEmbeddedRungeKuttaNyströmIntegrator::Solve(
    RightHandSideComputation<Position> compute_acceleration,
    SystemState<Position, Variation<Position>> const& initial_value,
    Time const& t_final,
    Time const& first_time_step,
    Difference<Position> const& position_tolerance,
    Variation<Position> const& velocity_tolerance,
    double const safety_factor,
    not_null<Solution<Position, Variation<Position>>*> const solution) const {
  for (int i = 0; i < stages_; ++i) {
    for (int j = 0; j < stages_; ++j) {
    }
  }
}

}  // namespace integrators
}  // namespace principia
