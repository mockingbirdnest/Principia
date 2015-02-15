#pragma once

#undef TRACE_SYMPLECTIC_PARTITIONED_RUNGE_KUTTA_INTEGRATOR

#include <vector>

#include "base/not_null.hpp"
#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"
#include "quantities/named_quantities.hpp"

namespace principia {

using base::not_null;
using integrators::SPRKIntegrator;
using quantities::Length;
using quantities::Momentum;

namespace benchmarks {

inline void SolveHarmonicOscillator(
    not_null<std::vector<
        SPRKIntegrator<Length, Momentum>::SystemState>*> const solution);

}  // namespace benchmarks
}  // namespace principia

#include "benchmarks/symplectic_partitioned_runge_kutta_integrator_body.hpp"
