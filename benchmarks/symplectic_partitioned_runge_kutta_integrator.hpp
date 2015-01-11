#pragma once

#undef TRACE_SYMPLECTIC_PARTITIONED_RUNGE_KUTTA_INTEGRATOR

#include <vector>

#include "base/not_null.hpp"
#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"
#include "quantities/named_quantities.hpp"

using principia::base::not_null;
using principia::integrators::SPRKIntegrator;
using principia::quantities::Length;
using principia::quantities::Momentum;

namespace principia {
namespace benchmarks {

inline void SolveHarmonicOscillator(
    not_null<std::vector<
        SPRKIntegrator<Length, Momentum>::SystemState>*> const solution);

}  // namespace benchmarks
}  // namespace principia

#include "benchmarks/symplectic_partitioned_runge_kutta_integrator_body.hpp"
