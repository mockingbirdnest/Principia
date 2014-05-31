#pragma once

#undef TRACE_SYMPLECTIC_PARTITIONED_RUNGE_KUTTA_INTEGRATOR

#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"

using principia::integrators::SPRKIntegrator;

namespace principia {
namespace benchmarks {

inline void SolveHarmonicOscillator(SPRKIntegrator::Solution* solution);

}  // namespace benchmarks
}  // namespace principia

#include "benchmarks/symplectic_partitioned_runge_kutta_integrator_body.hpp"
