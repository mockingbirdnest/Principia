#pragma once

#include <vector>

#undef TRACE_SYMPLECTIC_PARTITIONED_RUNGE_KUTTA_INTEGRATOR

#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"
#include "testing_utilities/numerical_analysis.hpp"

using principia::integrators::SPRKIntegrator;

namespace principia {
namespace benchmarks {

inline void SolveHarmonicOscillator(SPRKIntegrator::Solution* solution) {
  using principia::testing_utilities::ComputeHarmonicOscillatorForce;
  using principia::testing_utilities::ComputeHarmonicOscillatorVelocity;
  SPRKIntegrator integrator;
  SPRKIntegrator::Parameters parameters;

  integrator.Initialize(integrator.Order5Optimal());

  parameters.q0 = {1.0};
  parameters.p0 = {0.0};
  parameters.t0 = 0.0;
#ifdef _DEBUG
  parameters.tmax = 100.0;
#else
  parameters.tmax = 1000.0;
#endif
  parameters.Δt = 1.0E-4;
  parameters.sampling_period = 1;
  integrator.Solve(&ComputeHarmonicOscillatorForce,
                   &ComputeHarmonicOscillatorVelocity,
                   parameters,
                   solution);
}

}  // namespace benchmarks
}  // namespace principia
