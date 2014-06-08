#pragma once

#include <vector>

#undef TRACE_SYMPLECTIC_PARTITIONED_RUNGE_KUTTA_INTEGRATOR

#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"
#include "testing_utilities/numerical_analysis.hpp"

using principia::integrators::SPRKIntegrator;
using principia::quantities::Length;
using principia::quantities::Momentum;
using principia::quantities::SIUnit;

namespace principia {
namespace benchmarks {

inline void SolveHarmonicOscillator(
    SPRKIntegrator<Length, Momentum>::Solution* solution) {
  using principia::testing_utilities::ComputeHarmonicOscillatorForce;
  using principia::testing_utilities::ComputeHarmonicOscillatorVelocity;
  SPRKIntegrator<Length, Momentum> integrator;
  SPRKIntegrator<Length, Momentum>::Parameters parameters;

  integrator.Initialize(integrator.Order5Optimal());

  parameters.q0 = {SIUnit<Length>()};
  parameters.p0 = {0.0 * SIUnit<Momentum>()};
  parameters.t0 = 0.0 * SIUnit<Time>();
#ifdef _DEBUG
  parameters.tmax = 100.0 * SIUnit<Time>();
#else
  parameters.tmax = 1000.0 * SIUnit<Time>();
#endif
  parameters.Δt = 1.0E-4 * SIUnit<Time>();
  parameters.sampling_period = 1;
  integrator.Solve(&ComputeHarmonicOscillatorForce,
                   &ComputeHarmonicOscillatorVelocity,
                   parameters,
                   solution);
}

}  // namespace benchmarks
}  // namespace principia
