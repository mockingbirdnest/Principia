#pragma once

#include <vector>

#undef TRACE_SYMPLECTIC_PARTITIONED_RUNGE_KUTTA_INTEGRATOR

#include "base/not_null.hpp"
#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"
#include "testing_utilities/numerical_analysis.hpp"

namespace principia {

using base::not_null;
using integrators::SPRKIntegrator;
using quantities::Length;
using quantities::Momentum;
using quantities::SIUnit;
using testing_utilities::ComputeHarmonicOscillatorForce;
using testing_utilities::ComputeHarmonicOscillatorVelocity;

namespace benchmarks {

inline void SolveHarmonicOscillator(
    not_null<std::vector<
        SPRKIntegrator<Length, Momentum>::SystemState>*> const solution) {
  SPRKIntegrator<Length, Momentum> integrator;
  SPRKIntegrator<Length, Momentum>::Parameters parameters;

  integrator.Initialize(integrator.McLachlanAtela1992Order5Optimal());

  parameters.initial.positions.emplace_back(SIUnit<Length>());
  parameters.initial.momenta.emplace_back(Momentum());
  parameters.initial.time = Time();
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
