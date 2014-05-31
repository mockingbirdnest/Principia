#include "clr_benchmarks_adapter/symplectic_partitioned_runge_kutta_integrator.hpp"

#undef TRACE_SYMPLECTIC_PARTITIONED_RUNGE_KUTTA_INTEGRATOR
#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"

using principia::integrators::SPRKIntegrator;

namespace principia {
namespace clr_benchmarks {

namespace {

inline void compute_harmonic_oscillator_force(double const t,
                                              std::vector<double> const& q,
                                              std::vector<double>* result) {
  (*result)[0] = -q[0];
}

inline void compute_harmonic_oscillator_velocity(std::vector<double> const& p,
                                                 std::vector<double>* result) {
  (*result)[0] = p[0];
}

}  // namespace

void SPRKIntegratorCLRBenchmark::SolveHarmonicOscillator() {
  SPRKIntegrator integrator;
  SPRKIntegrator::Parameters parameters;
  SPRKIntegrator::Solution solution;

  parameters.q0 = {1.0};
  parameters.p0 = {0.0};
  parameters.t0 = 0.0;
#ifdef _DEBUG
  parameters.tmax = 100.0;
#else
  parameters.tmax = 1000.0;
#endif
  parameters.Δt = 1.0E-4;
  parameters.coefficients = integrator.Order5Optimal();
  parameters.sampling_period = 1;
  integrator.Solve(&compute_harmonic_oscillator_force,
                   &compute_harmonic_oscillator_velocity,
                   parameters,
                   &solution);
}

}  // namespace clr_benchmarks
}  // namespace principia
