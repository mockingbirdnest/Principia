#include "clr_benchmarks_adapter/symplectic_partitioned_runge_kutta_integrator.hpp"  // NOLINT(whitespace/line_length)

#include "benchmarks/symplectic_partitioned_runge_kutta_integrator.hpp"
#undef TRACE_SYMPLECTIC_PARTITIONED_RUNGE_KUTTA_INTEGRATOR
#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"

using principia::integrators::SPRKIntegrator;

namespace principia {
namespace clr_benchmarks_adapter {

void SPRKIntegratorCLRBenchmark::SolveHarmonicOscillator() {
  SPRKIntegrator::Solution solution;
  principia::benchmarks::SolveHarmonicOscillator(&solution);
}

}  // namespace clr_benchmarks_adapter
}  // namespace principia
