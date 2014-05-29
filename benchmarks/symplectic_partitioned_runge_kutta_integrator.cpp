#define GLOG_NO_ABBREVIATED_SEVERITIES

#include <algorithm>

#include "benchmark/benchmark.h"
#include "integrators/symplectic_partitioned_runge_kutta_integrator.hpp"

using principia::integrators::SPRKIntegrator;

namespace principia {
namespace benchmarks {

namespace {

inline void compute_harmonic_oscillator_force(double const t,
                                              std::vector<double> const& q,
                                              std::vector<double>* result) {
  (*result)[0] = -q[0];
}

inline void compute_harmonice_oscillator_velocity(std::vector<double> const& p,
                                                  std::vector<double>* result) {
  (*result)[0] = p[0];
}

}  // namespace

void SolveHarmonicOscillator() {
  std::unique_ptr<SPRKIntegrator> integrator(
      new SPRKIntegrator);
  std::unique_ptr<SPRKIntegrator::Parameters> parameters(
      new SPRKIntegrator::Parameters);
  std::unique_ptr<SPRKIntegrator::Solution> solution(
      new SPRKIntegrator::Solution);

  parameters->q0 = {1.0};
  parameters->p0 = {0.0};
  parameters->t0 = 0.0;
#ifdef _DEBUG
  parameters->tmax = 100.0;
#else
  parameters->tmax = 1000.0;
#endif
  parameters->Δt = 1.0E-4;
  parameters->coefficients = integrator->Order5Optimal();
  parameters->sampling_period = 1;
  integrator->Solve(&compute_harmonic_oscillator_force,
                         &compute_harmonice_oscillator_velocity,
                         *parameters,
                         solution.get());
  double q_error = 0;
  double p_error = 0;
  for (size_t i = 0; i < solution->time.quantities.size(); ++i) {
    q_error = std::max(q_error,
                       std::abs(solution->position[0].quantities[i] -
                                std::cos(solution->time.quantities[i])));
    p_error = std::max(p_error,
                       std::abs(solution->momentum[0].quantities[i] +
                                std::sin(solution->time.quantities[i])));
  }
  //LOG(ERROR) << "q_error = " << q_error;
  //LOG(ERROR) << "p_error = " << p_error;
  //EXPECT_THAT(AbsoluteError(0, q_error), Lt(2E-16 * parameters_->tmax));
  //EXPECT_THAT(AbsoluteError(0, p_error), Lt(2E-16 * parameters_->tmax));
}

static void BM_SolveHarmonicOscillator(benchmark::State& state) {
  while (state.KeepRunning()) {
    SolveHarmonicOscillator();
  }
}
BENCHMARK(BM_SolveHarmonicOscillator);

}  // namespace benchmarks
}  // namespace principia

int main(int argc, const char* argv[]) {
  benchmark::Initialize(&argc, argv);

  benchmark::RunSpecifiedBenchmarks();
}
