
#include <random>

#include "quantities/si.hpp"
#include "numerics/чебышёв_series.hpp"

// Must come last to avoid conflicts when defining the CHECK macros.
#include "benchmark/benchmark.h"

namespace principia {

using numerics::ЧебышёвSeries;
using si::Second;

namespace benchmarks {

void BM_EvaluateDouble(benchmark::State& state) {
  state.PauseTiming();
  int const degree = state.range_x();
  std::mt19937_64 random(42);
  std::vector<double> coefficients;
  for (int i = 0; i <= degree; ++i) {
    coefficients.push_back(static_cast<double>(random()));
  }
  Instant const t_min(random() * Second);
  Instant const t_max = t_min + random() * Second;
  ЧебышёвSeries<double> const series(coefficients, t_min, t_max);

  Instant t = t_min;
  Time const ∆t = (t_max - t_min) * 1E-9;
  double result = 0.0;

  state.ResumeTiming();
  while (state.KeepRunning()) {
    result += series.Evaluate(t);
    t += ∆t;
  }
  state.PauseTiming();

  // This weird call to |SetLabel| has no effect except that it uses |result|
  // and therefore prevents the loop from being optimized away.
  state.SetLabel(std::to_string(result).substr(0, 0));
}

BENCHMARK(BM_EvaluateDouble)->Arg(4)->Arg(8)->Arg(16);

}  // namespace benchmarks
}  // namespace principia
