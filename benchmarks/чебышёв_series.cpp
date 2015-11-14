
// .\Release\benchmarks.exe --benchmark_repetitions=3 --benchmark_filter=Evaluate  // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 3310 MHz CPU
// 2015/06/09-22:18:17
// Benchmark                    Time(ns)    CPU(ns) Iterations
// -----------------------------------------------------------
// BM_EvaluateDouble/4              8031       8139      63254
// BM_EvaluateDouble/8             15548      15602      32995
// BM_EvaluateDouble/15            35945      36121      14252
// BM_EvaluateDouble/16            36019      36104      14259
// BM_EvaluateDouble/17            39483      39558      13014
// BM_EvaluateDouble/18            42357      42413      12138
// BM_EvaluateDouble/19            45989      44785      11495
// BM_EvaluateDisplacement/4       52921      53034       9707
// BM_EvaluateDisplacement/8      101641     101800       5057
// BM_EvaluateDisplacement/15     186835     187269       2749
// BM_EvaluateDisplacement/16     199863     200781       2564
// BM_EvaluateDisplacement/17     211855     212203       2426
// BM_EvaluateDisplacement/18     225181     225494       2283
// BM_EvaluateDisplacement/19     236922     237237       2170

// .\Release\benchmarks.exe --benchmark_repetitions=10 --benchmark_min_time=2 --benchmark_filter=Newhall  // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 3310 MHz CPU
// 2015/05/24-13:16:32
// Benchmark                    Time(ns)    CPU(ns) Iterations
// -----------------------------------------------------------
// BM_NewhallApproximation/4         589        562    2000000
// BM_NewhallApproximation/8         657        624    2000000
// BM_NewhallApproximation/16        754        741    2000000

#include <random>
#include <vector>

#include "astronomy/frames.hpp"
#include "geometry/named_quantities.hpp"
#include "numerics/чебышёв_series.hpp"
#include "quantities/si.hpp"

// Must come last to avoid conflicts when defining the CHECK macros.
#include "benchmark/benchmark.h"

namespace principia {

using astronomy::ICRFJ2000Ecliptic;
using geometry::Displacement;
using quantities::si::Metre;
using quantities::si::Second;

namespace numerics {

namespace {
int const kEvaluationsPerIteration = 1000;
}  // namespace

void BM_EvaluateDouble(benchmark::State& state) {  // NOLINT(runtime/references)
  int const degree = state.range_x();
  std::mt19937_64 random(42);
  std::vector<double> coefficients;
  for (int i = 0; i <= degree; ++i) {
    coefficients.push_back(static_cast<double>(random()));
  }
  Instant const t0;
  Instant const t_min = t0 + static_cast<double>(random()) * Second;
  Instant const t_max = t_min + static_cast<double>(random()) * Second;
  ЧебышёвSeries<double> const series(coefficients, t_min, t_max);

  Instant t = t_min;
  Time const Δt = (t_max - t_min) * 1E-9;
  double result = 0.0;

  while (state.KeepRunning()) {
    for (int i = 0; i < kEvaluationsPerIteration; ++i) {
      result += series.Evaluate(t);
      t += Δt;
    }
  }

  // This weird call to |SetLabel| has no effect except that it uses |result|
  // and therefore prevents the loop from being optimized away.
  state.SetLabel(std::to_string(result).substr(0, 0));
}

void BM_EvaluateDisplacement(
    benchmark::State& state) {  // NOLINT(runtime/references)
  int const degree = state.range_x();
  std::mt19937_64 random(42);
  std::vector<Displacement<ICRFJ2000Ecliptic>> coefficients;
  for (int i = 0; i <= degree; ++i) {
    coefficients.push_back(
        Displacement<ICRFJ2000Ecliptic>(
            {static_cast<double>(random()) * Metre,
             static_cast<double>(random()) * Metre,
             static_cast<double>(random()) * Metre}));
  }
  Instant const t0;
  Instant const t_min = t0 + static_cast<double>(random()) * Second;
  Instant const t_max = t_min + static_cast<double>(random()) * Second;
  ЧебышёвSeries<Displacement<ICRFJ2000Ecliptic>> const series(
      coefficients, t_min, t_max);

  Instant t = t_min;
  Time const Δt = (t_max - t_min) * 1E-9;
  Displacement<ICRFJ2000Ecliptic> result{};

  while (state.KeepRunning()) {
    for (int i = 0; i < kEvaluationsPerIteration; ++i) {
      result += series.Evaluate(t);
      t += Δt;
    }
  }

  // This weird call to |SetLabel| has no effect except that it uses |result|
  // and therefore prevents the loop from being optimized away.
  std::stringstream ss;
  ss << result;
  state.SetLabel(ss.str().substr(0, 0));
}

void BM_NewhallApproximation(
    benchmark::State& state) {  // NOLINT(runtime/references)
  int const degree = state.range_x();
  std::mt19937_64 random(42);
  std::vector<double> p;
  std::vector<Variation<double>> v;
  Instant const t0;
  Instant const t_min = t0 + static_cast<double>(random()) * Second;
  Instant const t_max = t_min + static_cast<double>(random()) * Second;

  while (state.KeepRunning()) {
    state.PauseTiming();
    p.clear();
    v.clear();
    for (int i = 0; i <= 8; ++i) {
      p.push_back(static_cast<double>(static_cast<double>(random())));
      v.push_back(static_cast<double>(static_cast<double>(random())) / Second);
    }
    state.ResumeTiming();
    auto const series =
        ЧебышёвSeries<double>::NewhallApproximation(degree, p, v, t_min, t_max);
  }
}

BENCHMARK(BM_EvaluateDouble)->
    Arg(4)->Arg(8)->Arg(15)->Arg(16)->Arg(17)->Arg(18)->Arg(19);
BENCHMARK(BM_EvaluateDisplacement)->
    Arg(4)->Arg(8)->Arg(15)->Arg(16)->Arg(17)->Arg(18)->Arg(19);
BENCHMARK(BM_NewhallApproximation)->
    Arg(4)->Arg(8)->Arg(16);

}  // namespace benchmarks
}  // namespace principia
