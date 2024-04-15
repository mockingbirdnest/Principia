// .\Release\x64\benchmarks.exe --benchmark_repetitions=10 --benchmark_filter=Newhall  // NOLINT(whitespace/line_length)

#include <memory>
#include <random>
#include <vector>

#include "astronomy/frames.hpp"
#include "base/not_null.hpp"
#include "benchmark/benchmark.h"
#include "geometry/instant.hpp"
#include "geometry/space.hpp"
#include "numerics/newhall.hpp"
#include "numerics/polynomial.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "numerics/polynomial_in_monomial_basis.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace numerics {

using namespace principia::astronomy::_frames;
using namespace principia::base::_not_null;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_space;
using namespace principia::numerics::_newhall;
using namespace principia::numerics::_polynomial;
using namespace principia::numerics::_polynomial_evaluators;
using namespace principia::numerics::_polynomial_in_monomial_basis;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_si;

void BM_NewhallApproximationDouble(benchmark::State& state) {
  int const degree = state.range(0);
  std::mt19937_64 random(42);
  std::vector<double> p;
  std::vector<Variation<double>> v;
  Instant const t0;
  Instant const t_min = t0 + static_cast<double>(random()) * Second;
  Instant const t_max = t_min + static_cast<double>(random()) * Second;

  double error_estimate;
  for (auto _ : state) {
    state.PauseTiming();
    p.clear();
    v.clear();
    for (int i = 0; i <= 8; ++i) {
      p.push_back(static_cast<double>(static_cast<double>(random())));
      v.push_back(static_cast<double>(static_cast<double>(random())) / Second);
    }
    state.ResumeTiming();
    auto const series = NewhallApproximationInMonomialBasis<double>(
        degree, p, v, t_min, t_max, Policy::AlwaysEstrin(), error_estimate);
  }
}

void BM_NewhallApproximationDisplacement(benchmark::State& state) {
  int const degree = state.range(0);
  std::mt19937_64 random(42);
  std::vector<Displacement<ICRS>> p;
  std::vector<Variation<Displacement<ICRS>>> v;
  Instant const t0;
  Instant const t_min = t0 + static_cast<double>(random()) * Second;
  Instant const t_max = t_min + static_cast<double>(random()) * Second;

  Displacement<ICRS> error_estimate;
  for (auto _ : state) {
    state.PauseTiming();
    p.clear();
    v.clear();
    for (int i = 0; i <= 8; ++i) {
      p.push_back(Displacement<ICRS>({static_cast<double>(random()) * Metre,
                                      static_cast<double>(random()) * Metre,
                                      static_cast<double>(random()) * Metre}));
      v.push_back(Variation<Displacement<ICRS>>(
          {static_cast<double>(random()) * Metre / Second,
           static_cast<double>(random()) * Metre / Second,
           static_cast<double>(random()) * Metre / Second}));
    }
    state.ResumeTiming();
    auto const series = NewhallApproximationInMonomialBasis<Displacement<ICRS>>(
        degree, p, v, t_min, t_max, Policy::AlwaysEstrin(), error_estimate);
  }
}

using ResultMonomialDouble =
    not_null<std::unique_ptr<Polynomial<double, Instant>>>;
using ResultMonomialDisplacement =
    not_null<std::unique_ptr<Polynomial<Displacement<ICRS>, Instant>>>;

BENCHMARK(BM_NewhallApproximationDouble)->Arg(4)->Arg(8)->Arg(16);
BENCHMARK(BM_NewhallApproximationDisplacement)->Arg(4)->Arg(8)->Arg(16);

}  // namespace numerics
}  // namespace principia
