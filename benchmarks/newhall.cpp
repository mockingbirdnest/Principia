
// .\Release\x64\benchmarks.exe --benchmark_repetitions=10 --benchmark_min_time=2 --benchmark_filter=Newhall  // NOLINT(whitespace/line_length)

#include <memory>
#include <random>
#include <vector>

#include "astronomy/frames.hpp"
#include "base/not_null.hpp"
#include "benchmark/benchmark.h"
#include "geometry/named_quantities.hpp"
#include "numerics/newhall.hpp"
#include "numerics/polynomial.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/si.hpp"

namespace principia {

using astronomy::ICRS;
using base::not_null;
using geometry::Displacement;
using geometry::Instant;
using quantities::Variation;
using quantities::si::Metre;
using quantities::si::Second;

namespace numerics {

template<typename Result,
         Result (*newhall)(int degree,
                           std::vector<double> const& q,
                           std::vector<Variation<double>> const& v,
                           Instant const& t_min,
                           Instant const& t_max,
                           double& error_estimate)>
void BM_NewhallApproximationDouble(benchmark::State& state) {
  int const degree = state.range_x();
  std::mt19937_64 random(42);
  std::vector<double> p;
  std::vector<Variation<double>> v;
  Instant const t0;
  Instant const t_min = t0 + static_cast<double>(random()) * Second;
  Instant const t_max = t_min + static_cast<double>(random()) * Second;

  double error_estimate;
  while (state.KeepRunning()) {
    state.PauseTiming();
    p.clear();
    v.clear();
    for (int i = 0; i <= 8; ++i) {
      p.push_back(static_cast<double>(static_cast<double>(random())));
      v.push_back(static_cast<double>(static_cast<double>(random())) / Second);
    }
    state.ResumeTiming();
    auto const series = newhall(degree, p, v, t_min, t_max, error_estimate);
  }
}

template<typename Result,
         Result (*newhall)(int degree,
                           std::vector<Displacement<ICRS>> const& q,
                           std::vector<Variation<Displacement<ICRS>>> const& v,
                           Instant const& t_min,
                           Instant const& t_max,
                           Displacement<ICRS>& error_estimate)>
void BM_NewhallApproximationDisplacement(benchmark::State& state) {
  int const degree = state.range_x();
  std::mt19937_64 random(42);
  std::vector<Displacement<ICRS>> p;
  std::vector<Variation<Displacement<ICRS>>> v;
  Instant const t0;
  Instant const t_min = t0 + static_cast<double>(random()) * Second;
  Instant const t_max = t_min + static_cast<double>(random()) * Second;

  Displacement<ICRS> error_estimate;
  while (state.KeepRunning()) {
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
    auto const series = newhall(degree, p, v, t_min, t_max, error_estimate);
  }
}

using ResultЧебышёвDouble = ЧебышёвSeries<double>;
using ResultЧебышёвDisplacement = ЧебышёвSeries<Displacement<ICRS>>;
using ResultMonomialDouble =
    not_null<std::unique_ptr<Polynomial<double, Instant>>>;
using ResultMonomialDisplacement =
    not_null<std::unique_ptr<Polynomial<Displacement<ICRS>, Instant>>>;

BENCHMARK_TEMPLATE(
    BM_NewhallApproximationDouble,
    ResultЧебышёвDouble,
    &NewhallApproximationInЧебышёвBasis<double>)
    ->Arg(4)->Arg(8)->Arg(16);
BENCHMARK_TEMPLATE(
    BM_NewhallApproximationDisplacement,
    ResultЧебышёвDisplacement,
    &NewhallApproximationInЧебышёвBasis<Displacement<ICRS>>)
    ->Arg(4)->Arg(8)->Arg(16);
BENCHMARK_TEMPLATE(
    BM_NewhallApproximationDouble,
    ResultMonomialDouble,
    (&NewhallApproximationInMonomialBasis<double, EstrinEvaluator>))
    ->Arg(4)->Arg(8)->Arg(16);
BENCHMARK_TEMPLATE(
    BM_NewhallApproximationDisplacement,
    ResultMonomialDisplacement,
    (&NewhallApproximationInMonomialBasis<Displacement<ICRS>,
                                          EstrinEvaluator>))
    ->Arg(4)->Arg(8)->Arg(16);

}  // namespace numerics
}  // namespace principia
