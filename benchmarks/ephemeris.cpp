
// .\Release\benchmarks.exe  --benchmark_repetitions=3 --benchmark_filter=Solar

#include <memory>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "physics/ephemeris.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/solar_system.hpp"

// Must come last to avoid conflicts when defining the CHECK macros.
#include "benchmark/benchmark.h"

namespace principia {

using astronomy::JulianYear;
using base::not_null;
using geometry::Position;
using integrators::McLachlanAtela1992Order5Optimal;
using physics::Ephemeris;
using quantities::DebugString;
using si::AstronomicalUnit;
using si::Milli;
using si::Minute;
using testing_utilities::ICRFJ2000Ecliptic;
using testing_utilities::SolarSystem;

namespace benchmarks {

namespace {

void EphemerisSolarSystemBenchmark(SolarSystem::Accuracy const accuracy,
                                   not_null<benchmark::State*> const state) {
  Length error;
  state->PauseTiming();
  while (state->KeepRunning()) {
    not_null<std::unique_ptr<SolarSystem>> const at_спутник_1_launch =
        SolarSystem::AtСпутник1Launch(accuracy);
    not_null<std::unique_ptr<SolarSystem>> const at_спутник_2_launch =
        SolarSystem::AtСпутник2Launch(accuracy);
    std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies =
        at_спутник_1_launch->massive_bodies();
    std::vector<DegreesOfFreedom<ICRFJ2000Ecliptic>> const initial_state =
        at_спутник_1_launch->initial_state();

    std::vector<not_null<MassiveBody const*>> unowned_bodies;
    for (auto const& body : bodies) {
      unowned_bodies.push_back(body.get());
    }

    Ephemeris<ICRFJ2000Ecliptic>
        ephemeris(
            std::move(bodies),
            initial_state,
            at_спутник_1_launch->time(),
            McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Ecliptic>>(),
            45 * Minute,
            0.1 * Milli(Metre),
            5 * Milli(Metre));

    Instant final_time = at_спутник_1_launch->time() + 100 * JulianYear;
    state->ResumeTiming();
    ephemeris.Prolong(final_time);
    state->PauseTiming();
    error = (ephemeris.trajectory(unowned_bodies[SolarSystem::kSun]).
                 EvaluatePosition(final_time, nullptr) -
             ephemeris.trajectory(unowned_bodies[SolarSystem::kEarth]).
                 EvaluatePosition(final_time, nullptr)).
                 Norm();
  }
  state->SetLabel(DebugString(error / AstronomicalUnit) + " ua");
}

}  // namespace

void BM_EphemerisSolarSystemMajorBodiesOnly(
    benchmark::State& state) {  // NOLINT(runtime/references)
  EphemerisSolarSystemBenchmark(SolarSystem::Accuracy::kMajorBodiesOnly,
                                &state);
}

void BM_EphemerisSolarSystemMinorAndMajorBodies(
    benchmark::State& state) {  // NOLINT(runtime/references)
  EphemerisSolarSystemBenchmark(SolarSystem::Accuracy::kMinorAndMajorBodies,
                                &state);
}

void BM_EphemerisSolarSystemAllBodiesAndOblateness(
    benchmark::State& state) {  // NOLINT(runtime/references)
  EphemerisSolarSystemBenchmark(SolarSystem::Accuracy::kAllBodiesAndOblateness,
                                &state);
}

BENCHMARK(BM_EphemerisSolarSystemMajorBodiesOnly);
BENCHMARK(BM_EphemerisSolarSystemMinorAndMajorBodies);
BENCHMARK(BM_EphemerisSolarSystemAllBodiesAndOblateness);

}  // namespace benchmarks
}  // namespace principia
