
// .\Release\benchmarks.exe  --benchmark_repetitions=3 --benchmark_filter=Ephemeris                                                                     // NOLINT(whitespace/line_length)
// Benchmarking on 1 X 3310 MHz CPU
// 2015/06/06-20:50:07
// Benchmark                                              Time(ns)    CPU(ns) Iterations                                                                // NOLINT(whitespace/line_length)
// -------------------------------------------------------------------------------------                                                                // NOLINT(whitespace/line_length)
// BM_EphemerisSolarSystemMajorBodiesOnly               32679763294 32198606400          1                                 +1.00027592626789170e+00 ua  // NOLINT(whitespace/line_length)
// BM_EphemerisSolarSystemMajorBodiesOnly               32417363858 32105005800          1                                 +1.00027592626789170e+00 ua  // NOLINT(whitespace/line_length)
// BM_EphemerisSolarSystemMajorBodiesOnly               32999366021 32775810100          1                                 +1.00027592626789170e+00 ua  // NOLINT(whitespace/line_length)
// BM_EphemerisSolarSystemMajorBodiesOnly_mean          32698831058 32359807433          1                                 +1.00027592626789170e+00 ua  // NOLINT(whitespace/line_length)
// BM_EphemerisSolarSystemMajorBodiesOnly_stddev          237983633   296629884          0                                 +1.00027592626789170e+00 ua  // NOLINT(whitespace/line_length)
// BM_EphemerisSolarSystemMinorAndMajorBodies           62917942977 62291199300          1                                 +1.00027592631378660e+00 ua  // NOLINT(whitespace/line_length)
// BM_EphemerisSolarSystemMinorAndMajorBodies           62459306455 61791996100          1                                 +1.00027592631378660e+00 ua  // NOLINT(whitespace/line_length)
// BM_EphemerisSolarSystemMinorAndMajorBodies           62429231382 61901196800          1                                 +1.00027592631378660e+00 ua  // NOLINT(whitespace/line_length)
// BM_EphemerisSolarSystemMinorAndMajorBodies_mean      62602160271 61994797400          1                                 +1.00027592631378660e+00 ua  // NOLINT(whitespace/line_length)
// BM_EphemerisSolarSystemMinorAndMajorBodies_stddev      223629404   214276710          0                                 +1.00027592631378660e+00 ua  // NOLINT(whitespace/line_length)
// BM_EphemerisSolarSystemAllBodiesAndOblateness        76531402564 75816486000          1                                 +1.00027592630007800e+00 ua  // NOLINT(whitespace/line_length)
// BM_EphemerisSolarSystemAllBodiesAndOblateness        77371666885 75551284300          1                                 +1.00027592630007800e+00 ua  // NOLINT(whitespace/line_length)
// BM_EphemerisSolarSystemAllBodiesAndOblateness        77830555654 76128488000          1                                 +1.00027592630007800e+00 ua  // NOLINT(whitespace/line_length)
// BM_EphemerisSolarSystemAllBodiesAndOblateness_mean   77244541701 75832086100          1                                 +1.00027592630007800e+00 ua  // NOLINT(whitespace/line_length)
// BM_EphemerisSolarSystemAllBodiesAndOblateness_stddev   537940701   235900474          0                                 +1.00027592630007800e+00 ua  // NOLINT(whitespace/line_length)
#include <memory>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/rotation.hpp"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "physics/ephemeris.hpp"
#include "physics/massless_body.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/solar_system.hpp"

// Must come last to avoid conflicts when defining the CHECK macros.
#include "benchmark/benchmark.h"

namespace principia {

using astronomy::JulianYear;
using base::not_null;
using geometry::Position;
using geometry::Quaternion;
using geometry::Rotation;
using integrators::DormandElMikkawyPrince1986RKN434FM;
using integrators::McLachlanAtela1992Order5Optimal;
using physics::Ephemeris;
using physics::MasslessBody;
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
  while (state->KeepRunning()) {
    state->PauseTiming();
    not_null<std::unique_ptr<SolarSystem>> const at_спутник_1_launch =
        SolarSystem::AtСпутник1Launch(accuracy);
    std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies =
        at_спутник_1_launch->massive_bodies();
    std::vector<DegreesOfFreedom<ICRFJ2000Ecliptic>> const initial_state =
        at_спутник_1_launch->initial_state();

    Instant const initial_time = at_спутник_1_launch->time();
    Instant const final_time = initial_time + 100 * JulianYear;

    Ephemeris<ICRFJ2000Ecliptic>
        ephemeris(
            std::move(bodies),
            initial_state,
            at_спутник_1_launch->time(),
            McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Ecliptic>>(),
            45 * Minute,
            0.1 * Milli(Metre),
            5 * Milli(Metre));

    state->ResumeTiming();
    ephemeris.Prolong(final_time);
    state->PauseTiming();
    error = (ephemeris.trajectory(ephemeris.bodies()[SolarSystem::kSun]).
                 EvaluatePosition(final_time, nullptr) -
             ephemeris.trajectory(ephemeris.bodies()[SolarSystem::kEarth]).
                 EvaluatePosition(final_time, nullptr)).
                 Norm();
    state->ResumeTiming();
  }
  state->SetLabel(DebugString(error / AstronomicalUnit) + " ua");
}

void EphemerisL4ProbeBenchmark(SolarSystem::Accuracy const accuracy,
                               not_null<benchmark::State*> const state) {
  Length sun_error;
  Length earth_error;
  int steps;

  not_null<std::unique_ptr<SolarSystem>> const at_спутник_1_launch =
      SolarSystem::AtСпутник1Launch(accuracy);
  std::vector<not_null<std::unique_ptr<MassiveBody const>>> bodies =
      at_спутник_1_launch->massive_bodies();
  std::vector<DegreesOfFreedom<ICRFJ2000Ecliptic>> const initial_state =
      at_спутник_1_launch->initial_state();

  Instant const initial_time = at_спутник_1_launch->time();
  Instant const final_time = initial_time + 100 * JulianYear;

  Ephemeris<ICRFJ2000Ecliptic>
      ephemeris(
          std::move(bodies),
          initial_state,
          at_спутник_1_launch->time(),
          McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Ecliptic>>(),
          45 * Minute,
          0.1 * Milli(Metre),
          5 * Milli(Metre));

  ephemeris.Prolong(final_time);

  while (state->KeepRunning()) {
    state->PauseTiming();
    // A probe near the L4 point of the Sun-Earth system.
    MasslessBody probe;
    Trajectory<ICRFJ2000Ecliptic> trajectory(&probe);
    DegreesOfFreedom<ICRFJ2000Ecliptic> const sun_degrees_of_freedom =
        initial_state[SolarSystem::kSun];
    DegreesOfFreedom<ICRFJ2000Ecliptic> const earth_degrees_of_freedom =
        initial_state[SolarSystem::kEarth];
    Displacement<ICRFJ2000Ecliptic> const sun_earth_displacement =
        earth_degrees_of_freedom.position() -
        sun_degrees_of_freedom.position();
    Rotation<ICRFJ2000Ecliptic, ICRFJ2000Ecliptic> const l4_rotation(
        Quaternion(cos(π / 6), {0, 0, sin(π / 6)}));
    Displacement<ICRFJ2000Ecliptic> const sun_l4_displacement =
        l4_rotation(sun_earth_displacement);
    Velocity<ICRFJ2000Ecliptic> const sun_earth_velocity =
        earth_degrees_of_freedom.velocity() -
        sun_degrees_of_freedom.velocity();
    Velocity<ICRFJ2000Ecliptic> const sun_l4_velocity =
        l4_rotation(sun_earth_velocity);
    trajectory.Append(initial_time,
                      DegreesOfFreedom<ICRFJ2000Ecliptic>(
                          sun_degrees_of_freedom.position() +
                              sun_l4_displacement,
                          sun_degrees_of_freedom.velocity() + sun_l4_velocity));

    state->ResumeTiming();
    ephemeris.Flow(&trajectory,
                   1 * Metre,
                   1 * Metre / Second,
                   DormandElMikkawyPrince1986RKN434FM<
                       Position<ICRFJ2000Ecliptic>>(),
                   final_time);
    state->PauseTiming();

    sun_error = (ephemeris.trajectory(ephemeris.bodies()[SolarSystem::kSun]).
                     EvaluatePosition(final_time, nullptr) -
                 trajectory.last().degrees_of_freedom().position()).
                     Norm();
    earth_error = (ephemeris.trajectory(
                       ephemeris.bodies()[SolarSystem::kEarth]).
                           EvaluatePosition(final_time, nullptr) -
                   trajectory.last().degrees_of_freedom().position()).
                       Norm();
    steps = trajectory.Times().size();
    state->ResumeTiming();
  }
  std::stringstream ss;
  ss << steps;
  state->SetLabel(ss.str() + " steps, " +
                  DebugString(sun_error / AstronomicalUnit) + " ua, " + 
                  DebugString(earth_error / AstronomicalUnit) + " ua");
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

void BM_EphemerisL4ProbeMajorBodiesOnly(
    benchmark::State& state) {  // NOLINT(runtime/references)
  EphemerisL4ProbeBenchmark(SolarSystem::Accuracy::kMajorBodiesOnly,
                            &state);
}

void BM_EphemerisL4ProbeMinorAndMajorBodies(
    benchmark::State& state) {  // NOLINT(runtime/references)
  EphemerisL4ProbeBenchmark(SolarSystem::Accuracy::kMinorAndMajorBodies,
                            &state);
}

void BM_EphemerisL4ProbeAllBodiesAndOblateness(
    benchmark::State& state) {  // NOLINT(runtime/references)
  EphemerisL4ProbeBenchmark(SolarSystem::Accuracy::kAllBodiesAndOblateness,
                            &state);
}

BENCHMARK(BM_EphemerisSolarSystemMajorBodiesOnly);
BENCHMARK(BM_EphemerisSolarSystemMinorAndMajorBodies);
BENCHMARK(BM_EphemerisSolarSystemAllBodiesAndOblateness);
BENCHMARK(BM_EphemerisL4ProbeMajorBodiesOnly);
BENCHMARK(BM_EphemerisL4ProbeMinorAndMajorBodies);
BENCHMARK(BM_EphemerisL4ProbeAllBodiesAndOblateness);

}  // namespace benchmarks
}  // namespace principia
