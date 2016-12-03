
// .\Release\x64\benchmarks.exe --benchmark_repetitions=3 --benchmark_filter=Ephemeris                                                                     // NOLINT(whitespace/line_length)

#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "astronomy/frames.hpp"
#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/rotation.hpp"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/massless_body.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/bipm.hpp"
#include "quantities/elementary_functions.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/solar_system_factory.hpp"

// Must come last to avoid conflicts when defining the CHECK macros.
#include "benchmark/benchmark.h"

namespace principia {

using astronomy::ICRFJ2000Ecliptic;
using astronomy::ICRFJ2000Equator;
using astronomy::equatorial_to_ecliptic;
using base::not_null;
using geometry::Displacement;
using geometry::Position;
using geometry::Quaternion;
using geometry::Rotation;
using geometry::Velocity;
using integrators::DormandElMikkawyPrince1986RKN434FM;
using integrators::McLachlanAtela1992Order5Optimal;
using quantities::DebugString;
using quantities::Length;
using quantities::Speed;
using quantities::Sqrt;
using quantities::astronomy::JulianYear;
using quantities::bipm::NauticalMile;
using quantities::si::AstronomicalUnit;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Second;
using testing_utilities::SolarSystemFactory;

namespace physics {

namespace {

void EphemerisSolarSystemBenchmark(SolarSystemFactory::Accuracy const accuracy,
                                   benchmark::State& state) {
  Length const fitting_tolerance = 5 * std::pow(10.0, state.range_x()) * Metre;
  Length error;
  while (state.KeepRunning()) {
    state.PauseTiming();
    auto const at_спутник_1_launch =
        SolarSystemFactory::AtСпутник1Launch(accuracy);
    Instant const final_time = at_спутник_1_launch->epoch() + 100 * JulianYear;

    auto const ephemeris =
        at_спутник_1_launch->MakeEphemeris(
            fitting_tolerance,
            Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
                McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Equator>>(),
                /*step=*/45 * Minute));

    state.ResumeTiming();
    ephemeris->Prolong(final_time);
    state.PauseTiming();
    error = (at_спутник_1_launch->trajectory(
                 *ephemeris,
                 SolarSystemFactory::name(SolarSystemFactory::Sun)).
                     EvaluatePosition(final_time, nullptr) -
             at_спутник_1_launch->trajectory(
                 *ephemeris,
                 SolarSystemFactory::name(SolarSystemFactory::Earth)).
                     EvaluatePosition(final_time, nullptr)).
                 Norm();
    state.ResumeTiming();
  }
  state.SetLabel(quantities::DebugString(error / AstronomicalUnit) + " ua");
}

void EphemerisL4ProbeBenchmark(SolarSystemFactory::Accuracy const accuracy,
                               benchmark::State& state) {
  Length const fitting_tolerance = 5 * std::pow(10.0, state.range_x()) * Metre;
  Length sun_error;
  Length earth_error;
  int steps;

  auto const at_спутник_1_launch =
      SolarSystemFactory::AtСпутник1Launch(accuracy);
  Instant const final_time = at_спутник_1_launch->epoch() + 100 * JulianYear;

  auto const ephemeris =
      at_спутник_1_launch->MakeEphemeris(
          fitting_tolerance,
          Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
              McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Equator>>(),
              /*step=*/45 * Minute));

  ephemeris->Prolong(final_time);

  // Compute the total degree of the underlying polynomials.  Useful for
  // benchmarking the effect of the fitting tolerance.
  double total_degree = 0;
  for (auto const& body : ephemeris->bodies()) {
    total_degree += ephemeris->trajectory(body)->average_degree();
  }

  while (state.KeepRunning()) {
    state.PauseTiming();
    // A probe near the L4 point of the Sun-Earth system.
    MasslessBody probe;
    DiscreteTrajectory<ICRFJ2000Equator> trajectory;
    DegreesOfFreedom<ICRFJ2000Equator> const sun_degrees_of_freedom =
        at_спутник_1_launch->initial_state(
            SolarSystemFactory::name(SolarSystemFactory::Sun));
    DegreesOfFreedom<ICRFJ2000Equator> const earth_degrees_of_freedom =
        at_спутник_1_launch->initial_state(
            SolarSystemFactory::name(SolarSystemFactory::Earth));
    Displacement<ICRFJ2000Ecliptic> const sun_earth_displacement =
        equatorial_to_ecliptic(earth_degrees_of_freedom.position() -
                              sun_degrees_of_freedom.position());
    Rotation<ICRFJ2000Ecliptic, ICRFJ2000Ecliptic> const l4_rotation(
        Quaternion(cos(π / 6), {0, 0, sin(π / 6)}));
    Displacement<ICRFJ2000Ecliptic> const sun_l4_displacement =
        l4_rotation(sun_earth_displacement);
    Velocity<ICRFJ2000Ecliptic> const sun_earth_velocity =
        equatorial_to_ecliptic(earth_degrees_of_freedom.velocity() -
                              sun_degrees_of_freedom.velocity());
    Velocity<ICRFJ2000Ecliptic> const sun_l4_velocity =
        l4_rotation(sun_earth_velocity);
    trajectory.Append(at_спутник_1_launch->epoch(),
                      DegreesOfFreedom<ICRFJ2000Equator>(
                          sun_degrees_of_freedom.position() +
                              equatorial_to_ecliptic.Inverse()(
                                  sun_l4_displacement),
                          sun_degrees_of_freedom.velocity() +
                              equatorial_to_ecliptic.Inverse()(
                                  sun_l4_velocity)));

    state.ResumeTiming();
    ephemeris->FlowWithAdaptiveStep(
        &trajectory,
        Ephemeris<ICRFJ2000Equator>::NoIntrinsicAcceleration,
        final_time,
        Ephemeris<ICRFJ2000Equator>::AdaptiveStepParameters(
            DormandElMikkawyPrince1986RKN434FM<Position<ICRFJ2000Equator>>(),
            /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
            /*length_integration_tolerance=*/1 * Metre,
            /*speed_integration_tolerance=*/1 * Metre / Second),
        Ephemeris<ICRFJ2000Equator>::unlimited_max_ephemeris_steps);
    state.PauseTiming();

    sun_error = (at_спутник_1_launch->trajectory(
                     *ephemeris,
                     SolarSystemFactory::name(SolarSystemFactory::Sun)).
                         EvaluatePosition(final_time, nullptr) -
                 trajectory.last().degrees_of_freedom().position()).
                     Norm();
    earth_error = (at_спутник_1_launch->trajectory(
                       *ephemeris,
                       SolarSystemFactory::name(SolarSystemFactory::Earth)).
                           EvaluatePosition(final_time, nullptr) -
                   trajectory.last().degrees_of_freedom().position()).
                       Norm();
    steps = trajectory.Size();
    state.ResumeTiming();
  }
  std::stringstream ss;
  ss << steps;
  state.SetLabel(ss.str() + " steps, " +
                  quantities::DebugString(sun_error / AstronomicalUnit) +
                  " ua, " +
                  quantities::DebugString(earth_error / AstronomicalUnit) +
                  " ua, degree " +
                  std::to_string(total_degree));
}

void EphemerisLEOProbeBenchmark(SolarSystemFactory::Accuracy const accuracy,
                                benchmark::State& state) {
  Length const fitting_tolerance = 5 * std::pow(10.0, state.range_x()) * Metre;
  Length sun_error;
  Length earth_error;
  int steps;

  auto const at_спутник_1_launch =
      SolarSystemFactory::AtСпутник1Launch(accuracy);
  Instant const final_time = at_спутник_1_launch->epoch() + 1 * JulianYear;

  auto const ephemeris =
      at_спутник_1_launch->MakeEphemeris(
          fitting_tolerance,
          Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
              McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Equator>>(),
              /*step=*/45 * Minute));

  ephemeris->Prolong(final_time);

  while (state.KeepRunning()) {
    state.PauseTiming();
    // A probe in low earth orbit.
    MasslessBody probe;
    DiscreteTrajectory<ICRFJ2000Equator> trajectory;
    DegreesOfFreedom<ICRFJ2000Equator> const earth_degrees_of_freedom =
        at_спутник_1_launch->initial_state(
            SolarSystemFactory::name(SolarSystemFactory::Earth));
    Displacement<ICRFJ2000Equator> const earth_probe_displacement(
        {6371 * Kilo(Metre) + 100 * NauticalMile, 0 * Metre, 0 * Metre});
    Speed const earth_probe_speed =
        Sqrt(at_спутник_1_launch->gravitational_parameter(
                 SolarSystemFactory::name(SolarSystemFactory::Earth)) /
                     earth_probe_displacement.Norm());
    Velocity<ICRFJ2000Equator> const earth_probe_velocity(
        {0 * Metre / Second, earth_probe_speed, 0 * Metre / Second});
    trajectory.Append(at_спутник_1_launch->epoch(),
                      DegreesOfFreedom<ICRFJ2000Equator>(
                          earth_degrees_of_freedom.position() +
                              earth_probe_displacement,
                          earth_degrees_of_freedom.velocity() +
                              earth_probe_velocity));

    state.ResumeTiming();
    ephemeris->FlowWithAdaptiveStep(
        &trajectory,
        Ephemeris<ICRFJ2000Equator>::NoIntrinsicAcceleration,
        final_time,
        Ephemeris<ICRFJ2000Equator>::AdaptiveStepParameters(
            DormandElMikkawyPrince1986RKN434FM<Position<ICRFJ2000Equator>>(),
            /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
            /*length_integration_tolerance=*/1 * Metre,
            /*speed_integration_tolerance=*/1 * Metre / Second),
        Ephemeris<ICRFJ2000Equator>::unlimited_max_ephemeris_steps);
    state.PauseTiming();

    sun_error = (at_спутник_1_launch->trajectory(
                     *ephemeris,
                     SolarSystemFactory::name(SolarSystemFactory::Sun)).
                         EvaluatePosition(final_time, nullptr) -
                 trajectory.last().degrees_of_freedom().position()).
                     Norm();
    earth_error = (at_спутник_1_launch->trajectory(
                       *ephemeris,
                       SolarSystemFactory::name(SolarSystemFactory::Earth)).
                           EvaluatePosition(final_time, nullptr) -
                   trajectory.last().degrees_of_freedom().position()).
                       Norm();
    steps = trajectory.Size();
    state.ResumeTiming();
  }
  std::stringstream ss;
  ss << steps;
  state.SetLabel(ss.str() + " steps, " +
                 quantities::DebugString(sun_error / AstronomicalUnit) +
                 " ua, " +
                 quantities::DebugString((earth_error - 6371 * Kilo(Metre)) /
                                         NauticalMile) +
                 " nmi");
}

}  // namespace

void BM_EphemerisSolarSystemMajorBodiesOnly(
    benchmark::State& state) {  // NOLINT(runtime/references)
  EphemerisSolarSystemBenchmark(SolarSystemFactory::Accuracy::MajorBodiesOnly,
                                state);
}

void BM_EphemerisSolarSystemMinorAndMajorBodies(
    benchmark::State& state) {  // NOLINT(runtime/references)
  EphemerisSolarSystemBenchmark(
      SolarSystemFactory::Accuracy::MinorAndMajorBodies,
      state);
}

void BM_EphemerisSolarSystemAllBodiesAndOblateness(
    benchmark::State& state) {  // NOLINT(runtime/references)
  EphemerisSolarSystemBenchmark(
      SolarSystemFactory::Accuracy::AllBodiesAndOblateness,
      state);
}

void BM_EphemerisL4ProbeMajorBodiesOnly(
    benchmark::State& state) {  // NOLINT(runtime/references)
  EphemerisL4ProbeBenchmark(SolarSystemFactory::Accuracy::MajorBodiesOnly,
                            state);
}

void BM_EphemerisL4ProbeMinorAndMajorBodies(
    benchmark::State& state) {  // NOLINT(runtime/references)
  EphemerisL4ProbeBenchmark(SolarSystemFactory::Accuracy::MinorAndMajorBodies,
                            state);
}

void BM_EphemerisL4ProbeAllBodiesAndOblateness(
    benchmark::State& state) {  // NOLINT(runtime/references)
  EphemerisL4ProbeBenchmark(
      SolarSystemFactory::Accuracy::AllBodiesAndOblateness,
      state);
}

void BM_EphemerisLEOProbeMajorBodiesOnly(
    benchmark::State& state) {  // NOLINT(runtime/references)
  EphemerisLEOProbeBenchmark(SolarSystemFactory::Accuracy::MajorBodiesOnly,
                             state);
}

void BM_EphemerisLEOProbeMinorAndMajorBodies(
    benchmark::State& state) {  // NOLINT(runtime/references)
  EphemerisLEOProbeBenchmark(SolarSystemFactory::Accuracy::MinorAndMajorBodies,
                             state);
}

void BM_EphemerisLEOProbeAllBodiesAndOblateness(
    benchmark::State& state) {  // NOLINT(runtime/references)
  EphemerisLEOProbeBenchmark(
      SolarSystemFactory::Accuracy::AllBodiesAndOblateness,
      state);
}

void BM_EphemerisFittingTolerance(
    benchmark::State& state) {  // NOLINT(runtime/references)
  EphemerisL4ProbeBenchmark(SolarSystemFactory::Accuracy::MajorBodiesOnly,
                            state);
}

BENCHMARK(BM_EphemerisSolarSystemMajorBodiesOnly)->Arg(-3);
BENCHMARK(BM_EphemerisSolarSystemMinorAndMajorBodies)->Arg(-3);
BENCHMARK(BM_EphemerisSolarSystemAllBodiesAndOblateness)->Arg(-3);
BENCHMARK(BM_EphemerisL4ProbeMajorBodiesOnly)->Arg(-3);
BENCHMARK(BM_EphemerisL4ProbeMinorAndMajorBodies)->Arg(-3);
BENCHMARK(BM_EphemerisL4ProbeAllBodiesAndOblateness)->Arg(-3);
BENCHMARK(BM_EphemerisLEOProbeMajorBodiesOnly)->Arg(-3);
BENCHMARK(BM_EphemerisLEOProbeMinorAndMajorBodies)->Arg(-3);
BENCHMARK(BM_EphemerisLEOProbeAllBodiesAndOblateness)->Arg(-3);

BENCHMARK(BM_EphemerisFittingTolerance)->DenseRange(-4, 4);

}  // namespace physics
}  // namespace principia
