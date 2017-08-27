
// .\Release\x64\benchmarks.exe --benchmark_repetitions=3 --benchmark_filter=Ephemeris                                                                     // NOLINT(whitespace/line_length)

#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "astronomy/frames.hpp"
#include "base/not_null.hpp"
#include "base/thread_pool.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/rotation.hpp"
#include "integrators/integrators.hpp"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
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
using astronomy::ICRFJ200EquatorialToEcliptic;
using base::not_null;
using base::ThreadPool;
using geometry::Displacement;
using geometry::Instant;
using geometry::Position;
using geometry::Quaternion;
using geometry::Rotation;
using geometry::Velocity;
using integrators::Integrator;
using integrators::DormandElMikkawyPrince1986RKN434FM;
using integrators::McLachlanAtela1992Order5Optimal;
using integrators::Quinlan1999Order8A;
using integrators::QuinlanTremaine1990Order12;
using quantities::DebugString;
using quantities::Length;
using quantities::Speed;
using quantities::Sqrt;
using quantities::Time;
using quantities::astronomy::JulianYear;
using quantities::bipm::NauticalMile;
using quantities::si::AstronomicalUnit;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::SolarSystemFactory;

namespace physics {

namespace {

using Flow =
    void(not_null<DiscreteTrajectory<ICRFJ2000Equator>*> const trajectory,
         Instant const& t,
         Ephemeris<ICRFJ2000Equator>& ephemeris);

Length FittingTolerance(int const scale) {
  return 5 * std::pow(10.0, scale) * Metre;
}

Ephemeris<ICRFJ2000Equator>::FixedStepParameters EphemerisParameters() {
  return Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
      QuinlanTremaine1990Order12<Position<ICRFJ2000Equator>>(),
      /*step=*/10 * Minute);
}

void EphemerisSolarSystemBenchmark(SolarSystemFactory::Accuracy const accuracy,
                                   benchmark::State& state) {
  Length error;
  while (state.KeepRunning()) {
    state.PauseTiming();
    auto const at_спутник_1_launch =
        SolarSystemFactory::AtСпутник1Launch(accuracy);
    Instant const final_time = at_спутник_1_launch->epoch() + 100 * JulianYear;

    auto const ephemeris =
        at_спутник_1_launch->MakeEphemeris(FittingTolerance(state.range_x()),
                                           EphemerisParameters());

    state.ResumeTiming();
    ephemeris->Prolong(final_time);
    state.PauseTiming();
    error = (at_спутник_1_launch->trajectory(
                 *ephemeris,
                 SolarSystemFactory::name(SolarSystemFactory::Sun)).
                     EvaluatePosition(final_time) -
             at_спутник_1_launch->trajectory(
                 *ephemeris,
                 SolarSystemFactory::name(SolarSystemFactory::Earth)).
                     EvaluatePosition(final_time)).
                 Norm();
    state.ResumeTiming();
  }
  state.SetLabel(quantities::DebugString(error / AstronomicalUnit) + " ua");
}

template<Flow* flow>
void EphemerisL4ProbeBenchmark(SolarSystemFactory::Accuracy const accuracy,
                               Time const integration_duration,
                               benchmark::State& state) {
  Length sun_error;
  Length earth_error;
  int steps;

  auto const at_спутник_1_launch =
      SolarSystemFactory::AtСпутник1Launch(accuracy);
  Instant const final_time = at_спутник_1_launch->epoch() +
                             integration_duration;

  auto const ephemeris =
      at_спутник_1_launch->MakeEphemeris(FittingTolerance(state.range_x()),
                                         EphemerisParameters());

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
        at_спутник_1_launch->degrees_of_freedom(
            SolarSystemFactory::name(SolarSystemFactory::Sun));
    DegreesOfFreedom<ICRFJ2000Equator> const earth_degrees_of_freedom =
        at_спутник_1_launch->degrees_of_freedom(
            SolarSystemFactory::name(SolarSystemFactory::Earth));
    Displacement<ICRFJ2000Ecliptic> const sun_earth_displacement =
        ICRFJ200EquatorialToEcliptic(earth_degrees_of_freedom.position() -
                              sun_degrees_of_freedom.position());
    Rotation<ICRFJ2000Ecliptic, ICRFJ2000Ecliptic> const l4_rotation(
        Quaternion(cos(π / 6), {0, 0, sin(π / 6)}));
    Displacement<ICRFJ2000Ecliptic> const sun_l4_displacement =
        l4_rotation(sun_earth_displacement);
    Velocity<ICRFJ2000Ecliptic> const sun_earth_velocity =
        ICRFJ200EquatorialToEcliptic(earth_degrees_of_freedom.velocity() -
                              sun_degrees_of_freedom.velocity());
    Velocity<ICRFJ2000Ecliptic> const sun_l4_velocity =
        l4_rotation(sun_earth_velocity);
    trajectory.Append(at_спутник_1_launch->epoch(),
                      DegreesOfFreedom<ICRFJ2000Equator>(
                          sun_degrees_of_freedom.position() +
                              ICRFJ200EquatorialToEcliptic.Inverse()(
                                  sun_l4_displacement),
                          sun_degrees_of_freedom.velocity() +
                              ICRFJ200EquatorialToEcliptic.Inverse()(
                                  sun_l4_velocity)));

    state.ResumeTiming();
    flow(&trajectory, final_time, *ephemeris);
    state.PauseTiming();

    sun_error = (at_спутник_1_launch->trajectory(
                     *ephemeris,
                     SolarSystemFactory::name(SolarSystemFactory::Sun)).
                         EvaluatePosition(final_time) -
                 trajectory.last().degrees_of_freedom().position()).
                     Norm();
    earth_error = (at_спутник_1_launch->trajectory(
                       *ephemeris,
                       SolarSystemFactory::name(SolarSystemFactory::Earth)).
                           EvaluatePosition(final_time) -
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

template<Flow* flow>
void EphemerisLEOProbeBenchmark(SolarSystemFactory::Accuracy const accuracy,
                                benchmark::State& state) {
  Length sun_error;
  Length earth_error;
  int steps;

  auto const at_спутник_1_launch =
      SolarSystemFactory::AtСпутник1Launch(accuracy);
  Instant const final_time = at_спутник_1_launch->epoch() + 1 * JulianYear;

  auto const ephemeris =
      at_спутник_1_launch->MakeEphemeris(FittingTolerance(state.range_x()),
                                         EphemerisParameters());

  ephemeris->Prolong(final_time);

  while (state.KeepRunning()) {
    state.PauseTiming();
    // A probe in low earth orbit.
    MasslessBody probe;
    DiscreteTrajectory<ICRFJ2000Equator> trajectory;
    DegreesOfFreedom<ICRFJ2000Equator> const earth_degrees_of_freedom =
        at_спутник_1_launch->degrees_of_freedom(
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
    flow(&trajectory, final_time, *ephemeris);
    state.PauseTiming();

    sun_error = (at_спутник_1_launch->trajectory(
                     *ephemeris,
                     SolarSystemFactory::name(SolarSystemFactory::Sun)).
                         EvaluatePosition(final_time) -
                 trajectory.last().degrees_of_freedom().position()).
                     Norm();
    earth_error = (at_спутник_1_launch->trajectory(
                       *ephemeris,
                       SolarSystemFactory::name(SolarSystemFactory::Earth)).
                           EvaluatePosition(final_time) -
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

void BM_EphemerisMultithreadingBenchmark(benchmark::State& state) {
  auto const at_спутник_1_launch =
      SolarSystemFactory::AtСпутник1Launch(
          SolarSystemFactory::Accuracy::AllBodiesAndOblateness);
  Instant const epoch = at_спутник_1_launch->epoch();
  auto const ephemeris =
      at_спутник_1_launch->MakeEphemeris(1 * Milli(Metre),
                                         EphemerisParameters());
  std::string const& earth_name =
      SolarSystemFactory::name(SolarSystemFactory::Earth);
  auto const earth_massive_body =
      at_спутник_1_launch->massive_body(*ephemeris, earth_name);
  auto const earth_degrees_of_freedom =
      at_спутник_1_launch->degrees_of_freedom(earth_name);

  MasslessBody probe;
  std::list<DiscreteTrajectory<ICRFJ2000Equator>> trajectories;
  for (int i = 0; i < state.range_x(); ++i) {
    KeplerianElements<ICRFJ2000Equator> elements;
    elements.eccentricity = 0;
    elements.semimajor_axis = (i + 1) * 10'000 * Kilo(Metre);
    elements.inclination = 0 * Radian;
    elements.longitude_of_ascending_node = 0 * Radian;
    elements.argument_of_periapsis = 0 * Radian;
    elements.true_anomaly = 0 * Radian;
    KeplerOrbit<ICRFJ2000Equator> const orbit(
        *earth_massive_body, probe, elements, epoch);
    trajectories.emplace_back();
    auto& trajectory = trajectories.back();
    trajectory.Append(epoch,
                      earth_degrees_of_freedom + orbit.StateVectors(epoch));
  }

  ThreadPool<void> pool(/*pool_size=*/state.range_y());
  Instant final_time = epoch;
  while (state.KeepRunning()) {
    state.PauseTiming();
    std::vector<not_null<std::unique_ptr<Integrator<
        Ephemeris<ICRFJ2000Equator>::NewtonianMotionEquation>::Instance>>>
        instances;
    for (auto& trajectory : trajectories) {
      instances.push_back(ephemeris->NewInstance(
          {&trajectory},
          Ephemeris<ICRFJ2000Equator>::NoIntrinsicAccelerations,
          Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
              Quinlan1999Order8A<Position<ICRFJ2000Equator>>(),
              /*step=*/10 * Second)));
    }
    final_time += 100'000 * Second;
    state.ResumeTiming();

    std::vector<std::future<void>> futures;
    for (auto& instance : instances) {
      futures.push_back(pool.Add([&ephemeris, &instance, final_time]() {
        ephemeris->FlowWithFixedStep(final_time, *instance);
      }));
    }
    for (auto const& future : futures) {
      future.wait();
    }
  }

  std::stringstream ss;
  for (auto const& trajectory : trajectories) {
    Length const earth_distance =
        (at_спутник_1_launch->trajectory(*ephemeris, earth_name).
             EvaluatePosition(final_time) -
         trajectory.last().degrees_of_freedom().position()).Norm();
    ss << earth_distance << " ";
  }
  state.SetLabel(ss.str());
}

}  // namespace

void BM_EphemerisSolarSystemMajorBodiesOnly(benchmark::State& state) {
  EphemerisSolarSystemBenchmark(SolarSystemFactory::Accuracy::MajorBodiesOnly,
                                state);
}

void BM_EphemerisSolarSystemMinorAndMajorBodies(benchmark::State& state) {
  EphemerisSolarSystemBenchmark(
      SolarSystemFactory::Accuracy::MinorAndMajorBodies,
      state);
}

void BM_EphemerisSolarSystemAllBodiesAndOblateness(benchmark::State& state) {
  EphemerisSolarSystemBenchmark(
      SolarSystemFactory::Accuracy::AllBodiesAndOblateness,
      state);
}

template<Flow* flow>
void BM_EphemerisL4ProbeMajorBodiesOnly(benchmark::State& state) {
  EphemerisL4ProbeBenchmark<flow>(SolarSystemFactory::Accuracy::MajorBodiesOnly,
                                  /*integration_duration=*/100 * JulianYear,
                                  state);
}

template<Flow* flow>
void BM_EphemerisL4ProbeMinorAndMajorBodies(benchmark::State& state) {
  EphemerisL4ProbeBenchmark<flow>(
      SolarSystemFactory::Accuracy::MinorAndMajorBodies,
      /*integration_duration=*/100 * JulianYear,
      state);
}

template<Flow* flow>
void BM_EphemerisL4ProbeAllBodiesAndOblateness(benchmark::State& state) {
  EphemerisL4ProbeBenchmark<flow>(
      SolarSystemFactory::Accuracy::AllBodiesAndOblateness,
      /*integration_duration=*/100 * JulianYear,
      state);
}

template<Flow* flow>
void BM_EphemerisLEOProbeMajorBodiesOnly(benchmark::State& state) {
  EphemerisLEOProbeBenchmark<flow>(
      SolarSystemFactory::Accuracy::MajorBodiesOnly, state);
}

template<Flow* flow>
void BM_EphemerisLEOProbeMinorAndMajorBodies(benchmark::State& state) {
  EphemerisLEOProbeBenchmark<flow>(
      SolarSystemFactory::Accuracy::MinorAndMajorBodies, state);
}

template<Flow* flow>
void BM_EphemerisLEOProbeAllBodiesAndOblateness(benchmark::State& state) {
  EphemerisLEOProbeBenchmark<flow>(
      SolarSystemFactory::Accuracy::AllBodiesAndOblateness, state);
}

template<Flow* flow>
void BM_EphemerisFittingTolerance(benchmark::State& state) {
  EphemerisL4ProbeBenchmark<flow>(SolarSystemFactory::Accuracy::MajorBodiesOnly,
                                  /*integration_duration=*/100 * JulianYear,
                                  state);
}

template<Flow* flow>
void BM_EphemerisStartup(benchmark::State& state) {
  EphemerisL4ProbeBenchmark<flow>(SolarSystemFactory::Accuracy::MajorBodiesOnly,
                                  /*integration_duration=*/100 * Second,
                                  state);
}

void FlowEphemerisWithAdaptiveStep(
    not_null<DiscreteTrajectory<ICRFJ2000Equator>*> const trajectory,
    Instant const& t,
    Ephemeris<ICRFJ2000Equator>& ephemeris) {
  CHECK(ephemeris.FlowWithAdaptiveStep(
      trajectory,
      Ephemeris<ICRFJ2000Equator>::NoIntrinsicAcceleration,
      t,
      Ephemeris<ICRFJ2000Equator>::AdaptiveStepParameters(
          DormandElMikkawyPrince1986RKN434FM<Position<ICRFJ2000Equator>>(),
          /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
          /*length_integration_tolerance=*/1 * Metre,
          /*speed_integration_tolerance=*/1 * Metre / Second),
      Ephemeris<ICRFJ2000Equator>::unlimited_max_ephemeris_steps,
      /*last_point_only=*/false));
}

void FlowEphemerisWithFixedStepSLMS(
    not_null<DiscreteTrajectory<ICRFJ2000Equator>*> const trajectory,
    Instant const& t,
    Ephemeris<ICRFJ2000Equator>& ephemeris) {
  auto const instance = ephemeris.NewInstance(
      {trajectory},
      Ephemeris<ICRFJ2000Equator>::NoIntrinsicAccelerations,
      Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
          Quinlan1999Order8A<Position<ICRFJ2000Equator>>(),
          /*step=*/10 * Second));
  ephemeris.FlowWithFixedStep(t, *instance);
}

void FlowEphemerisWithFixedStepSRKN(
    not_null<DiscreteTrajectory<ICRFJ2000Equator>*> const trajectory,
    Instant const& t,
    Ephemeris<ICRFJ2000Equator>& ephemeris) {
  auto const instance = ephemeris.NewInstance(
      {trajectory},
      Ephemeris<ICRFJ2000Equator>::NoIntrinsicAccelerations,
      Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
          McLachlanAtela1992Order5Optimal<Position<ICRFJ2000Equator>>(),
          /*step=*/10 * Second));
  ephemeris.FlowWithFixedStep(t, *instance);
}

BENCHMARK(BM_EphemerisMultithreadingBenchmark)
    ->ArgPair(3, 1)
    ->ArgPair(3, 2)
    ->ArgPair(3, 3)
    ->ArgPair(3, 4)
    ->ArgPair(3, 5);
BENCHMARK(BM_EphemerisSolarSystemMajorBodiesOnly)->Arg(-3);
BENCHMARK(BM_EphemerisSolarSystemMinorAndMajorBodies)->Arg(-3);
BENCHMARK(BM_EphemerisSolarSystemAllBodiesAndOblateness)->Arg(-3);
BENCHMARK_TEMPLATE1(BM_EphemerisL4ProbeMajorBodiesOnly,
                    &FlowEphemerisWithAdaptiveStep)->Arg(-3);
BENCHMARK_TEMPLATE1(BM_EphemerisL4ProbeMinorAndMajorBodies,
                    &FlowEphemerisWithAdaptiveStep)->Arg(-3);
BENCHMARK_TEMPLATE1(BM_EphemerisL4ProbeAllBodiesAndOblateness,
                    &FlowEphemerisWithAdaptiveStep)->Arg(-3);
BENCHMARK_TEMPLATE1(BM_EphemerisLEOProbeMajorBodiesOnly,
                    &FlowEphemerisWithAdaptiveStep)->Arg(-3);
BENCHMARK_TEMPLATE1(BM_EphemerisLEOProbeMajorBodiesOnly,
                    &FlowEphemerisWithFixedStepSLMS)->Arg(-3);
BENCHMARK_TEMPLATE1(BM_EphemerisLEOProbeMajorBodiesOnly,
                    &FlowEphemerisWithFixedStepSRKN)->Arg(-3);
BENCHMARK_TEMPLATE1(BM_EphemerisLEOProbeMinorAndMajorBodies,
                    &FlowEphemerisWithAdaptiveStep)->Arg(-3);
BENCHMARK_TEMPLATE1(BM_EphemerisLEOProbeMinorAndMajorBodies,
                    &FlowEphemerisWithFixedStepSLMS)->Arg(-3);
BENCHMARK_TEMPLATE1(BM_EphemerisLEOProbeMinorAndMajorBodies,
                    &FlowEphemerisWithFixedStepSRKN)->Arg(-3);
BENCHMARK_TEMPLATE1(BM_EphemerisLEOProbeAllBodiesAndOblateness,
                    &FlowEphemerisWithAdaptiveStep)->Arg(-3);
BENCHMARK_TEMPLATE1(BM_EphemerisLEOProbeAllBodiesAndOblateness,
                    &FlowEphemerisWithFixedStepSLMS)->Arg(-3);
BENCHMARK_TEMPLATE1(BM_EphemerisLEOProbeAllBodiesAndOblateness,
                    &FlowEphemerisWithFixedStepSRKN)->Arg(-3);

BENCHMARK_TEMPLATE1(BM_EphemerisFittingTolerance,
                    &FlowEphemerisWithAdaptiveStep)->DenseRange(-4, 4);

BENCHMARK_TEMPLATE1(BM_EphemerisStartup,
                    &FlowEphemerisWithFixedStepSLMS)->Arg(3);
BENCHMARK_TEMPLATE1(BM_EphemerisStartup,
                    &FlowEphemerisWithFixedStepSRKN)->Arg(3);

}  // namespace physics
}  // namespace principia
