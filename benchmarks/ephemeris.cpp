// .\Release\x64\benchmarks.exe --benchmark_repetitions=3 --benchmark_filter=Ephemeris                                                                     // NOLINT(whitespace/line_length)

#include <cmath>
#include <limits>
#include <list>
#include <memory>
#include <string>
#include <vector>

#include "astronomy/frames.hpp"
#include "astronomy/stabilize_ksp.hpp"
#include "base/not_null.hpp"
#include "base/status_utilities.hpp"
#include "base/thread_pool.hpp"
#include "benchmark/benchmark.h"
#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/quaternion.hpp"
#include "geometry/rotation.hpp"
#include "integrators/integrators.hpp"
#include "integrators/embedded_explicit_runge_kutta_nyström_integrator.hpp"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "ksp_plugin/frames.hpp"
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

namespace principia {
namespace physics {

using namespace principia::astronomy::_frames;
using namespace principia::base::_not_null;
using namespace principia::base::_thread_pool;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_identity;
using namespace principia::geometry::_named_quantities;
using namespace principia::geometry::_quaternion;
using namespace principia::geometry::_rotation;
using namespace principia::integrators::
    _embedded_explicit_runge_kutta_nyström_integrator;
using namespace principia::integrators::_integrators;
using namespace principia::integrators::_methods;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::integrators::
using namespace principia::ksp_plugin::_frames;
    _symplectic_runge_kutta_nyström_integrator;
using namespace principia::quantities::_astronomy;
using namespace principia::quantities::_bipm;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_solar_system_factory;

namespace {

using Flow = void(not_null<DiscreteTrajectory<Barycentric>*> const trajectory,
                  Instant const& t,
                  Ephemeris<Barycentric>& ephemeris);

Length FittingTolerance(int const scale) {
  return 5 * std::pow(10.0, scale) * Metre;
}

Ephemeris<Barycentric>::FixedStepParameters EphemerisParameters() {
  return Ephemeris<Barycentric>::FixedStepParameters(
      SymmetricLinearMultistepIntegrator<
          QuinlanTremaine1990Order12,
          Ephemeris<Barycentric>::NewtonianMotionEquation>(),
      /*step=*/10 * Minute);
}

not_null<std::unique_ptr<SolarSystem<Barycentric>>> SolarSystemAtСпутник1Launch(
    SolarSystemFactory::Accuracy const accuracy) {
  auto at_спутник_1_launch =
      make_not_null_unique<SolarSystem<Barycentric>>(
          SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
          SOLUTION_DIR / "astronomy" /
              "sol_initial_state_jd_2436145_604166667.proto.txt",
          /*ignore_frame=*/true);
  SolarSystemFactory::AdjustAccuracy(accuracy, *at_спутник_1_launch);
  return at_спутник_1_launch;
}

void BM_EphemerisKSPSystem(benchmark::State& state) {
  Length error;
  for (auto _ : state) {
    state.PauseTiming();

    auto at_origin = make_not_null_unique<SolarSystem<Barycentric>>(
        SOLUTION_DIR / "astronomy" / "kerbol_gravity_model.proto.txt",
        SOLUTION_DIR / "astronomy" / "kerbol_initial_state_0_0.proto.txt",
        /*ignore_frame=*/true);
    astronomy::StabilizeKSP(*at_origin);
    Instant const final_time = at_origin->epoch() + 100 * JulianYear;
    auto const ephemeris = at_origin->MakeEphemeris(
        /*accuracy_parameters=*/{FittingTolerance(state.range(0)),
                                 /*geopotential_tolerance=*/0x1p-24},
        Ephemeris<Barycentric>::FixedStepParameters(
            SymplecticRungeKuttaNyströmIntegrator<
                BlanesMoan2002SRKN14A,
                Ephemeris<Barycentric>::NewtonianMotionEquation>(),
            /*step=*/35 * Minute));

    state.ResumeTiming();
    CHECK_OK(ephemeris->Prolong(final_time));
    state.PauseTiming();
    error = (at_origin->trajectory(
                 *ephemeris,
                 SolarSystemFactory::name(SolarSystemFactory::Sun)).
                     EvaluatePosition(final_time) -
             at_origin->trajectory(
                 *ephemeris,
                 SolarSystemFactory::name(SolarSystemFactory::Earth)).
                     EvaluatePosition(final_time)).
                 Norm();
    state.ResumeTiming();
  }
  state.SetLabel(quantities::DebugString(error / AstronomicalUnit) + " ua");
}

template<SolarSystemFactory::Accuracy accuracy>
void BM_EphemerisSolarSystem(benchmark::State& state) {
  Length error;
  for (auto _ : state) {
    state.PauseTiming();

    auto const at_спутник_1_launch = SolarSystemAtСпутник1Launch(accuracy);
    Instant const final_time = at_спутник_1_launch->epoch() + 100 * JulianYear;
    auto const ephemeris =
        at_спутник_1_launch->MakeEphemeris(
            SolarSystemFactory::MakeAccuracyParameters<Barycentric>(
                FittingTolerance(state.range(0)),
                accuracy),
            EphemerisParameters());

    state.ResumeTiming();
    CHECK_OK(ephemeris->Prolong(final_time));
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

template<SolarSystemFactory::Accuracy accuracy, Flow* flow>
void BM_EphemerisLEOProbe(benchmark::State& state) {
  Length sun_error;
  Length earth_error;
  int steps;

  auto const at_спутник_1_launch = SolarSystemAtСпутник1Launch(accuracy);
  Instant const final_time = at_спутник_1_launch->epoch() + 1 * JulianYear;

  auto const ephemeris =
      at_спутник_1_launch->MakeEphemeris(
          SolarSystemFactory::MakeAccuracyParameters<Barycentric>(
              FittingTolerance(state.range(0)),
              accuracy),
          EphemerisParameters());

  CHECK_OK(ephemeris->Prolong(final_time));

  for (auto _ : state) {
    state.PauseTiming();
    // A probe in low earth orbit.
    MasslessBody probe;
    DiscreteTrajectory<Barycentric> trajectory;
    DegreesOfFreedom<Barycentric> const earth_degrees_of_freedom =
        at_спутник_1_launch->degrees_of_freedom(
            SolarSystemFactory::name(SolarSystemFactory::Earth));
    Displacement<Barycentric> const earth_probe_displacement(
        {6371 * Kilo(Metre) + 100 * NauticalMile, 0 * Metre, 0 * Metre});
    Speed const earth_probe_speed =
        Sqrt(at_спутник_1_launch->gravitational_parameter(
                 SolarSystemFactory::name(SolarSystemFactory::Earth)) /
                     earth_probe_displacement.Norm());
    Velocity<Barycentric> const earth_probe_velocity(
        {0 * Metre / Second, earth_probe_speed, 0 * Metre / Second});
    CHECK_OK(trajectory.Append(
        at_спутник_1_launch->epoch(),
        DegreesOfFreedom<Barycentric>(
            earth_degrees_of_freedom.position() + earth_probe_displacement,
            earth_degrees_of_freedom.velocity() + earth_probe_velocity)));

    state.ResumeTiming();
    flow(&trajectory, final_time, *ephemeris);
    state.PauseTiming();

    sun_error = (at_спутник_1_launch->trajectory(
                     *ephemeris,
                     SolarSystemFactory::name(SolarSystemFactory::Sun)).
                         EvaluatePosition(final_time) -
                 trajectory.back().degrees_of_freedom.position()).
                     Norm();
    earth_error = (at_спутник_1_launch->trajectory(
                       *ephemeris,
                       SolarSystemFactory::name(SolarSystemFactory::Earth)).
                           EvaluatePosition(final_time) -
                   trajectory.back().degrees_of_freedom.position()).
                       Norm();
    steps = trajectory.size();
    state.ResumeTiming();
  }
  std::stringstream ss;
  ss << steps;
  state.SetLabel(ss.str() + " steps, " +
                 quantities::DebugString(sun_error / AstronomicalUnit) +
                 " au, " +
                 quantities::DebugString((earth_error - 6371 * Kilo(Metre)) /
                                         NauticalMile) +
                 " nmi");
}

template<SolarSystemFactory::Accuracy accuracy, Flow* flow>
void BM_EphemerisTranslunarSpaceProbe(benchmark::State& state) {
  Length sun_error;
  Length earth_error;
  int steps;

  auto const at_спутник_1_launch = SolarSystemAtСпутник1Launch(accuracy);
  Instant const final_time = at_спутник_1_launch->epoch() + 1 * JulianYear;

  auto const ephemeris =
      at_спутник_1_launch->MakeEphemeris(
          SolarSystemFactory::MakeAccuracyParameters<Barycentric>(
              FittingTolerance(state.range(0)),
              accuracy),
          EphemerisParameters());

  CHECK_OK(ephemeris->Prolong(final_time));

  for (auto _ : state) {
    state.PauseTiming();
    // A probe orbiting the Earth beyond the orbit of the Moon.
    MasslessBody probe;
    DiscreteTrajectory<Barycentric> trajectory;
    DegreesOfFreedom<Barycentric> const earth_degrees_of_freedom =
        at_спутник_1_launch->degrees_of_freedom(
            SolarSystemFactory::name(SolarSystemFactory::Earth));
    Displacement<Barycentric> const earth_probe_displacement(
        {6371 * Kilo(Metre) + 500'000 * Kilo(Metre), 0 * Metre, 0 * Metre});
    Speed const earth_probe_speed =
        Sqrt(at_спутник_1_launch->gravitational_parameter(
                 SolarSystemFactory::name(SolarSystemFactory::Earth)) /
                     earth_probe_displacement.Norm());
    Velocity<Barycentric> const earth_probe_velocity(
        {0 * Metre / Second, earth_probe_speed, 0 * Metre / Second});
    CHECK_OK(trajectory.Append(
        at_спутник_1_launch->epoch(),
        DegreesOfFreedom<Barycentric>(
            earth_degrees_of_freedom.position() + earth_probe_displacement,
            earth_degrees_of_freedom.velocity() + earth_probe_velocity)));

    state.ResumeTiming();
    flow(&trajectory, final_time, *ephemeris);
    state.PauseTiming();

    sun_error = (at_спутник_1_launch->trajectory(
                     *ephemeris,
                     SolarSystemFactory::name(SolarSystemFactory::Sun)).
                         EvaluatePosition(final_time) -
                 trajectory.back().degrees_of_freedom.position()).
                     Norm();
    earth_error = (at_спутник_1_launch->trajectory(
                       *ephemeris,
                       SolarSystemFactory::name(SolarSystemFactory::Earth)).
                           EvaluatePosition(final_time) -
                   trajectory.back().degrees_of_freedom.position()).
                       Norm();
    steps = trajectory.size();
    state.ResumeTiming();
  }
  std::stringstream ss;
  ss << steps;
  state.SetLabel(ss.str() + " steps, " +
                 quantities::DebugString(sun_error / AstronomicalUnit) +
                 " au, " +
                 quantities::DebugString((earth_error - 6371 * Kilo(Metre)) /
                                         Kilo(Metre)) +
                 " km");
}

void BM_EphemerisMultithreadingBenchmark(benchmark::State& state) {
  auto const at_спутник_1_launch =
      SolarSystemAtСпутник1Launch(
          SolarSystemFactory::Accuracy::AllBodiesAndDampedOblateness);
  Instant const epoch = at_спутник_1_launch->epoch();
  auto const ephemeris =
      at_спутник_1_launch->MakeEphemeris(
          /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                                   /*geopotential_tolerance=*/0x1p-24},
          EphemerisParameters());
  std::string const& earth_name =
      SolarSystemFactory::name(SolarSystemFactory::Earth);
  auto const earth_massive_body =
      at_спутник_1_launch->massive_body(*ephemeris, earth_name);
  auto const earth_degrees_of_freedom =
      at_спутник_1_launch->degrees_of_freedom(earth_name);

  MasslessBody probe;
  std::list<DiscreteTrajectory<Barycentric>> trajectories;
  for (int i = 0; i < state.range(0); ++i) {
    KeplerianElements<Barycentric> elements;
    elements.eccentricity = 0;
    elements.semimajor_axis = (i + 1) * 10'000 * Kilo(Metre);
    elements.inclination = 0 * Radian;
    elements.longitude_of_ascending_node = 0 * Radian;
    elements.argument_of_periapsis = 0 * Radian;
    elements.true_anomaly = 0 * Radian;
    KeplerOrbit<Barycentric> const orbit(
        *earth_massive_body, probe, elements, epoch);
    trajectories.emplace_back();
    auto& trajectory = trajectories.back();
    CHECK_OK(trajectory.Append(
        epoch, earth_degrees_of_freedom + orbit.StateVectors(epoch)));
  }

  ThreadPool<void> pool(/*pool_size=*/state.range(1));
  static constexpr int warp_factor = 6E6;
  static constexpr Frequency refresh_frequency = 50 * Hertz;
  static constexpr Time step = warp_factor / refresh_frequency;
  Instant final_time = epoch;
  for (auto _ : state) {
    state.PauseTiming();
    std::vector<not_null<std::unique_ptr<Integrator<Ephemeris<
        Barycentric>::NewtonianMotionEquation>::Instance>>>
        instances;
    for (auto& trajectory : trajectories) {
      instances.push_back(ephemeris->NewInstance(
          {&trajectory},
          Ephemeris<Barycentric>::NoIntrinsicAccelerations,
          Ephemeris<Barycentric>::FixedStepParameters(
              SymmetricLinearMultistepIntegrator<
                  Quinlan1999Order8A,
                  Ephemeris<Barycentric>::NewtonianMotionEquation>(),
              /*step=*/10 * Second)));
    }
    final_time += step;
    state.ResumeTiming();

    std::vector<std::future<void>> futures;
    for (auto& instance : instances) {
      futures.push_back(pool.Add([&ephemeris, &instance, final_time]() {
        CHECK_OK(ephemeris->FlowWithFixedStep(final_time, *instance));
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
         trajectory.back().degrees_of_freedom.position()).Norm();
    ss << earth_distance << " ";
  }
  state.SetLabel(ss.str());
}

template<SolarSystemFactory::Accuracy accuracy, Flow* flow>
void EphemerisL4ProbeBenchmark(Time const integration_duration,
                               benchmark::State& state) {
  Length sun_error;
  Length earth_error;
  int steps;

  auto const at_спутник_1_launch = SolarSystemAtСпутник1Launch(accuracy);
  Instant const final_time = at_спутник_1_launch->epoch() +
                             integration_duration;

  auto const ephemeris =
      at_спутник_1_launch->MakeEphemeris(
          SolarSystemFactory::MakeAccuracyParameters<Barycentric>(
              FittingTolerance(state.range(0)),
              accuracy),
          EphemerisParameters());

  CHECK_OK(ephemeris->Prolong(final_time));

  auto make_l4_probe_trajectory = [&at_спутник_1_launch]() {
    // A probe near the L4 point of the Sun-Earth system.
    Identity<ICRS, Barycentric> to_barycentric;
    Identity<Barycentric, ICRS> from_barycentric;
    MasslessBody probe;
    auto trajectory = std::make_unique<DiscreteTrajectory<Barycentric>>();
    DegreesOfFreedom<Barycentric> const sun_degrees_of_freedom =
        at_спутник_1_launch->degrees_of_freedom(
            SolarSystemFactory::name(SolarSystemFactory::Sun));
    DegreesOfFreedom<Barycentric> const earth_degrees_of_freedom =
        at_спутник_1_launch->degrees_of_freedom(
            SolarSystemFactory::name(SolarSystemFactory::Earth));

    // The BCRS, but with non-ICRS axes; the x axis is the ICRS x axis, the xy
    // plane is the ecliptic realized by the obliquity given by the IAU in 1976
    // (16th general assembly, resolution 10, commission 4, recommendation 1),
    // identifying the ICRS xy plane with the mean equator of J2000.0.
    using Ecliptic = Frame<struct EclipticTag>;

    Rotation<ICRS, Ecliptic> const equatorial_to_ecliptic(
        23 * Degree + 26 * ArcMinute + 21.448 * ArcSecond,
        Bivector<double, ICRS>({1, 0, 0}),
        DefinesFrame<Ecliptic>{});
    auto const ecliptic_to_equatorial = equatorial_to_ecliptic.Inverse();

    Displacement<Ecliptic> const sun_earth_displacement =
        equatorial_to_ecliptic(
            from_barycentric(earth_degrees_of_freedom.position() -
                             sun_degrees_of_freedom.position()));
    Rotation<Ecliptic, Ecliptic> const l4_rotation(
        Quaternion(cos(π / 6), {0, 0, sin(π / 6)}));
    Displacement<Ecliptic> const sun_l4_displacement =
        l4_rotation(sun_earth_displacement);
    Velocity<Ecliptic> const sun_earth_velocity = equatorial_to_ecliptic(
        from_barycentric(earth_degrees_of_freedom.velocity() -
                         sun_degrees_of_freedom.velocity()));
    Velocity<Ecliptic> const sun_l4_velocity = l4_rotation(sun_earth_velocity);
    CHECK_OK(trajectory->Append(
        at_спутник_1_launch->epoch(),
        DegreesOfFreedom<Barycentric>(
            sun_degrees_of_freedom.position() +
                to_barycentric(ecliptic_to_equatorial(sun_l4_displacement)),
            sun_degrees_of_freedom.velocity() +
                to_barycentric(ecliptic_to_equatorial(sun_l4_velocity)))));
    return trajectory;
  };

  // Compute the total degree of the underlying polynomials.  Useful for
  // benchmarking the effect of the fitting tolerance.
  double total_degree = 0;
  for (auto const& body : ephemeris->bodies()) {
    total_degree += ephemeris->trajectory(body)->average_degree();
  }

  for (auto _ : state) {
    not_null<std::unique_ptr<DiscreteTrajectory<Barycentric>>> trajectory =
        make_l4_probe_trajectory();
    state.PauseTiming();

    state.ResumeTiming();
    flow(trajectory.get(), final_time, *ephemeris);
    state.PauseTiming();

    sun_error =
        (at_спутник_1_launch->trajectory(
             *ephemeris,
             SolarSystemFactory::name(
                 SolarSystemFactory::Sun)).EvaluatePosition(final_time) -
         trajectory->back().degrees_of_freedom.position()).Norm();
    earth_error =
        (at_спутник_1_launch->trajectory(
             *ephemeris,
             SolarSystemFactory::name(
                 SolarSystemFactory::Earth)).EvaluatePosition(final_time) -
         trajectory->back().degrees_of_freedom.position()).Norm();
    steps = trajectory->size();
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

}  // namespace

template<SolarSystemFactory::Accuracy accuracy, Flow* flow>
void BM_EphemerisL4Probe(benchmark::State& state) {
  EphemerisL4ProbeBenchmark<accuracy, flow>(
      /*integration_duration=*/100 * JulianYear, state);
}

template<SolarSystemFactory::Accuracy accuracy, Flow* flow>
void BM_EphemerisL4Probe1Year(benchmark::State& state) {
  EphemerisL4ProbeBenchmark<accuracy, flow>(
      /*integration_duration=*/1 * JulianYear, state);
}

template<Flow* flow>
void BM_EphemerisFittingTolerance(benchmark::State& state) {
  EphemerisL4ProbeBenchmark<SolarSystemFactory::Accuracy::MajorBodiesOnly,
                            flow>(/*integration_duration=*/100 * JulianYear,
                                  state);
}

template<Flow* flow>
void BM_EphemerisStartup(benchmark::State& state) {
  EphemerisL4ProbeBenchmark<SolarSystemFactory::Accuracy::MajorBodiesOnly,
                            flow>(/*integration_duration=*/100 * Second,
                                  state);
}

void FlowEphemerisWithAdaptiveStep(
    not_null<DiscreteTrajectory<Barycentric>*> const trajectory,
    Instant const& t,
    Ephemeris<Barycentric>& ephemeris) {
  CHECK_OK(ephemeris.FlowWithAdaptiveStep(
      trajectory,
      Ephemeris<Barycentric>::NoIntrinsicAcceleration,
      t,
      Ephemeris<Barycentric>::AdaptiveStepParameters(
          EmbeddedExplicitRungeKuttaNyströmIntegrator<
              DormandالمكاوىPrince1986RKN434FM,
              Ephemeris<Barycentric>::NewtonianMotionEquation>(),
          /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
          /*length_integration_tolerance=*/1 * Metre,
          /*speed_integration_tolerance=*/1 * Metre / Second),
      Ephemeris<Barycentric>::unlimited_max_ephemeris_steps));
}

void FlowEphemerisWithFixedStepSLMS(
    not_null<DiscreteTrajectory<Barycentric>*> const trajectory,
    Instant const& t,
    Ephemeris<Barycentric>& ephemeris) {
  auto const instance = ephemeris.NewInstance(
      {trajectory},
      Ephemeris<Barycentric>::NoIntrinsicAccelerations,
      Ephemeris<Barycentric>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<
              Quinlan1999Order8A,
              Ephemeris<Barycentric>::NewtonianMotionEquation>(),
          /*step=*/10 * Second));
  CHECK_OK(ephemeris.FlowWithFixedStep(t, *instance));
}

void FlowEphemerisWithFixedStepSRKN(
    not_null<DiscreteTrajectory<Barycentric>*> const trajectory,
    Instant const& t,
    Ephemeris<Barycentric>& ephemeris) {
  auto const instance = ephemeris.NewInstance(
      {trajectory},
      Ephemeris<Barycentric>::NoIntrinsicAccelerations,
      Ephemeris<Barycentric>::FixedStepParameters(
          SymplecticRungeKuttaNyströmIntegrator<
              McLachlanAtela1992Order5Optimal,
              Ephemeris<Barycentric>::NewtonianMotionEquation>(),
          /*step=*/10 * Second));
  CHECK_OK(ephemeris.FlowWithFixedStep(t, *instance));
}

BENCHMARK(BM_EphemerisMultithreadingBenchmark)
    ->ArgPair(3, 1)
    ->ArgPair(3, 2)
    ->ArgPair(3, 3)
    ->ArgPair(3, 4)
    ->ArgPair(3, 5)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_EphemerisKSPSystem)->Arg(-3)->Unit(benchmark::kSecond);
BENCHMARK_TEMPLATE(BM_EphemerisSolarSystem,
                   SolarSystemFactory::Accuracy::MajorBodiesOnly)
    ->Arg(-3)
    ->Unit(benchmark::kSecond);
BENCHMARK_TEMPLATE(BM_EphemerisSolarSystem,
                   SolarSystemFactory::Accuracy::MinorAndMajorBodies)
    ->Arg(-3)
    ->Unit(benchmark::kSecond);
BENCHMARK_TEMPLATE(BM_EphemerisSolarSystem,
                   SolarSystemFactory::Accuracy::AllBodiesAndDampedOblateness)
    ->Arg(-3)
    ->Unit(benchmark::kSecond);
BENCHMARK_TEMPLATE(BM_EphemerisL4Probe,
                   SolarSystemFactory::Accuracy::MajorBodiesOnly,
                   &FlowEphemerisWithAdaptiveStep)
    ->Arg(-3)
    ->Unit(benchmark::kSecond);
BENCHMARK_TEMPLATE(BM_EphemerisL4Probe,
                   SolarSystemFactory::Accuracy::MinorAndMajorBodies,
                   &FlowEphemerisWithAdaptiveStep)
    ->Arg(-3)
    ->Unit(benchmark::kSecond);
BENCHMARK_TEMPLATE(BM_EphemerisL4Probe,
                   SolarSystemFactory::Accuracy::AllBodiesAndDampedOblateness,
                   &FlowEphemerisWithAdaptiveStep)
    ->Arg(-3)
    ->Unit(benchmark::kSecond);
BENCHMARK_TEMPLATE(BM_EphemerisLEOProbe,
                   SolarSystemFactory::Accuracy::MajorBodiesOnly,
                   &FlowEphemerisWithAdaptiveStep)
    ->Arg(-3)
    ->Unit(benchmark::kSecond);
BENCHMARK_TEMPLATE(BM_EphemerisLEOProbe,
                   SolarSystemFactory::Accuracy::MajorBodiesOnly,
                   &FlowEphemerisWithFixedStepSLMS)
    ->Arg(-3)
    ->Unit(benchmark::kSecond);
BENCHMARK_TEMPLATE(BM_EphemerisLEOProbe,
                   SolarSystemFactory::Accuracy::MajorBodiesOnly,
                   &FlowEphemerisWithFixedStepSRKN)
    ->Arg(-3)
    ->Unit(benchmark::kSecond);
BENCHMARK_TEMPLATE(BM_EphemerisLEOProbe,
                   SolarSystemFactory::Accuracy::MinorAndMajorBodies,
                   &FlowEphemerisWithAdaptiveStep)
    ->Arg(-3)
    ->Unit(benchmark::kSecond);
BENCHMARK_TEMPLATE(BM_EphemerisLEOProbe,
                   SolarSystemFactory::Accuracy::MinorAndMajorBodies,
                   &FlowEphemerisWithFixedStepSLMS)
    ->Arg(-3)
    ->Unit(benchmark::kSecond);
BENCHMARK_TEMPLATE(BM_EphemerisLEOProbe,
                   SolarSystemFactory::Accuracy::MinorAndMajorBodies,
                   &FlowEphemerisWithFixedStepSRKN)
    ->Arg(-3)
    ->Unit(benchmark::kSecond);
BENCHMARK_TEMPLATE(BM_EphemerisLEOProbe,
                   SolarSystemFactory::Accuracy::AllBodiesAndDampedOblateness,
                   &FlowEphemerisWithAdaptiveStep)
    ->Arg(-3)
    ->Unit(benchmark::kSecond);
BENCHMARK_TEMPLATE(BM_EphemerisLEOProbe,
                   SolarSystemFactory::Accuracy::AllBodiesAndDampedOblateness,
                   &FlowEphemerisWithFixedStepSLMS)
    ->Arg(-3)
    ->Unit(benchmark::kSecond);
BENCHMARK_TEMPLATE(BM_EphemerisLEOProbe,
                   SolarSystemFactory::Accuracy::AllBodiesAndDampedOblateness,
                   &FlowEphemerisWithFixedStepSRKN)
    ->Arg(-3)
    ->Unit(benchmark::kSecond);
BENCHMARK_TEMPLATE(BM_EphemerisTranslunarSpaceProbe,
                   SolarSystemFactory::Accuracy::MajorBodiesOnly,
                   &FlowEphemerisWithFixedStepSLMS)
    ->Arg(-3)
    ->Unit(benchmark::kSecond);
BENCHMARK_TEMPLATE(BM_EphemerisTranslunarSpaceProbe,
                   SolarSystemFactory::Accuracy::MinorAndMajorBodies,
                   &FlowEphemerisWithFixedStepSLMS)
    ->Arg(-3)
    ->Unit(benchmark::kSecond);
BENCHMARK_TEMPLATE(BM_EphemerisTranslunarSpaceProbe,
                   SolarSystemFactory::Accuracy::AllBodiesAndDampedOblateness,
                   &FlowEphemerisWithFixedStepSLMS)
    ->Arg(-3)
    ->Unit(benchmark::kSecond);
BENCHMARK_TEMPLATE(BM_EphemerisL4Probe1Year,
                   SolarSystemFactory::Accuracy::MajorBodiesOnly,
                   &FlowEphemerisWithFixedStepSLMS)
    ->Arg(-3)
    ->Unit(benchmark::kSecond);
BENCHMARK_TEMPLATE(BM_EphemerisL4Probe1Year,
                   SolarSystemFactory::Accuracy::MinorAndMajorBodies,
                   &FlowEphemerisWithFixedStepSLMS)
    ->Arg(-3)
    ->Unit(benchmark::kSecond);
BENCHMARK_TEMPLATE(BM_EphemerisL4Probe1Year,
                   SolarSystemFactory::Accuracy::AllBodiesAndDampedOblateness,
                   &FlowEphemerisWithFixedStepSLMS)
    ->Arg(-3)
    ->Unit(benchmark::kSecond);

BENCHMARK_TEMPLATE(BM_EphemerisFittingTolerance, &FlowEphemerisWithAdaptiveStep)
    ->DenseRange(-4, 4)
    ->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(BM_EphemerisStartup, &FlowEphemerisWithFixedStepSLMS)
    ->Arg(3)
    ->Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_EphemerisStartup, &FlowEphemerisWithFixedStepSRKN)
    ->Arg(3)
    ->Unit(benchmark::kMicrosecond);

}  // namespace physics
}  // namespace principia
