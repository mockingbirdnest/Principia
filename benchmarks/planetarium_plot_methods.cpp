
// .\Release\x64\benchmarks.exe --benchmark_repetitions=3 --benchmark_filter=Planetarium  // NOLINT(whitespace/line_length)

#include "ksp_plugin/planetarium.hpp"

#include <algorithm>

#include "astronomy/time_scales.hpp"
#include "benchmark/benchmark.h"
#include "physics/body_centred_non_rotating_dynamic_frame.hpp"
#include "physics/solar_system.hpp"
#include "testing_utilities/solar_system_factory.hpp"

namespace principia {
namespace geometry {

using astronomy::operator""_UT1;
using astronomy::operator""_TT;
using base::make_not_null_unique;
using base::not_null;
using geometry::Bivector;
using geometry::Perspective;
using geometry::RigidTransformation;
using geometry::RP2Lines;
using geometry::Vector;
using geometry::Velocity;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::Quinlan1999Order8A;
using integrators::methods::QuinlanTremaine1990Order12;
using ksp_plugin::Camera;
using ksp_plugin::Barycentric;
using ksp_plugin::Navigation;
using ksp_plugin::NavigationFrame;
using ksp_plugin::Planetarium;
using physics::BodyCentredNonRotatingDynamicFrame;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using physics::Ephemeris;
using physics::KeplerianElements;
using physics::KeplerOrbit;
using physics::MassiveBody;
using physics::MasslessBody;
using physics::SolarSystem;
using quantities::Angle;
using quantities::Cos;
using quantities::Infinity;
using quantities::Length;
using quantities::Sin;
using quantities::Time;
using quantities::si::ArcMinute;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::SolarSystemFactory;

namespace {

constexpr Length near = 40'000 * Kilo(Metre);
constexpr Length far = 400'000 * Kilo(Metre);
constexpr Length focal = 1 * Metre;

Perspective<Navigation, Camera> PolarPerspective(
    Length const distance_from_earth) {
  return {
      RigidTransformation<Navigation, Camera>(
          Navigation::origin + Displacement<Navigation>(
                                   {0 * Metre, 0 * Metre, distance_from_earth}),
          Camera::origin,
          Rotation<Navigation, Camera>(Vector<double, Navigation>({1, 0, 0}),
                                       Vector<double, Navigation>({0, -1, 0}),
                                       Bivector<double, Navigation>({0, 0, -1}))
              .Forget()),
      focal};
}

Perspective<Navigation, Camera> EquatorialPerspective(
    Length const distance_from_earth) {
  return {
      RigidTransformation<Navigation, Camera>(
          Navigation::origin + Displacement<Navigation>(
                                   {0 * Metre, distance_from_earth, 0 * Metre}),
          Camera::origin,
          Rotation<Navigation, Camera>(Vector<double, Navigation>({1, 0, 0}),
                                       Vector<double, Navigation>({0, 0, 1}),
                                       Bivector<double, Navigation>({0, -1, 0}))
              .Forget()),
      focal};
}

class Satellites {
 public:
  Satellites()
      : solar_system_(make_not_null_unique<SolarSystem<Barycentric>>(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2451545_000000000.proto.txt",
            /*ignore_frame=*/true)),
        ephemeris_(
            solar_system_->MakeEphemeris(/*fitting_tolerance=*/1 * Milli(Metre),
                                         EphemerisParameters())),
        earth_(solar_system_->massive_body(
            *ephemeris_,
            SolarSystemFactory::name(SolarSystemFactory::Earth))),
        earth_centred_inertial_(
            make_not_null_unique<
                BodyCentredNonRotatingDynamicFrame<Barycentric, Navigation>>(
                ephemeris_.get(),
                earth_)) {
    // Two-line elements for GOES-8:
    // 1 23051U 94022A   00004.06628221 -.00000243  00000-0  00000-0 0  9630
    // 2 23051   0.4232  97.7420 0004776 192.8349 121.5613  1.00264613 28364
    constexpr Instant goes_8_epoch = "JD2451548.56628221"_UT1;
    KeplerianElements<Barycentric> goes_8_elements;
    goes_8_elements.inclination = 0.4232 * Degree;
    goes_8_elements.longitude_of_ascending_node = 97.7420 * Degree;
    goes_8_elements.eccentricity = 0.0004776;
    goes_8_elements.argument_of_periapsis = 192.8349 * Degree;
    goes_8_elements.mean_anomaly = 121.5613 * Degree;
    goes_8_elements.mean_motion = 1.00264613 * (2 * π * Radian / Day);

    ephemeris_->Prolong(goes_8_epoch);
    KeplerOrbit<Barycentric> const goes_8_orbit(
        *earth_, MasslessBody{}, goes_8_elements, goes_8_epoch);
    goes_8_trajectory_.Append(
        goes_8_epoch,
        ephemeris_->trajectory(earth_)->EvaluateDegreesOfFreedom(goes_8_epoch) +
            goes_8_orbit.StateVectors(goes_8_epoch));
    auto goes_8_instance = ephemeris_->NewInstance(
        {&goes_8_trajectory_},
        Ephemeris<Barycentric>::NoIntrinsicAccelerations,
        HistoryParameters());
    CHECK_OK(ephemeris_->FlowWithFixedStep(goes_8_epoch + 100 * Day,
                                           *goes_8_instance));
  }

  DiscreteTrajectory<Barycentric> const& goes_8_trajectory() const {
    return goes_8_trajectory_;
  }

  Planetarium MakePlanetarium(
      Perspective<Navigation, Camera> const& perspective) const {
    // No dark area, human visual acuity, wide field of view.
    Planetarium::Parameters parameters(
        /*sphere_radius_multiplier=*/1,
        /*angular_resolution=*/0.4 * ArcMinute,
        /*field_of_view=*/90 * Degree);
    return Planetarium(parameters,
                       perspective,
                       ephemeris_.get(),
                       earth_centred_inertial_.get());
  }

 private:
  Ephemeris<Barycentric>::FixedStepParameters EphemerisParameters() {
    return Ephemeris<Barycentric>::FixedStepParameters(
        SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                           Position<Barycentric>>(),
        /*step=*/10 * Minute);
  }

  Ephemeris<Barycentric>::FixedStepParameters HistoryParameters() {
    return Ephemeris<Barycentric>::FixedStepParameters(
        SymmetricLinearMultistepIntegrator<Quinlan1999Order8A,
                                           Position<Barycentric>>(),
        /*step=*/10 * Second);
  }

  not_null<std::unique_ptr<SolarSystem<Barycentric>>> const solar_system_;
  not_null<std::unique_ptr<Ephemeris<Barycentric>>> const ephemeris_;
  not_null<MassiveBody const*> const earth_;
  not_null<std::unique_ptr<NavigationFrame>> const earth_centred_inertial_;
  DiscreteTrajectory<Barycentric> goes_8_trajectory_;
};

}  // namespace

void RunBenchmark(benchmark::State& state,
                  Perspective<Navigation, Camera> const& perspective) {
  Satellites satellites;
  Planetarium planetarium = satellites.MakePlanetarium(perspective);
  RP2Lines<Length, Camera> lines;
  int total_lines = 0;
  int iterations = 0;
  // This is the time of a lunar eclipse in January 2000.
  constexpr Instant now = "2000-01-21T04:41:30,5"_TT;
  while (state.KeepRunning()) {
    lines = planetarium.PlotMethod2(satellites.goes_8_trajectory().Begin(),
                                    satellites.goes_8_trajectory().End(),
                                    now,
                                    /*reverse=*/false);
    total_lines += lines.size();
    ++iterations;
  }
  Length min_x = Infinity<Length>();
  Length min_y = Infinity<Length>();
  Length max_x = -Infinity<Length>();
  Length max_y = -Infinity<Length>();
  int points = 0;
  for (auto const& line : lines) {
    points += line.size();
    for (auto const& point : line) {
      min_x = std::min(min_x, point.x());
      min_y = std::min(min_y, point.y());
      max_x = std::max(max_x, point.x());
      max_y = std::max(max_y, point.y());
    }
  }
  state.SetLabel(std::to_string(points) + " points in " +
                 std::to_string(total_lines / iterations) + " lines within [" +
                 DebugString(min_x) + ", " + DebugString(max_x) + "] × [" +
                 DebugString(min_y) + ", " + DebugString(max_y) + "]");
}

void BM_PlanetariumPlotMethod2NearPolarPerspective(benchmark::State& state) {
  RunBenchmark(state, PolarPerspective(near));
}

void BM_PlanetariumPlotMethod2FarPolarPerspective(benchmark::State& state) {
  RunBenchmark(state, PolarPerspective(far));
}

void BM_PlanetariumPlotMethod2NearEquatorialPerspective(
    benchmark::State& state) {
  RunBenchmark(state, EquatorialPerspective(near));
}

void BM_PlanetariumPlotMethod2FarEquatorialPerspective(
    benchmark::State& state) {
  RunBenchmark(state, EquatorialPerspective(far));
}

BENCHMARK(BM_PlanetariumPlotMethod2NearPolarPerspective);
BENCHMARK(BM_PlanetariumPlotMethod2FarPolarPerspective);
BENCHMARK(BM_PlanetariumPlotMethod2NearEquatorialPerspective);
BENCHMARK(BM_PlanetariumPlotMethod2FarEquatorialPerspective);

}  // namespace geometry
}  // namespace principia
