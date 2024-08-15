// .\Release\x64\benchmarks.exe --benchmark_repetitions=3 --benchmark_filter=Planetarium  // NOLINT(whitespace/line_length)

#include <limits>
#include <memory>
#include <vector>

#include "astronomy/frames.hpp"
#include "astronomy/time_scales.hpp"
#include "base/not_null.hpp"
#include "base/status_utilities.hpp"  // 🧙 For CHECK_OK.
#include "benchmark/benchmark.h"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/instant.hpp"
#include "geometry/interval.hpp"
#include "geometry/orthogonal_map.hpp"
#include "geometry/perspective.hpp"
#include "geometry/rotation.hpp"
#include "geometry/sign.hpp"
#include "geometry/signature.hpp"
#include "geometry/space.hpp"
#include "geometry/space_transformations.hpp"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/planetarium.hpp"
#include "physics/body_centred_body_direction_reference_frame.hpp"
#include "physics/body_centred_non_rotating_reference_frame.hpp"
#include "physics/body_surface_reference_frame.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/massive_body.hpp"
#include "physics/massless_body.hpp"
#include "physics/reference_frame.hpp"
#include "physics/rotating_body.hpp"
#include "physics/rotating_pulsating_reference_frame.hpp"
#include "physics/solar_system.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/solar_system_factory.hpp"

namespace principia {
namespace geometry {

using namespace principia::astronomy::_frames;
using namespace principia::astronomy::_time_scales;
using namespace principia::base::_not_null;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_instant;
using namespace principia::geometry::_interval;
using namespace principia::geometry::_orthogonal_map;
using namespace principia::geometry::_perspective;
using namespace principia::geometry::_rotation;
using namespace principia::geometry::_sign;
using namespace principia::geometry::_signature;
using namespace principia::geometry::_space;
using namespace principia::geometry::_space_transformations;
using namespace principia::integrators::_methods;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::ksp_plugin::_frames;
using namespace principia::ksp_plugin::_planetarium;
using namespace principia::physics::_body_centred_body_direction_reference_frame;  // NOLINT
using namespace principia::physics::_body_centred_non_rotating_reference_frame;
using namespace principia::physics::_body_surface_reference_frame;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_kepler_orbit;
using namespace principia::physics::_massive_body;
using namespace principia::physics::_massless_body;
using namespace principia::physics::_reference_frame;
using namespace principia::physics::_rotating_body;
using namespace principia::physics::_rotating_pulsating_reference_frame;
using namespace principia::physics::_solar_system;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_solar_system_factory;

namespace {

constexpr Length near = 40'000 * Kilo(Metre);
constexpr Length far = 400'000 * Kilo(Metre);
constexpr Length focal = 1 * Metre;

Perspective<Navigation, Camera> PolarPerspective(
    Similarity<Navigation, GCRS> const& navigation_to_gcrs_at_epoch,
    Length const distance_from_earth) {
  using LeftGCRS = Frame<struct LeftGCRSTag, Arbitrary, Handedness::Left>;
  return {RigidTransformation<GCRS, Camera>(
              GCRS::origin + Displacement<GCRS>(
                                 {0 * Metre, 0 * Metre, distance_from_earth}),
              Camera::origin,
              Rotation<LeftGCRS, Camera>(Vector<double, LeftGCRS>({1, 0, 0}),
                                         Vector<double, LeftGCRS>({0, -1, 0}),
                                         Bivector<double, LeftGCRS>({0, 0, -1}))
                      .Forget<OrthogonalMap>() *
                  Signature<GCRS, LeftGCRS>(Sign::Positive(),
                                            Sign::Positive(),
                                            DeduceSignReversingOrientation{})
                      .Forget<OrthogonalMap>())
                  .Forget<Similarity>() *
              navigation_to_gcrs_at_epoch,
          focal};
}

Perspective<Navigation, Camera> EquatorialPerspective(
    Similarity<Navigation, GCRS> const& navigation_to_gcrs_at_epoch,
    Length const distance_from_earth) {
  using LeftGCRS = Frame<struct LeftGCRSTag, Arbitrary, Handedness::Left>;
  return {RigidTransformation<GCRS, Camera>(
              GCRS::origin + Displacement<GCRS>(
                                 {0 * Metre, distance_from_earth, 0 * Metre}),
              Camera::origin,
              Rotation<LeftGCRS, Camera>(Vector<double, LeftGCRS>({1, 0, 0}),
                                         Vector<double, LeftGCRS>({0, 0, 1}),
                                         Bivector<double, LeftGCRS>({0, -1, 0}))
                      .Forget<OrthogonalMap>() *
                  Signature<GCRS, LeftGCRS>(Sign::Positive(),
                                            Sign::Positive(),
                                            DeduceSignReversingOrientation{})
                      .Forget<OrthogonalMap>())
                  .Forget<Similarity>() *
              navigation_to_gcrs_at_epoch,
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
        ephemeris_(solar_system_->MakeEphemeris(
            /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                                     /*geopotential_tolerance=*/0x1p-24},
            EphemerisParameters())),
        sun_(solar_system_->massive_body(
            *ephemeris_,
            SolarSystemFactory::name(SolarSystemFactory::Sun))),
        mercury_(solar_system_->massive_body(
            *ephemeris_,
            SolarSystemFactory::name(SolarSystemFactory::Mercury))),
        venus_(solar_system_->massive_body(
            *ephemeris_,
            SolarSystemFactory::name(SolarSystemFactory::Venus))),
        earth_(solar_system_->massive_body(
            *ephemeris_,
            SolarSystemFactory::name(SolarSystemFactory::Earth))),
        moon_(solar_system_->massive_body(
            *ephemeris_,
            SolarSystemFactory::name(SolarSystemFactory::Moon))),
        gcrs_(
            make_not_null_unique<
                BodyCentredNonRotatingReferenceFrame<Barycentric, GCRS>>(
                ephemeris_.get(),
            earth_)),
        earth_centred_inertial_(
            make_not_null_unique<
                BodyCentredNonRotatingReferenceFrame<Barycentric, Navigation>>(
                ephemeris_.get(),
                earth_)),
        earth_centred_earth_fixed_(
            make_not_null_unique<
                BodySurfaceReferenceFrame<Barycentric, Navigation>>(
                ephemeris_.get(),
                earth_)),
        geocentric_solar_ecliptic_(
            make_not_null_unique<
                BodyCentredBodyDirectionReferenceFrame<Barycentric,
                                                       Navigation>>(
                ephemeris_.get(),
                earth_,
                sun_)),
        earth_moon_lagrange_(
            make_not_null_unique<
                RotatingPulsatingReferenceFrame<Barycentric, Navigation>>(
                ephemeris_.get(),
                std::vector<not_null<MassiveBody const*>>{earth_},
                std::vector<not_null<MassiveBody const*>>{moon_})),
        sun_earth_lagrange_(
            make_not_null_unique<
                RotatingPulsatingReferenceFrame<Barycentric, Navigation>>(
                ephemeris_.get(),
                std::vector<not_null<MassiveBody const*>>{sun_,
                                                          mercury_,
                                                          venus_},
                std::vector<not_null<MassiveBody const*>>{earth_, moon_})) {
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

    CHECK_OK(ephemeris_->Prolong(goes_8_epoch));
    KeplerOrbit<Barycentric> const goes_8_orbit(
        *earth_, MasslessBody{}, goes_8_elements, goes_8_epoch);
    CHECK_OK(goes_8_trajectory_.Append(
        goes_8_epoch,
        ephemeris_->trajectory(earth_)->EvaluateDegreesOfFreedom(goes_8_epoch) +
            goes_8_orbit.StateVectors(goes_8_epoch)));
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
      Perspective<Navigation, Camera> const& perspective,
      not_null<PlottingFrame const*> plotting_frame) const {
    // No dark area, human visual acuity, wide field of view.
    Planetarium::Parameters parameters(
        /*sphere_radius_multiplier=*/1,
        /*angular_resolution=*/0.4 * ArcMinute,
        /*field_of_view=*/90 * Degree);
    Instant const t = goes_8_trajectory().front().time;
    Similarity<Navigation, GCRS> const plotting_to_gcrs =
        gcrs().ToThisFrameAtTimeSimilarly(t).similarity() *
        plotting_frame->FromThisFrameAtTimeSimilarly(t).similarity();
    return Planetarium(
        parameters,
        perspective,
        ephemeris_.get(),
        plotting_frame,
        [plotting_to_gcrs](Position<Navigation> const& plotted_point) {
          constexpr auto inverse_scale_factor = 1 / (6000 * Metre);
          return ScaledSpacePoint::FromCoordinates(
              ((plotting_to_gcrs(plotted_point) - GCRS::origin) *
               inverse_scale_factor)
                  .coordinates());
        });
  }

  ReferenceFrame<Barycentric, GCRS> const& gcrs() const {
    return *gcrs_;
  }

  PlottingFrame const& earth_centred_inertial() const {
    return *earth_centred_inertial_;
  }

  PlottingFrame const& earth_centred_earth_fixed() const {
    return *earth_centred_earth_fixed_;
  }

  PlottingFrame const& geocentric_solar_ecliptic() const {
    return *geocentric_solar_ecliptic_;
  }

  PlottingFrame const& earth_moon_lagrange() const {
    return *earth_moon_lagrange_;
  }

  PlottingFrame const& sun_earth_lagrange() const {
    return *sun_earth_lagrange_;
  }

 private:
  Ephemeris<Barycentric>::FixedStepParameters EphemerisParameters() {
    return Ephemeris<Barycentric>::FixedStepParameters(
        SymmetricLinearMultistepIntegrator<
            QuinlanTremaine1990Order12,
            Ephemeris<Barycentric>::NewtonianMotionEquation>(),
        /*step=*/10 * Minute);
  }

  Ephemeris<Barycentric>::FixedStepParameters HistoryParameters() {
    return Ephemeris<Barycentric>::FixedStepParameters(
        SymmetricLinearMultistepIntegrator<
            Quinlan1999Order8A,
            Ephemeris<Barycentric>::NewtonianMotionEquation>(),
        /*step=*/10 * Second);
  }

  not_null<std::unique_ptr<SolarSystem<Barycentric>>> const solar_system_;
  not_null<std::unique_ptr<Ephemeris<Barycentric>>> const ephemeris_;
  not_null<MassiveBody const*> const sun_;
  not_null<MassiveBody const*> const mercury_;
  not_null<MassiveBody const*> const venus_;
  not_null<RotatingBody<Barycentric> const*> const earth_;
  not_null<MassiveBody const*> const moon_;
  not_null<std::unique_ptr<ReferenceFrame<Barycentric, GCRS>>> const gcrs_;
  not_null<std::unique_ptr<PlottingFrame>> const earth_centred_inertial_;
  not_null<std::unique_ptr<PlottingFrame>> const earth_centred_earth_fixed_;
  not_null<std::unique_ptr<PlottingFrame>> const geocentric_solar_ecliptic_;
  not_null<std::unique_ptr<PlottingFrame>> const earth_moon_lagrange_;
  not_null<std::unique_ptr<PlottingFrame>> const sun_earth_lagrange_;
  DiscreteTrajectory<Barycentric> goes_8_trajectory_;
};

}  // namespace

void BM_PlanetariumPlotMethod3(
    benchmark::State& state,
    Perspective<Navigation, Camera> (*const perspective)(
        Similarity<Navigation, GCRS> const& navigation_to_gcrs_at_epoch,
        Length const distance_from_earth),
    Length const distance_from_earth,
    PlottingFrame const& (Satellites::*const plotting_frame)() const) {
  static Satellites satellites;
  Instant const t = satellites.goes_8_trajectory().front().time;
  PlottingFrame const& plotting = (satellites.*plotting_frame)();
  Planetarium planetarium = satellites.MakePlanetarium(
      perspective((satellites.gcrs().ToThisFrameAtTimeSimilarly(t) *
                   plotting.FromThisFrameAtTimeSimilarly(t)).similarity(),
                  distance_from_earth),
      &plotting);
  std::vector<ScaledSpacePoint> line;
  int iterations = 0;
  // This is the time of a lunar eclipse in January 2000.
  constexpr Instant now = "2000-01-21T04:41:30,5"_TT;
  for (auto _ : state) {
    line.clear();
    planetarium.PlotMethod3(
        satellites.goes_8_trajectory(),
        satellites.goes_8_trajectory().begin(),
        satellites.goes_8_trajectory().end(),
        now,
        /*t_max=*/InfiniteFuture,
        /*reverse=*/false,
        /*add_point=*/
        [&line](ScaledSpacePoint const& point) { line.push_back(point); },
        /*max_points=*/std::numeric_limits<int>::max());
    ++iterations;
  }
  Interval<double> x;
  Interval<double> y;
  Interval<double> z;
  for (auto const& point : line) {
    x.Include(point.x);
    y.Include(point.y);
    z.Include(point.z);
  }
  state.SetLabel((std::stringstream() << line.size() << " points within " << x
                                      << " × " << y << " × " << z)
                     .str());
}

#define PRINCIPIA_BENCHMARK_PLANETARIUM_PLOT_METHODS_NEAR_AND_FAR( \
    name, perspective, plotting_frame)                             \
  BENCHMARK_CAPTURE(BM_PlanetariumPlotMethod3,                     \
                    Near##name,                                    \
                    (perspective),                                 \
                    near,                                          \
                    (plotting_frame))                              \
      ->Unit(benchmark::kMillisecond);                             \
  BENCHMARK_CAPTURE(BM_PlanetariumPlotMethod3,                     \
                    Far##name,                                     \
                    (perspective),                                 \
                    far,                                           \
                    (plotting_frame))                              \
      ->Unit(benchmark::kMillisecond)

#define PRINCIPIA_BENCHMARK_PLANETARIUM_PLOT_METHODS_POLAR_AND_EQUATORIAL( \
    name, plotting_frame)                                                  \
  PRINCIPIA_BENCHMARK_PLANETARIUM_PLOT_METHODS_NEAR_AND_FAR(               \
      PolarPerspective##name, &PolarPerspective, (plotting_frame));        \
  PRINCIPIA_BENCHMARK_PLANETARIUM_PLOT_METHODS_NEAR_AND_FAR(               \
      EquatorialPerspective##name, &EquatorialPerspective, (plotting_frame))

PRINCIPIA_BENCHMARK_PLANETARIUM_PLOT_METHODS_POLAR_AND_EQUATORIAL(
    ECI,
    &Satellites::earth_centred_inertial);
PRINCIPIA_BENCHMARK_PLANETARIUM_PLOT_METHODS_POLAR_AND_EQUATORIAL(
    ECEF,
    &Satellites::earth_centred_earth_fixed);
PRINCIPIA_BENCHMARK_PLANETARIUM_PLOT_METHODS_POLAR_AND_EQUATORIAL(
    GSE,
    &Satellites::geocentric_solar_ecliptic);
PRINCIPIA_BENCHMARK_PLANETARIUM_PLOT_METHODS_POLAR_AND_EQUATORIAL(
    EML,
    &Satellites::earth_moon_lagrange);
PRINCIPIA_BENCHMARK_PLANETARIUM_PLOT_METHODS_POLAR_AND_EQUATORIAL(
    SEL,
    &Satellites::sun_earth_lagrange);

}  // namespace geometry
}  // namespace principia
