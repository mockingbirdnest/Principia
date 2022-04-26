
#include <algorithm>
#include <optional>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "astronomy/frames.hpp"
#include "base/not_null.hpp"
#include "geometry/named_quantities.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/integrators.hpp"
#include "integrators/methods.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "mathematica/mathematica.hpp"
#include "physics/continuous_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/rigid_motion.hpp"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/solar_system_factory.hpp"

namespace principia {

using base::dynamic_cast_not_null;
using geometry::AngleBetween;
using geometry::AngularVelocity;
using geometry::BarycentreCalculator;
using geometry::Bivector;
using geometry::Commutator;
using geometry::DefinesFrame;
using geometry::Displacement;
using geometry::Frame;
using geometry::Inertial;
using geometry::Instant;
using geometry::OrthogonalMap;
using geometry::Position;
using geometry::Rotation;
using geometry::Velocity;
using geometry::Wedge;
using integrators::FixedStepSizeIntegrator;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::SymplecticRungeKuttaNyströmIntegrator;
using integrators::methods::Quinlan1999Order8A;
using integrators::methods::QuinlanTremaine1990Order8;
using integrators::methods::QuinlanTremaine1990Order10;
using integrators::methods::QuinlanTremaine1990Order12;
using integrators::methods::BlanesMoan2002SRKN11B;
using integrators::methods::BlanesMoan2002SRKN14A;
using integrators::methods::McLachlanAtela1992Order5Optimal;
using physics::ContinuousTrajectory;
using physics::DegreesOfFreedom;
using physics::Ephemeris;
using physics::KeplerianElements;
using physics::KeplerOrbit;
using physics::RelativeDegreesOfFreedom;
using physics::RigidMotion;
using physics::RotatingBody;
using physics::SolarSystem;
using quantities::Angle;
using quantities::AngularMomentum;
using quantities::GravitationalParameter;
using quantities::Length;
using quantities::Time;
using quantities::astronomy::JulianYear;
using quantities::si::ArcSecond;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Hour;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AbsoluteError;
using testing_utilities::ApproximateQuantity;
using testing_utilities::IsNear;
using testing_utilities::RelativeError;
using testing_utilities::SolarSystemFactory;
using testing_utilities::operator""_;
using ::testing::Eq;
using ::testing::Lt;
using ::testing::Gt;

namespace astronomy {

class SolarSystemDynamicsTest : public ::testing::Test {
 protected:
  struct ActualOrbitError final {
    Angle separation_per_orbit;
    Angle inclination_drift_per_orbit;
    Angle longitude_of_ascending_node_drift_per_orbit;
    Angle argument_of_periapsis_drift_per_orbit;
  };

  struct ExpectedOrbitError final {
    ApproximateQuantity<Angle> separation_per_orbit;
    ApproximateQuantity<Angle> inclination_drift_per_orbit;
    ApproximateQuantity<Angle> longitude_of_ascending_node_drift_per_orbit;
    ApproximateQuantity<Angle> argument_of_periapsis_drift_per_orbit;
  };

  SolarSystemDynamicsTest() {
    google::LogToStderr();
    for (int primary = 0; primary <= SolarSystemFactory::LastBody; ++primary) {
      for (int i = SolarSystemFactory::Sun + 1;
           i <= SolarSystemFactory::LastBody;
           ++i) {
        if (SolarSystemFactory::parent(i) == primary) {
          bodies_orbiting_[primary].emplace_back(i);
        }
      }
    }
  }

  ActualOrbitError CompareOrbits(int const index,
                                 Ephemeris<ICRS> const& ephemeris,
                                 SolarSystem<ICRS> const& system,
                                 SolarSystem<ICRS> const& expected_system) {
    Instant const& epoch = expected_system.epoch();
    Time const duration = epoch - system.epoch();

    auto const name = SolarSystemFactory::name(index);
    auto const parent_name =
        SolarSystemFactory::name(SolarSystemFactory::parent(index));
    auto const body = system.massive_body(ephemeris, name);
    auto const parent = dynamic_cast_not_null<RotatingBody<ICRS> const*>(
        system.massive_body(ephemeris, parent_name));

    BarycentreCalculator<DegreesOfFreedom<ICRS>,
                         GravitationalParameter> actual_subsystem_barycentre;
    BarycentreCalculator<DegreesOfFreedom<ICRS>,
                         GravitationalParameter> expected_subsystem_barycentre;

    actual_subsystem_barycentre.Add(
        ephemeris.trajectory(body)->EvaluateDegreesOfFreedom(epoch),
        body->gravitational_parameter());
    expected_subsystem_barycentre.Add(expected_system.degrees_of_freedom(name),
                                      body->gravitational_parameter());
    for (int const moon_index : bodies_orbiting_[index]) {
      std::string  moon_name = SolarSystemFactory::name(moon_index);
      auto const moon = system.massive_body(ephemeris, moon_name);
      actual_subsystem_barycentre.Add(
          ephemeris.trajectory(moon)->EvaluateDegreesOfFreedom(epoch),
          moon->gravitational_parameter());
      expected_subsystem_barycentre.Add(
          expected_system.degrees_of_freedom(moon_name),
          moon->gravitational_parameter());
    }

    auto const actual_dof = actual_subsystem_barycentre.Get();
    auto const expected_dof = expected_subsystem_barycentre.Get();

    auto const actual_parent_dof =
        ephemeris.trajectory(parent)->EvaluateDegreesOfFreedom(epoch);
    auto const expected_parent_dof =
        expected_system.degrees_of_freedom(parent_name);

    // We transform to a frame in which |parent| has the z-axis as its rotation
    // axis by rotating around the normal to Earth's and |parent|'s rotation
    // axes.
    // If |parent| is the Sun, we use the normal to the invariable plane instead
    // of the Sun's axis.
    // TODO(egg): perhaps rotating bodies should export a rotation to their
    // celestial reference frame, we'll use that in the plugin too.
    using ParentEquator = Frame<enum class ParentEquatorTag, Inertial>;
    auto const z = Bivector<double, ICRS>({0, 0, 1});
    std::optional<Rotation<ICRS, ParentEquator>> rotation;

    if (SolarSystemFactory::parent(index) == SolarSystemFactory::Sun) {
      Bivector<AngularMomentum, ICRS> solar_system_angular_momentum;
      for (int i = SolarSystemFactory::Sun + 1;
           i <= SolarSystemFactory::LastBody;
           ++i) {
        auto const body_name = SolarSystemFactory::name(i);
        auto const body = system.massive_body(ephemeris, body_name);
        RelativeDegreesOfFreedom<ICRS> const from_solar_system_barycentre =
            expected_system.degrees_of_freedom(body_name) -
            DegreesOfFreedom<ICRS>(ICRS::origin, {});
        solar_system_angular_momentum +=
            Wedge(from_solar_system_barycentre.displacement(),
                  body->mass() * from_solar_system_barycentre.velocity()) *
            Radian;
      }
      Bivector<double, ICRS> const normal =
          Commutator(z, Normalize(solar_system_angular_momentum));
      // Check that we computed the invariable plane properly by computing its
      // angle to the Sun's equator.
      // The actual figure is "almost exactly 6 deg", see Bailey, Batygin, and
      // Brown, Solar Obliquity Induced by Planet Nine,
      // https://arxiv.org/pdf/1607.03963.pdf.
      EXPECT_THAT(AngleBetween(solar_system_angular_momentum,
                               parent->angular_velocity()),
                  IsNear(5.9_(1) * Degree));
      auto const declination_of_invariable_plane =
          OrientedAngleBetween(z, solar_system_angular_momentum, normal);
      EXPECT_THAT(declination_of_invariable_plane, Gt(0 * Radian));
      rotation = Rotation<ICRS, ParentEquator>(declination_of_invariable_plane,
                                               normal,
                                               DefinesFrame<ParentEquator>{});
    } else {
      auto const ω = parent->angular_velocity();
      Bivector<double, ICRS> const normal = Commutator(z, Normalize(ω));
      auto const parent_axis_declination = OrientedAngleBetween(z, ω, normal);
      EXPECT_THAT(parent_axis_declination, Gt(0 * Radian));
      rotation = Rotation<ICRS, ParentEquator>(
          parent_axis_declination, normal, DefinesFrame<ParentEquator>{});
    }
    RigidMotion<ICRS, ParentEquator> const to_parent_equator(
        {ICRS::origin,
         ParentEquator::origin,
         rotation->Forget<OrthogonalMap>()},
        /*angular_velocity_of_to_frame=*/ICRS::nonrotating,
        /*velocity_of_to_frame_origin=*/ICRS::unmoving);

    KeplerOrbit<ParentEquator> actual_osculating_orbit(
        /*primary=*/*parent,
        /*secondary=*/*body,
        /*state_vectors=*/to_parent_equator(actual_dof) -
            to_parent_equator(actual_parent_dof),
        epoch);
    KeplerOrbit<ParentEquator> expected_osculating_orbit(
        /*primary=*/*parent,
        /*secondary=*/*body,
        /*state_vectors=*/to_parent_equator(expected_dof) -
            to_parent_equator(expected_parent_dof),
        epoch);
    KeplerianElements<ParentEquator> const& actual_elements =
        actual_osculating_orbit.elements_at_epoch();
    KeplerianElements<ParentEquator> const& expected_elements =
        expected_osculating_orbit.elements_at_epoch();

    double const orbits = duration / *actual_elements.period;

    ActualOrbitError result;
    result.separation_per_orbit =
        AngleBetween(actual_dof.position() - actual_parent_dof.position(),
                     expected_dof.position() - expected_parent_dof.position()) /
        orbits;
    result.inclination_drift_per_orbit =
        AbsoluteError(expected_elements.inclination,
                      actual_elements.inclination) / orbits;
    result.longitude_of_ascending_node_drift_per_orbit =
        AbsoluteError(expected_elements.longitude_of_ascending_node,
                      actual_elements.longitude_of_ascending_node) / orbits;

    result.argument_of_periapsis_drift_per_orbit =
        AbsoluteError(*expected_elements.argument_of_periapsis,
                      *actual_elements.argument_of_periapsis) / orbits;
    return result;
  }

  std::map<int, std::vector<int>> bodies_orbiting_;
};

#if !_DEBUG
TEST_F(SolarSystemDynamicsTest, DISABLED_TenYearsFromJ2000) {
  SolarSystem<ICRS> solar_system_at_j2000(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2451545_000000000.proto.txt");

  SolarSystem<ICRS> ten_years_later(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2455200_500000000.proto.txt");

  // NOTE(phl): Keep these parameters aligned with
  // sol_numerics_blueprint.proto.txt.
  auto const ephemeris = solar_system_at_j2000.MakeEphemeris(
      /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<ICRS>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                             Position<ICRS>>(),
          /*step=*/10 * Minute));
  EXPECT_OK(ephemeris->Prolong(ten_years_later.epoch()));

  // NOTE(phl):
  // For Mercury and Venus the separation is about the order of magnitude we'd
  // expect from GR for either of those bodies; since it's a combination of
  // perihelion precession and change in anomaly, it's hard to get an exact
  // figure for that.
  // For Mercury the perihelion drift is what we expect from GR to within 1%.
  // For Pluto, WTF is wrong?
  std::map<int, ExpectedOrbitError> const expected_planet_orbit_errors{
      {SolarSystemFactory::Jupiter,
       {.separation_per_orbit = 0.033516_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.000013_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.0103_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 0.074_(1) * ArcSecond}},
      {SolarSystemFactory::Saturn,
       {.separation_per_orbit = 0.017580_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.000850_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.070520_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 0.168439_(1) * ArcSecond}},
      {SolarSystemFactory::Neptune,
       {.separation_per_orbit = 0.000356_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.000392_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.027661_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 0.465338_(1) * ArcSecond}},
      {SolarSystemFactory::Uranus,
       {.separation_per_orbit = 0.000137_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.000171_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.019118_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 0.007208_(1) * ArcSecond}},
      {SolarSystemFactory::Earth,
       {.separation_per_orbit = 0.085351_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.0000034_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.000101_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 0.036543_(1) * ArcSecond}},
      {SolarSystemFactory::Venus,
       {.separation_per_orbit = 0.104901_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.000001_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.000006_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 0.072214_(1) * ArcSecond}},
      {SolarSystemFactory::Mars,
       {.separation_per_orbit = 0.054834_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.0000027_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.000644_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 0.02552_(1) * ArcSecond}},
      {SolarSystemFactory::Mercury,
       {.separation_per_orbit = 0.189560_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.0000037_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.000020_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 0.102786_(1) * ArcSecond}},
      {SolarSystemFactory::Eris,
       {.separation_per_orbit = 0.000120_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.011912_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.005387_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 0.035167_(1) * ArcSecond}},
      {SolarSystemFactory::Pluto,
       {.separation_per_orbit = 26.521661_(1) * ArcSecond,
        .inclination_drift_per_orbit = 7.66074_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 42.3903_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 807.008_(1) * ArcSecond}},
      {SolarSystemFactory::Ceres,
       {.separation_per_orbit = 0.033067_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.000039_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.000956_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 0.001887_(1) * ArcSecond}},
      {SolarSystemFactory::Vesta,
       {.separation_per_orbit = 0.042217_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.000034_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.000103_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 0.021530_(1) * ArcSecond}},
  };

  for (int const planet_or_minor_planet :
       bodies_orbiting_[SolarSystemFactory::Sun]) {
    LOG(INFO) << "=== " << SolarSystemFactory::name(planet_or_minor_planet);
    auto const actual_orbit_error = CompareOrbits(planet_or_minor_planet,
                                                  *ephemeris,
                                                  solar_system_at_j2000,
                                                  ten_years_later);
    LOG(INFO) << "separation = " << std::fixed
              << actual_orbit_error.separation_per_orbit / ArcSecond
              << u8"″/orbit";
    LOG(INFO) << u8"Δi         = " << std::fixed
              << actual_orbit_error.inclination_drift_per_orbit / ArcSecond
              << u8"″/orbit";
    LOG(INFO)
        << u8"ΔΩ         = " << std::fixed
        << actual_orbit_error.longitude_of_ascending_node_drift_per_orbit /
               ArcSecond
        << u8"″/orbit";
    LOG(INFO) << u8"Δω         = " << std::fixed
              << actual_orbit_error.argument_of_periapsis_drift_per_orbit /
                     ArcSecond
              << u8"″/orbit";

    auto const& expected_orbit_error =
        expected_planet_orbit_errors.at(planet_or_minor_planet);
    EXPECT_THAT(actual_orbit_error.separation_per_orbit,
                IsNear(expected_orbit_error.separation_per_orbit));
    EXPECT_THAT(actual_orbit_error.inclination_drift_per_orbit,
                IsNear(expected_orbit_error.inclination_drift_per_orbit));
    EXPECT_THAT(
        actual_orbit_error.longitude_of_ascending_node_drift_per_orbit,
        IsNear(
            expected_orbit_error.longitude_of_ascending_node_drift_per_orbit));
    EXPECT_THAT(
        actual_orbit_error.argument_of_periapsis_drift_per_orbit,
        IsNear(expected_orbit_error.argument_of_periapsis_drift_per_orbit));
  }

  std::map<int, ExpectedOrbitError> const expected_moon_orbit_errors{
      {SolarSystemFactory::Ganymede,
       {.separation_per_orbit = 0.039263_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.001664_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.609354_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 0.63778_(1) * ArcSecond}},
      {SolarSystemFactory::Callisto,
       {.separation_per_orbit = 0.000271_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.000449_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.161792_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 0.16000_(1) * ArcSecond}},
      {SolarSystemFactory::Io,
       {.separation_per_orbit = 0.076865_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.000715_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 3.63026_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 3.6136_(1) * ArcSecond}},
      {SolarSystemFactory::Europa,
       {.separation_per_orbit = 0.065998_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.004845_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.100670_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 0.07751_(1) * ArcSecond}},

      {SolarSystemFactory::Titan,
       {.separation_per_orbit = 0.159615_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.000111_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.231945_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 0.123311_(1) * ArcSecond}},
      {SolarSystemFactory::Rhea,
       {.separation_per_orbit = 0.012813_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.004425_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.513725_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 0.496266_(1) * ArcSecond}},
      {SolarSystemFactory::Iapetus,
       {.separation_per_orbit = 0.039661_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.000035_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.001990_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 0.012217_(1) * ArcSecond}},
      {SolarSystemFactory::Dione,
       {.separation_per_orbit = 0.019140_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.000095_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.104560_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 0.117174_(1) * ArcSecond}},
      {SolarSystemFactory::Tethys,
       {.separation_per_orbit = 0.033548_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.000368_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.076115_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 0.045729_(1) * ArcSecond}},
      {SolarSystemFactory::Enceladus,
       {.separation_per_orbit = 0.029591_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.002010_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 6.84988_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 6.87334_(1) * ArcSecond}},
      {SolarSystemFactory::Mimas,
       {.separation_per_orbit = 0.021714_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.001474_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.113614_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 0.236115_(1) * ArcSecond}},

      {SolarSystemFactory::Triton,
       {.separation_per_orbit = 0.837046_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.038544_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.532907_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 0.52463_(1) * ArcSecond}},

      {SolarSystemFactory::Titania,
       {.separation_per_orbit = 0.070052_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.000030_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.716624_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 0.775136_(1) * ArcSecond}},
      {SolarSystemFactory::Oberon,
       {.separation_per_orbit = 0.008150_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.000407_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.130781_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 0.037069_(1) * ArcSecond}},
      {SolarSystemFactory::Ariel,
       {.separation_per_orbit = 0.004591_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.002529_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 20.1305_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 20.1339_(1) * ArcSecond}},
      {SolarSystemFactory::Umbriel,
       {.separation_per_orbit = 0.002270_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.000293_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 4.05936_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 4.06979_(1) * ArcSecond}},
      {SolarSystemFactory::Miranda,
       {.separation_per_orbit = 0.005774_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.001982_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.075551_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 0.078908_(1) * ArcSecond}},

      {SolarSystemFactory::Moon,
       {.separation_per_orbit = 0.167318_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.000776_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.000584_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 0.022414_(1) * ArcSecond}},

      {SolarSystemFactory::Phobos,
       {.separation_per_orbit = 1.023_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.0015_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.133_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 0.494_(1) * ArcSecond}},
      {SolarSystemFactory::Deimos,
       {.separation_per_orbit = 0.0068_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.002945_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.21581_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 0.213_(1) * ArcSecond}},

      {SolarSystemFactory::Charon,
       {.separation_per_orbit = 254.947_(1) * ArcSecond,
        .inclination_drift_per_orbit = 0.000085_(1) * ArcSecond,
        .longitude_of_ascending_node_drift_per_orbit = 0.1131_(1) * ArcSecond,
        .argument_of_periapsis_drift_per_orbit = 512.615_(1) * ArcSecond}},
  };

  for (int const planet_or_minor_planet :
       bodies_orbiting_[SolarSystemFactory::Sun]) {
    if (bodies_orbiting_[planet_or_minor_planet].empty()) {
      continue;
    }
    LOG(INFO) << ">>>>>> Moons of "
              << SolarSystemFactory::name(planet_or_minor_planet);
    for (int const moon : bodies_orbiting_[planet_or_minor_planet]) {
      LOG(INFO) << "=== " << SolarSystemFactory::name(moon);
      auto const actual_orbit_error = CompareOrbits(
          moon, *ephemeris, solar_system_at_j2000, ten_years_later);
      LOG(INFO) << "separation = " << std::fixed
                << actual_orbit_error.separation_per_orbit / ArcSecond
                << u8"″/orbit";
      LOG(INFO) << u8"Δi         = " << std::fixed
                << actual_orbit_error.inclination_drift_per_orbit / ArcSecond
                << u8"″/orbit";
      LOG(INFO)
          << u8"ΔΩ         = " << std::fixed
          << actual_orbit_error.longitude_of_ascending_node_drift_per_orbit /
                 ArcSecond
          << u8"″/orbit";
      LOG(INFO) << u8"Δω         = " << std::fixed
                << actual_orbit_error.argument_of_periapsis_drift_per_orbit /
                       ArcSecond
                << u8"″/orbit";

      auto const& expected_orbit_error = expected_moon_orbit_errors.at(moon);
      EXPECT_THAT(actual_orbit_error.separation_per_orbit,
                  IsNear(expected_orbit_error.separation_per_orbit));
      EXPECT_THAT(actual_orbit_error.inclination_drift_per_orbit,
                  IsNear(expected_orbit_error.inclination_drift_per_orbit));
      EXPECT_THAT(
          actual_orbit_error.longitude_of_ascending_node_drift_per_orbit,
          IsNear(expected_orbit_error
                     .longitude_of_ascending_node_drift_per_orbit));
      EXPECT_THAT(
          actual_orbit_error.argument_of_periapsis_drift_per_orbit,
          IsNear(expected_orbit_error.argument_of_periapsis_drift_per_orbit));
    }
  }
}

// This test produces the file phobos.generated.wl which is consumed by the
// notebook phobos.nb.
TEST(MarsTest, Phobos) {
  mathematica::Logger logger(TEMP_DIR / "phobos.generated.wl",
                             /*make_unique=*/false);

  SolarSystem<ICRS> solar_system_at_j2000(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2451545_000000000.proto.txt");
  // NOTE(phl): Keep these parameters aligned with
  // sol_numerics_blueprint.proto.txt.
  auto const ephemeris = solar_system_at_j2000.MakeEphemeris(
      /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<ICRS>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                             Position<ICRS>>(),
          /*step=*/10 * Minute));
  EXPECT_OK(ephemeris->Prolong(J2000 + 1 * JulianYear));

  ContinuousTrajectory<ICRS> const& mars_trajectory =
      solar_system_at_j2000.trajectory(*ephemeris, "Mars");

  std::vector<Position<ICRS>> mars_positions;
  std::vector<Velocity<ICRS>> mars_velocities;

  ContinuousTrajectory<ICRS> const& phobos_trajectory =
      solar_system_at_j2000.trajectory(*ephemeris, "Phobos");

  std::vector<Position<ICRS>> phobos_positions;
  std::vector<Velocity<ICRS>> phobos_velocities;

  for (Instant t = J2000; t < J2000 + 30 * Day; t += 5 * Minute) {
    mars_positions.push_back(mars_trajectory.EvaluatePosition(t));
    mars_velocities.push_back(mars_trajectory.EvaluateVelocity(t));
    phobos_positions.push_back(phobos_trajectory.EvaluatePosition(t));
    phobos_velocities.push_back(phobos_trajectory.EvaluateVelocity(t));
    logger.Append("ppaMarsPhobosDisplacements",
                  phobos_positions.back() - mars_positions.back());
    logger.Append("ppaMarsPhobosVelocities",
                  phobos_velocities.back() - mars_velocities.back());
  }
}

#endif

struct ConvergenceTestParameters {
  FixedStepSizeIntegrator<
      Ephemeris<ICRS>::NewtonianMotionEquation> const& integrator;
  int iterations;
  double first_step_in_seconds;
};

class SolarSystemDynamicsConvergenceTest
    : public ::testing::TestWithParam<ConvergenceTestParameters> {
 public:
  static void SetUpTestCase() {
    logger_ = new mathematica::Logger(
        SOLUTION_DIR / "mathematica" / "solar_system_convergence.generated.wl",
        /*make_unique=*/false);
  }

  static void TearDownTestCase() {
    delete logger_;
  }

 protected:
  FixedStepSizeIntegrator<Ephemeris<ICRS>::NewtonianMotionEquation> const&
  integrator() const {
    return GetParam().integrator;
  }

  int iterations() const {
    return GetParam().iterations;
  }

  int first_step_in_seconds() const {
    return GetParam().first_step_in_seconds;
  }

  static mathematica::Logger* logger_;
};

mathematica::Logger* SolarSystemDynamicsConvergenceTest::logger_ = nullptr;

// This takes 7-8 minutes to run.
TEST_P(SolarSystemDynamicsConvergenceTest, DISABLED_Convergence) {
  google::LogToStderr();
  Time const integration_duration = 1 * JulianYear;

  SolarSystem<ICRS> solar_system_at_j2000(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2451545_000000000.proto.txt");
  std::vector<std::string> const& body_names = solar_system_at_j2000.names();

  std::map<std::string, std::vector<DegreesOfFreedom<ICRS>>>
      name_to_degrees_of_freedom;
  std::vector<Time> steps;
  std::vector<std::chrono::duration<double>> durations;
  for (int i = 0; i < iterations(); ++i) {
    Time const step = first_step_in_seconds() * (1 << i) * Second;
    steps.push_back(step);
    LOG(INFO) << "Integrating with step " << step;

    auto const start = std::chrono::system_clock::now();
    auto const ephemeris = solar_system_at_j2000.MakeEphemeris(
        /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                                 /*geopotential_tolerance=*/0x1p-24},
        Ephemeris<ICRS>::FixedStepParameters(integrator(), step));
    EXPECT_OK(ephemeris->Prolong(solar_system_at_j2000.epoch() +
                                 integration_duration));
    auto const end = std::chrono::system_clock::now();
    durations.push_back(end - start);

    for (auto const& body_name : body_names) {
      auto const body =
          solar_system_at_j2000.massive_body(*ephemeris, body_name);
      name_to_degrees_of_freedom[body_name].emplace_back(
          ephemeris->trajectory(body)->EvaluateDegreesOfFreedom(
              solar_system_at_j2000.epoch() + integration_duration));
    }
  }

  std::map<std::string, std::vector<RelativeDegreesOfFreedom<ICRS>>>
      name_to_errors;
  for (auto const& [name, degrees_of_freedom] : name_to_degrees_of_freedom) {
    CHECK_EQ(degrees_of_freedom.size(), iterations());
    for (int i = 1; i < iterations(); ++i) {
      name_to_errors[name].emplace_back(degrees_of_freedom[i] -
                                        degrees_of_freedom[0]);
    }
  }

  std::vector<Length> position_errors(iterations() - 1);
  std::vector<std::string> worst_body(iterations() - 1);
  std::string const test_name(
      ::testing::UnitTest::GetInstance()->current_test_info()->name());
  for (int i = 0; i < iterations() - 1; ++i) {
    for (auto const& [name, errors] : name_to_errors) {
      if (position_errors[i] < errors[i].displacement().Norm()) {
        worst_body[i] = name;
      }
      position_errors[i] = std::max(position_errors[i],
                                    errors[i].displacement().Norm());
    }
    logger_->Append(std::string("ppaSolarSystemConvergence") +
                        test_name[test_name.size() - 1],
                    std::tuple(steps[i + 1],
                               position_errors[i],
                               worst_body[i],
                               durations[i + 1].count() * Second));
  }
}

INSTANTIATE_TEST_SUITE_P(
    AllSolarSystemDynamicsConvergenceTests,
    SolarSystemDynamicsConvergenceTest,
    ::testing::Values(
        ConvergenceTestParameters{
            .integrator =
                SymplecticRungeKuttaNyströmIntegrator<BlanesMoan2002SRKN11B,
                                                      Position<ICRS>>(),
            .iterations = 8,
            .first_step_in_seconds = 64},
        ConvergenceTestParameters{
            .integrator =
                SymplecticRungeKuttaNyströmIntegrator<
                    McLachlanAtela1992Order5Optimal,
                    Position<ICRS>>(),
             .iterations = 8,
             .first_step_in_seconds = 32},
        ConvergenceTestParameters{
            .integrator =
                SymmetricLinearMultistepIntegrator<Quinlan1999Order8A,
                                                   Position<ICRS>>(),
            .iterations = 6,
            .first_step_in_seconds = 64},
        ConvergenceTestParameters{
            .integrator =
                SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order8,
                                                   Position<ICRS>>(),
            .iterations = 6,
            .first_step_in_seconds = 64},
        ConvergenceTestParameters{
            .integrator =
                SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order10,
                                                   Position<ICRS>>(),
            .iterations = 6,
            .first_step_in_seconds = 64},

        // This is our favorite integrator.  For a step of 600 s it gives a
        // position error of about 2.3 m on Miranda and takes about 2.0 s of
        // elapsed time.  For steps larger than about 680 s, the errors explode.
        ConvergenceTestParameters{
            .integrator =
                SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                                   Position<ICRS>>(),
            .iterations = 5,
            .first_step_in_seconds = 75}));

}  // namespace astronomy
}  // namespace principia
