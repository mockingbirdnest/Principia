
#include <algorithm>
#include <experimental/filesystem>
#include <experimental/optional>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "astronomy/frames.hpp"
#include "base/file.hpp"
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
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/solar_system_factory.hpp"

namespace principia {

using base::dynamic_cast_not_null;
using base::OFStream;
using geometry::AngleBetween;
using geometry::BarycentreCalculator;
using geometry::Bivector;
using geometry::Commutator;
using geometry::DefinesFrame;
using geometry::Displacement;
using geometry::Frame;
using geometry::Instant;
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
using testing_utilities::IsNear;
using testing_utilities::SolarSystemFactory;
using ::testing::Lt;
using ::testing::Gt;

namespace astronomy {

class SolarSystemDynamicsTest : public testing::Test {
 protected:
  struct OrbitError final {
    Angle separation_per_orbit;
    Angle inclination_drift_per_orbit;
    std::experimental::optional<Angle>
        longitude_of_ascending_node_drift_per_orbit;
    std::experimental::optional<Angle>
        argument_of_periapsis_drift_per_orbit;
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

  OrbitError CompareOrbits(
      int const index,
      Ephemeris<ICRFJ2000Equator> const& ephemeris,
      SolarSystem<ICRFJ2000Equator> const& system,
      SolarSystem<ICRFJ2000Equator> const& expected_system) {
    Instant const& epoch = expected_system.epoch();
    Time const duration = epoch - system.epoch();

    auto const name = SolarSystemFactory::name(index);
    auto const parent_name =
        SolarSystemFactory::name(SolarSystemFactory::parent(index));
    auto const body = system.massive_body(ephemeris, name);
    auto const parent =
        dynamic_cast_not_null<RotatingBody<ICRFJ2000Equator> const*>(
            system.massive_body(ephemeris, parent_name));

    BarycentreCalculator<DegreesOfFreedom<ICRFJ2000Equator>,
                         GravitationalParameter> actual_subsystem_barycentre;
    BarycentreCalculator<DegreesOfFreedom<ICRFJ2000Equator>,
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
    enum LocalFrameTag { tag };
    using ParentEquator = Frame<LocalFrameTag, tag, /*frame_is_inertial=*/true>;
    auto const z = Bivector<double, ICRFJ2000Equator>({0, 0, 1});
    std::experimental::optional<Rotation<ICRFJ2000Equator, ParentEquator>>
        rotation;

    if (SolarSystemFactory::parent(index) == SolarSystemFactory::Sun) {
      Bivector<AngularMomentum, ICRFJ2000Equator>
          solar_system_angular_momentum;
      for (int i = SolarSystemFactory::Sun + 1;
           i <= SolarSystemFactory::LastBody;
           ++i) {
        auto const body_name = SolarSystemFactory::name(i);
        auto const body = system.massive_body(ephemeris, body_name);
        RelativeDegreesOfFreedom<ICRFJ2000Equator> const
            from_solar_system_barycentre =
                expected_system.degrees_of_freedom(body_name) -
                DegreesOfFreedom<ICRFJ2000Equator>(ICRFJ2000Equator::origin,
                                                   {});
        solar_system_angular_momentum +=
            Wedge(from_solar_system_barycentre.displacement(),
                  body->mass() * from_solar_system_barycentre.velocity()) *
            Radian;
      }
      Bivector<double, ICRFJ2000Equator> const normal =
          Commutator(z, Normalize(solar_system_angular_momentum));
      // Check that we computed the invariable plane properly by computing its
      // angle to the Sun's equator.
      // The actual figure is "almost exactly 6 deg", see Bailey, Batygin, and
      // Brown, Solar Obliquity Induced by Planet Nine,
      // https://arxiv.org/pdf/1607.03963.pdf.
      EXPECT_THAT(AngleBetween(solar_system_angular_momentum,
                               parent->angular_velocity()),
                  IsNear(5.9 * Degree));
      auto const declination_of_invariable_plane =
          OrientedAngleBetween(z, solar_system_angular_momentum, normal);
      EXPECT_THAT(declination_of_invariable_plane, Gt(0 * Radian));
      rotation =
          Rotation<ICRFJ2000Equator, ParentEquator>(
              declination_of_invariable_plane,
              normal,
              DefinesFrame<ParentEquator>{});
    } else {
      auto const ω = parent->angular_velocity();
      Bivector<double, ICRFJ2000Equator> const normal =
          Commutator(z, Normalize(ω));
      auto const parent_axis_declination = OrientedAngleBetween(z, ω, normal);
      EXPECT_THAT(parent_axis_declination, Gt(0 * Radian));
      rotation = Rotation<ICRFJ2000Equator, ParentEquator>(
          parent_axis_declination, normal, DefinesFrame<ParentEquator>{});
    }
    RigidMotion<ICRFJ2000Equator, ParentEquator> const
        to_parent_equator(
            {ICRFJ2000Equator::origin,
             ParentEquator::origin,
             rotation->Forget()},
            /*angular_velocity_of_to_frame=*/{},
            /*velocity_of_to_frame_origin=*/{});

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

    OrbitError result;
    result.separation_per_orbit =
        AngleBetween(actual_dof.position() - actual_parent_dof.position(),
                     expected_dof.position() - expected_parent_dof.position()) /
        orbits;
    result.inclination_drift_per_orbit =
        AbsoluteError(expected_elements.inclination,
                      actual_elements.inclination) / orbits;
    if (actual_elements.inclination > 0.1 * Degree) {
      result.longitude_of_ascending_node_drift_per_orbit =
          AbsoluteError(expected_elements.longitude_of_ascending_node,
                        actual_elements.longitude_of_ascending_node) / orbits;

      if (actual_elements.eccentricity > 0.1) {
        result.argument_of_periapsis_drift_per_orbit =
            AbsoluteError(*expected_elements.argument_of_periapsis,
                          *actual_elements.argument_of_periapsis) / orbits;
      }
    }
    return result;
  }

  std::map<int, std::vector<int>> bodies_orbiting_;
};

// This takes a minute to run.
TEST_F(SolarSystemDynamicsTest, DISABLED_TenYearsFromJ2000) {
  SolarSystem<ICRFJ2000Equator> solar_system_at_j2000(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2451545_000000000.proto.txt");

  SolarSystem<ICRFJ2000Equator> ten_years_later(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2455200_500000000.proto.txt");

  auto const ephemeris = solar_system_at_j2000.MakeEphemeris(
      /*fitting_tolerance=*/5 * Milli(Metre),
      Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
          SymplecticRungeKuttaNyströmIntegrator<BlanesMoan2002SRKN14A,
                                                Position<ICRFJ2000Equator>>(),
          /*step=*/45 * Minute));
  ephemeris->Prolong(ten_years_later.epoch());

  for (int const planet_or_minor_planet :
       bodies_orbiting_[SolarSystemFactory::Sun]) {
    LOG(INFO) << "=== " << SolarSystemFactory::name(planet_or_minor_planet);
    auto const error = CompareOrbits(planet_or_minor_planet,
                                     *ephemeris,
                                     solar_system_at_j2000,
                                     ten_years_later);
    LOG(INFO) << "separation = " << std::fixed
              << error.separation_per_orbit / ArcSecond << u8"″/orbit";
    LOG(INFO) << u8"Δi         = " << std::fixed
              << error.inclination_drift_per_orbit / ArcSecond << u8"″/orbit";
    if (error.longitude_of_ascending_node_drift_per_orbit) {
      LOG(INFO) << u8"ΔΩ         = " << std::fixed
                << *error.longitude_of_ascending_node_drift_per_orbit /
                       ArcSecond
                << u8"″/orbit";
    }
    if (error.argument_of_periapsis_drift_per_orbit) {
      LOG(INFO) << u8"Δω         = " << std::fixed
                << *error.argument_of_periapsis_drift_per_orbit /
                       ArcSecond
                << u8"″/orbit";
    }

    // This is about the order of magnitude we'd expect from GR for either of
    // those bodies; since it's a combination of perihelion precession and
    // change in anomaly, it's hard to get an exact figure for that.
    if (planet_or_minor_planet == SolarSystemFactory::Mercury ||
        planet_or_minor_planet == SolarSystemFactory::Venus) {
      EXPECT_THAT(error.separation_per_orbit,
                  IsNear(140 * Milli(ArcSecond), 2));
    } else if (planet_or_minor_planet != SolarSystemFactory::Pluto) {
      EXPECT_THAT(error.separation_per_orbit, Lt(100 * Milli(ArcSecond)));
    }

    if (error.argument_of_periapsis_drift_per_orbit) {
      switch (planet_or_minor_planet) {
        case SolarSystemFactory::Mercury:
          // This is what we expect from GR to the last sigfig.
          EXPECT_THAT(*error.argument_of_periapsis_drift_per_orbit,
                      IsNear(103 * Milli(ArcSecond), 1.01));
          break;
        case SolarSystemFactory::Eris:
          // I'm not sure what's going on with Eris; it's not clear what
          // ephemeris HORIZONS uses either.
          EXPECT_THAT(*error.argument_of_periapsis_drift_per_orbit,
                      Lt(50 * Milli(ArcSecond)));
          break;
        case SolarSystemFactory::Pluto:
          // WTF is wrong with Pluto?
          break;
        default:
          LOG(FATAL) << u8"Unexpected Δω for "
                     << SolarSystemFactory::name(planet_or_minor_planet);
      }
    }
    switch (planet_or_minor_planet) {
      case SolarSystemFactory::Pluto:
      case SolarSystemFactory::Eris:
        // Eris is likely from a non-integrated ephemeris; Pluto is mad.
        break;
      default:
        EXPECT_THAT(error.inclination_drift_per_orbit,
                    Lt(1 * Milli(ArcSecond)));
        break;
    }
  }

  // Moons.
  for (int const planet_or_minor_planet :
       bodies_orbiting_[SolarSystemFactory::Sun]) {
    if (bodies_orbiting_[planet_or_minor_planet].empty()) {
      continue;
    }
    LOG(INFO) << ">>>>>> Moons of "
              << SolarSystemFactory::name(planet_or_minor_planet);
    for (int const moon : bodies_orbiting_[planet_or_minor_planet]) {
      LOG(INFO) << "=== " << SolarSystemFactory::name(moon);
      auto const error = CompareOrbits(
          moon, *ephemeris, solar_system_at_j2000, ten_years_later);
      LOG(INFO) << "separation = " << std::fixed
                << error.separation_per_orbit / ArcSecond << u8"″/orbit";
      LOG(INFO) << u8"Δi         = " << std::fixed
                << error.inclination_drift_per_orbit / ArcSecond << u8"″/orbit";
      if (error.longitude_of_ascending_node_drift_per_orbit) {
        LOG(INFO) << u8"ΔΩ         = " << std::fixed
                  << *error.longitude_of_ascending_node_drift_per_orbit /
                         ArcSecond
                  << u8"″/orbit";
      }
      if (error.argument_of_periapsis_drift_per_orbit) {
        LOG(INFO) << u8"Δω         = " << std::fixed
                  << *error.argument_of_periapsis_drift_per_orbit / ArcSecond
                  << u8"″/orbit";
      }
    }
  }
}

#if !_DEBUG
// This test produces the file phobos.generated.wl which is consumed by the
// notebook phobos.nb.
TEST(MarsTest, Phobos) {
  SolarSystem<ICRFJ2000Equator> solar_system_at_j2000(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2451545_000000000.proto.txt");
  auto const ephemeris = solar_system_at_j2000.MakeEphemeris(
      /*fitting_tolerance=*/5 * Milli(Metre),
      Ephemeris<ICRFJ2000Equator>::FixedStepParameters(
          SymplecticRungeKuttaNyströmIntegrator<BlanesMoan2002SRKN14A,
                                                Position<ICRFJ2000Equator>>(),
          /*step=*/45 * Minute));
  ephemeris->Prolong(J2000 + 1 * JulianYear);

  ContinuousTrajectory<ICRFJ2000Equator> const& mars_trajectory =
      solar_system_at_j2000.trajectory(*ephemeris, "Mars");

  std::vector<Position<ICRFJ2000Equator>> mars_positions;
  std::vector<Velocity<ICRFJ2000Equator>> mars_velocities;

  ContinuousTrajectory<ICRFJ2000Equator> const& phobos_trajectory =
      solar_system_at_j2000.trajectory(*ephemeris, "Phobos");

  std::vector<Position<ICRFJ2000Equator>> phobos_positions;
  std::vector<Velocity<ICRFJ2000Equator>> phobos_velocities;

  std::vector<Displacement<ICRFJ2000Equator>> mars_phobos_displacements;
  std::vector<Velocity<ICRFJ2000Equator>> mars_phobos_velocities;
  for (Instant t = J2000; t < J2000 + 30 * Day; t += 5 * Minute) {
    mars_positions.push_back(mars_trajectory.EvaluatePosition(t));
    mars_velocities.push_back(mars_trajectory.EvaluateVelocity(t));
    phobos_positions.push_back(phobos_trajectory.EvaluatePosition(t));
    phobos_velocities.push_back(phobos_trajectory.EvaluateVelocity(t));
    mars_phobos_displacements.push_back(
        phobos_positions.back() - mars_positions.back());
    mars_phobos_velocities.push_back(phobos_velocities.back() -
                                     mars_velocities.back());
  }

  OFStream file(TEMP_DIR / "phobos.generated.wl");
  file << mathematica::Assign("ppaMarsPhobosDisplacements",
                              mars_phobos_displacements);
  file << mathematica::Assign("ppaMarsPhobosVelocities",
                              mars_phobos_velocities);
}

#endif

struct ConvergenceTestParameters {
  FixedStepSizeIntegrator<
      Ephemeris<ICRFJ2000Equator>::NewtonianMotionEquation> const& integrator;
  int iterations;
  double first_step_in_seconds;
};

class SolarSystemDynamicsConvergenceTest
    : public ::testing::TestWithParam<ConvergenceTestParameters> {
 public:
  static void SetUpTestCase() {
    file_ = OFStream(SOLUTION_DIR / "mathematica" /
                     "solar_system_convergence.generated.wl");
  }

  static void TearDownTestCase() {
    file_.~OFStream();
  }

 protected:
  FixedStepSizeIntegrator<
      Ephemeris<ICRFJ2000Equator>::NewtonianMotionEquation> const&
  integrator() const {
    return GetParam().integrator;
  }

  int iterations() const {
    return GetParam().iterations;
  }

  int first_step_in_seconds() const {
    return GetParam().first_step_in_seconds;
  }

  static OFStream file_;
};

OFStream SolarSystemDynamicsConvergenceTest::file_;

// This takes 7-8 minutes to run.
TEST_P(SolarSystemDynamicsConvergenceTest, DISABLED_Convergence) {
  google::LogToStderr();
  Time const integration_duration = 1 * JulianYear;

  SolarSystem<ICRFJ2000Equator> solar_system_at_j2000(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2451545_000000000.proto.txt");
  std::vector<std::string> const& body_names = solar_system_at_j2000.names();

  std::map<std::string, std::vector<DegreesOfFreedom<ICRFJ2000Equator>>>
      name_to_degrees_of_freedom;
  std::vector<Time> steps;
  std::vector<std::chrono::duration<double>> durations;
  for (int i = 0; i < iterations(); ++i) {
    Time const step = first_step_in_seconds() * (1 << i) * Second;
    steps.push_back(step);
    LOG(INFO) << "Integrating with step " << step;

    auto const start = std::chrono::system_clock::now();
    auto const ephemeris = solar_system_at_j2000.MakeEphemeris(
        /*fitting_tolerance=*/5 * Milli(Metre),
        Ephemeris<ICRFJ2000Equator>::FixedStepParameters(integrator(), step));
    ephemeris->Prolong(solar_system_at_j2000.epoch() + integration_duration);
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

  std::map<std::string, std::vector<RelativeDegreesOfFreedom<ICRFJ2000Equator>>>
      name_to_errors;
  for (auto const& pair : name_to_degrees_of_freedom) {
    auto const& name = pair.first;
    auto const& degrees_of_freedom = pair.second;
    CHECK_EQ(degrees_of_freedom.size(), iterations());
    for (int i = 1; i < iterations(); ++i) {
      name_to_errors[name].emplace_back(degrees_of_freedom[i] -
                                        degrees_of_freedom[0]);
    }
  }

  using MathematicaEntry = std::tuple<Time, Length, std::string, Time>;
  using MathematicaEntries = std::vector<MathematicaEntry>;

  std::vector<Length> position_errors(iterations() - 1);
  std::vector<std::string> worst_body(iterations() - 1);
  MathematicaEntries mathematica_entries;
  for (int i = 0; i < iterations() - 1; ++i) {
    for (auto const& pair : name_to_errors) {
      auto const& name = pair.first;
      auto const& errors = pair.second;
      if (position_errors[i] < errors[i].displacement().Norm()) {
        worst_body[i] = name;
      }
      position_errors[i] = std::max(position_errors[i],
                                    errors[i].displacement().Norm());
    }
    mathematica_entries.push_back({steps[i + 1],
                                   position_errors[i],
                                   mathematica::Escape(worst_body[i]),
                                   durations[i + 1].count() * Second});
  }

  std::string const test_name(
      ::testing::UnitTest::GetInstance()->current_test_info()->name());
  file_ << mathematica::Assign(std::string("ppaSolarSystemConvergence") +
                                   test_name[test_name.size() - 1],
                               mathematica::ToMathematica(mathematica_entries));
}

INSTANTIATE_TEST_CASE_P(
    AllSolarSystemDynamicsConvergenceTests,
    SolarSystemDynamicsConvergenceTest,
    ::testing::Values(
        ConvergenceTestParameters{
            SymplecticRungeKuttaNyströmIntegrator<BlanesMoan2002SRKN11B,
                                                  Position<ICRFJ2000Equator>>(),
            /*iterations=*/8,
            /*first_step_in_seconds=*/64},
        ConvergenceTestParameters{
            SymplecticRungeKuttaNyströmIntegrator<
                 McLachlanAtela1992Order5Optimal,
                 Position<ICRFJ2000Equator>>(),
             /*iterations=*/8,
             /*first_step_in_seconds=*/32},
        ConvergenceTestParameters{
            SymmetricLinearMultistepIntegrator<Quinlan1999Order8A,
                                               Position<ICRFJ2000Equator>>(),
            /*iterations=*/6,
            /*first_step_in_seconds=*/64},
        ConvergenceTestParameters{
            SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order8,
                                               Position<ICRFJ2000Equator>>(),
            /*iterations=*/6,
            /*first_step_in_seconds=*/64},
        ConvergenceTestParameters{
            SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order10,
                                               Position<ICRFJ2000Equator>>(),
            /*iterations=*/6,
            /*first_step_in_seconds=*/64},

        // This is our favorite integrator.  For a step of 600 s it gives a
        // position error of about 2.3 m on Miranda and takes about 2.0 s of
        // elapsed time.  For steps larger than about 680 s, the errors explode.
        ConvergenceTestParameters{
            SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                               Position<ICRFJ2000Equator>>(),
            /*iterations=*/5,
            /*first_step_in_seconds=*/75}));

}  // namespace astronomy
}  // namespace principia
