#include "orbital_elements.hpp"

#include "astronomy/frames.hpp"
#include "base/not_null.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/body_centred_non_rotating_dynamic_frame.hpp"
#include "physics/ephemeris.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/solar_system.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/is_near.hpp"

namespace principia {
namespace astronomy {

using base::make_not_null_unique;
using base::not_null;
using geometry::Instant;
using geometry::Position;
using geometry::Velocity;
using integrators::EmbeddedExplicitRungeKuttaNyströmIntegrator;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::DormandالمكاوىPrince1986RKN434FM;
using integrators::methods::QuinlanTremaine1990Order12;
using physics::BodyCentredNonRotatingDynamicFrame;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using physics::Ephemeris;
using physics::KeplerianElements;
using physics::KeplerOrbit;
using physics::MassiveBody;
using physics::MasslessBody;
using physics::RotatingBody;
using physics::SolarSystem;
using quantities::astronomy::JulianYear;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::IsNear;
using testing_utilities::IsOk;

class OrbitalElementsTest : public ::testing::Test {
 protected:
  OrbitalElementsTest()
      : solar_system_(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2451545_000000000.proto.txt"),
        real_ephemeris_(solar_system_.MakeEphemeris(
            /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                                     /*geopotential_tolerance=*/0x1p-24},
            Ephemeris<ICRS>::FixedStepParameters(
                SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                                   Position<ICRS>>(),
                /*step=*/10 * Minute))),
        real_earth_(*solar_system_.massive_body(*real_ephemeris_, "Earth")),
        oblate_earth_ephemeris_([this] {
          std::vector<std::string> const names = solar_system_.names();
          for (auto const& name : names) {
            if (name != "Earth") {
              solar_system_.RemoveMassiveBody(name);
            }
          }
          solar_system_.LimitOblatenessToDegree("Earth", 2);
          solar_system_.LimitOblatenessToZonal("Earth");
          return solar_system_.MakeEphemeris(
              /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                                       /*geopotential_tolerance=*/0x1p-24},
              Ephemeris<ICRS>::FixedStepParameters(
                  SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                                     Position<ICRS>>(),
                  /*step=*/1 * JulianYear));
        }()),
        oblate_earth_(
            *solar_system_.massive_body(*oblate_earth_ephemeris_, "Earth")),
        spherical_earth_ephemeris_([this] {
          solar_system_.LimitOblatenessToDegree("Earth", 0);
          return solar_system_.MakeEphemeris(
              /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                                       /*geopotential_tolerance=*/0x1p-24},
              Ephemeris<ICRS>::FixedStepParameters(
                  SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                                     Position<ICRS>>(),
                  /*step=*/1 * JulianYear));
        }()),
        spherical_earth_(
            *solar_system_.massive_body(*spherical_earth_ephemeris_, "Earth")) {
  }

  // Completes |initial_osculating_elements| and returns a GCRS trajectory
  // obtained by flowing the corresponding initial conditions in |ephemeris|.
  static not_null<std::unique_ptr<DiscreteTrajectory<GCRS>>> EarthCentredTrajectory(
      KeplerianElements<GCRS>& initial_osculating_elements,
      Instant const& initial_time,
      Instant const& final_time,
      Ephemeris<ICRS>& ephemeris) {
    MassiveBody const& earth = [&ephemeris]() -> MassiveBody const& {
      for (not_null<MassiveBody const*> const body : ephemeris.bodies()) {
        if (body->name() == "Earth") {
          return *body;
        }
      }
      LOG(FATAL) << "Ephemeris has no Earth";
    }();
    ephemeris.Prolong(final_time);
    BodyCentredNonRotatingDynamicFrame<ICRS, GCRS> gcrs{&ephemeris, &earth};
    DiscreteTrajectory<ICRS> icrs_trajectory;
    KeplerOrbit<GCRS> initial_osculating_orbit{earth,
                                               MasslessBody{},
                                               initial_osculating_elements,
                                               initial_time};
    initial_osculating_elements = initial_osculating_orbit.elements_at_epoch();
    icrs_trajectory.Append(
        initial_time,
        gcrs.FromThisFrameAtTime(initial_time)(
            DegreesOfFreedom<GCRS>{GCRS::origin, Velocity<GCRS>{}} +
            initial_osculating_orbit.StateVectors(initial_time)));
    ephemeris.FlowWithAdaptiveStep(
        &icrs_trajectory,
        Ephemeris<ICRS>::NoIntrinsicAcceleration,
        final_time,
        Ephemeris<ICRS>::AdaptiveStepParameters{
            EmbeddedExplicitRungeKuttaNyströmIntegrator<
                DormandالمكاوىPrince1986RKN434FM,
                Position<ICRS>>(),
            /*max_steps=*/std::numeric_limits<std::int16_t>::max(),
            /*length_integration_tolerance=*/1 * Milli(Metre),
            /*speed_integration_tolerance=*/1 * Milli(Metre) / Second
        },
        /*max_ephemeris_steps=*/std::numeric_limits<std::int16_t>::max(),
        /*last_point_only=*/false);
    auto result = make_not_null_unique<DiscreteTrajectory<GCRS>>();
    for (auto it = icrs_trajectory.Begin(); it != icrs_trajectory.End(); ++it) {
      result->Append(
          it.time(),
          gcrs.ToThisFrameAtTime(it.time())(it.degrees_of_freedom()));
    }
    return result;
  }

  SolarSystem<ICRS> solar_system_;

  // The solar system.
  not_null<std::unique_ptr<Ephemeris<ICRS>>> real_ephemeris_;
  MassiveBody const& real_earth_;

  // An oblate Earth with no other bodies.
  not_null<std::unique_ptr<Ephemeris<ICRS>>> oblate_earth_ephemeris_;
  MassiveBody const& oblate_earth_;

  // A spherical Earth with no other bodies.
  not_null<std::unique_ptr<Ephemeris<ICRS>>> spherical_earth_ephemeris_;
  MassiveBody const& spherical_earth_;
};


TEST_F(OrbitalElementsTest, CircularEquatorialKeplerOrbit) {
  KeplerianElements<GCRS> initial_osculating;
  initial_osculating.semimajor_axis = 7000 * Kilo(Metre);
  initial_osculating.eccentricity = 1e-9;
  initial_osculating.inclination = 1.0 / (60 * 60) * Degree;
  initial_osculating.longitude_of_ascending_node = 10 * Degree;
  initial_osculating.argument_of_periapsis = 20 * Degree;
  initial_osculating.mean_anomaly = 30 * Degree;
  auto const status_or_elements = OrbitalElements::ForTrajectory(
      *EarthCentredTrajectory(initial_osculating,
                              J2000,
                              J2000 + 10 * Day,
                              *spherical_earth_ephemeris_),
      spherical_earth_,
      MasslessBody{});
  ASSERT_THAT(status_or_elements, IsOk());
  OrbitalElements const& elements = status_or_elements.ValueOrDie();
  EXPECT_THAT(elements.anomalistic_period(),
              AlmostEquals(*initial_osculating.period, 0));
  EXPECT_THAT(elements.nodal_period(),
              AlmostEquals(*initial_osculating.period, 0));
  EXPECT_THAT(elements.sidereal_period(),
              AlmostEquals(*initial_osculating.period, 0));

  EXPECT_THAT(elements.nodal_precession(),
              AlmostEquals(0 * Degree / Second, 0));

  EXPECT_THAT(elements.mean_semimajor_axis().min,
              AlmostEquals(*initial_osculating.semimajor_axis, 0));
  EXPECT_THAT(elements.mean_eccentricity().min,
              AlmostEquals(*initial_osculating.eccentricity, 0));
  EXPECT_THAT(elements.mean_inclination().min,
              AlmostEquals(initial_osculating.inclination, 0));

  EXPECT_THAT(elements.mean_longitude_of_ascending_node().min,
              AlmostEquals(initial_osculating.longitude_of_ascending_node, 0));
  EXPECT_THAT(elements.mean_argument_of_periapsis().min,
              AlmostEquals(initial_osculating.longitude_of_ascending_node, 0));
}

}  // namespace astronomy
}  // namespace principia
