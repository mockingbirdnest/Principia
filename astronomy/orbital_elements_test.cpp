#include "astronomy/orbital_elements.hpp"

#include <memory>
#include <string>
#include <vector>

#include "astronomy/epoch.hpp"
#include "astronomy/frames.hpp"
#include "astronomy/time_scales.hpp"
#include "base/not_null.hpp"
#include "geometry/instant.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/symmetric_linear_multistep_integrator.hpp"
#include "mathematica/logger.hpp"
#include "mathematica/mathematica.hpp"
#include "physics/body_centred_non_rotating_reference_frame.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "physics/discrete_trajectory.hpp"
#include "physics/ephemeris.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/massive_body.hpp"
#include "physics/massless_body.hpp"
#include "physics/solar_system.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"  // 🧙 For EXPECT_OK.
#include "testing_utilities/numerics_matchers.hpp"

namespace principia {
namespace astronomy {

using ::testing::AnyOf;
using ::testing::Lt;
using namespace principia::astronomy::_epoch;
using namespace principia::astronomy::_frames;
using namespace principia::astronomy::_orbital_elements;
using namespace principia::astronomy::_time_scales;
using namespace principia::base::_not_null;
using namespace principia::geometry::_instant;
using namespace principia::integrators::_methods;
using namespace principia::integrators::_symmetric_linear_multistep_integrator;
using namespace principia::mathematica::_logger;
using namespace principia::mathematica::_mathematica;
using namespace principia::physics::_body_centred_non_rotating_reference_frame;
using namespace principia::physics::_degrees_of_freedom;
using namespace principia::physics::_discrete_trajectory;
using namespace principia::physics::_ephemeris;
using namespace principia::physics::_kepler_orbit;
using namespace principia::physics::_massive_body;
using namespace principia::physics::_massless_body;
using namespace principia::physics::_solar_system;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_approximate_quantity;
using namespace principia::testing_utilities::_is_near;
using namespace principia::testing_utilities::_matchers;
using namespace principia::testing_utilities::_numerics_matchers;

class OrbitalElementsTest : public ::testing::Test {
 protected:
  OrbitalElementsTest() {}

  // Completes `initial_osculating_elements` and returns a GCRS trajectory
  // obtained by flowing the corresponding initial conditions in `ephemeris`.
  static not_null<std::unique_ptr<DiscreteTrajectory<GCRS>>>
  EarthCentredTrajectory(
      KeplerianElements<GCRS>& initial_osculating_elements,
      Instant const& initial_time,
      Instant const& final_time,
      Ephemeris<ICRS>& ephemeris) {
    MassiveBody const& earth = FindEarthOrDie(ephemeris);
    EXPECT_OK(ephemeris.Prolong(final_time));
    BodyCentredNonRotatingReferenceFrame<ICRS, GCRS> gcrs{&ephemeris, &earth};
    DiscreteTrajectory<ICRS> icrs_trajectory;
    KeplerOrbit<GCRS> initial_osculating_orbit{earth,
                                               MasslessBody{},
                                               initial_osculating_elements,
                                               initial_time};
    initial_osculating_elements = initial_osculating_orbit.elements_at_epoch();
    icrs_trajectory.segments().front().SetDownsampling(
        {.max_dense_intervals = 10'000,
         .tolerance = 1 * Milli(Metre)});
    EXPECT_OK(icrs_trajectory.Append(
        initial_time,
        gcrs.FromThisFrameAtTime(initial_time)(
            DegreesOfFreedom<GCRS>{GCRS::origin, GCRS::unmoving} +
            initial_osculating_orbit.StateVectors(initial_time))));
    auto instance = ephemeris.NewInstance(
        {&icrs_trajectory},
        Ephemeris<ICRS>::NoIntrinsicAccelerations,
        Ephemeris<ICRS>::FixedStepParameters(
            SymmetricLinearMultistepIntegrator<
                Quinlan1999Order8A,
                Ephemeris<ICRS>::NewtonianMotionEquation>(),
            /*step=*/10 * Second));
    CHECK_OK(ephemeris.FlowWithFixedStep(final_time, *instance));
    auto result = make_not_null_unique<DiscreteTrajectory<GCRS>>();
    for (auto const& [time, degrees_of_freedom] : icrs_trajectory) {
      EXPECT_OK(result->Append(
          time, gcrs.ToThisFrameAtTime(time)(degrees_of_freedom)));
    }
    return result;
  }

 private:
  static MassiveBody const& FindEarthOrDie(Ephemeris<ICRS> const& ephemeris) {
    for (not_null<MassiveBody const*> const body : ephemeris.bodies()) {
      if (body->name() == "Earth") {
        return *body;
      }
    }
    LOG(FATAL) << "Ephemeris has no Earth";
  }
};

#if !defined(_DEBUG)

// TODO(phl): Fix these tests, maybe by using a tolerance-to-error-ratio
// function that considers multiple parameters.
#if 0
TEST_F(OrbitalElementsTest, KeplerOrbit) {
  // The satellite is under the influence of an isotropic Earth and no third
  // bodies.
  SolarSystem<ICRS> solar_system(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2451545_000000000.proto.txt");
  std::vector<std::string> const names = solar_system.names();
  for (auto const& name : names) {
    if (name != "Earth") {
      solar_system.RemoveMassiveBody(name);
    }
  }
  solar_system.LimitOblatenessToDegree("Earth", 0);
  auto const ephemeris = solar_system.MakeEphemeris(
      /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<ICRS>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<
              QuinlanTremaine1990Order12,
              Ephemeris<ICRS>::NewtonianMotionEquation>(),
          /*step=*/1 * JulianYear));
  MassiveBody const& spherical_earth =
      *solar_system.massive_body(*ephemeris, "Earth");

  KeplerianElements<GCRS> initial_osculating;
  initial_osculating.semimajor_axis = 7000 * Kilo(Metre);
  initial_osculating.eccentricity = 1e-6;
  initial_osculating.inclination = 10 * Milli(ArcSecond);
  initial_osculating.longitude_of_ascending_node = 10 * Degree;
  initial_osculating.argument_of_periapsis = 20 * Degree;
  initial_osculating.mean_anomaly = 30 * Degree;
  auto const status_or_elements = OrbitalElements::ForTrajectory(
      *EarthCentredTrajectory(initial_osculating,
                              J2000,
                              J2000 + 10 * Day, *ephemeris),
      spherical_earth,
      MasslessBody{},
      /*fill_osculating_equinoctial_elements=*/true);
  ASSERT_THAT(status_or_elements, IsOk());
  OrbitalElements const& elements = status_or_elements.value();
  EXPECT_THAT(
      elements.anomalistic_period(),
      AbsoluteErrorFrom(*initial_osculating.period, Lt(510 * Micro(Second))));
  EXPECT_THAT(
      elements.nodal_period(),
      AbsoluteErrorFrom(*initial_osculating.period, Lt(4.2 * Milli(Second))));
  EXPECT_THAT(
      elements.sidereal_period(),
      AbsoluteErrorFrom(*initial_osculating.period, Lt(1.9 * Micro(Second))));

  EXPECT_THAT(elements.nodal_precession(), Lt(1.4 * Degree / JulianYear));

  // Mean element values.
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              AbsoluteErrorFrom(*initial_osculating.semimajor_axis,
                                Lt(410 * Micro(Metre))));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              AbsoluteErrorFrom(*initial_osculating.eccentricity,
                                Lt(4.6e-11)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              AbsoluteErrorFrom(initial_osculating.inclination,
                                Lt(0.64 * Micro(ArcSecond))));
  EXPECT_THAT(elements.mean_longitude_of_ascending_node_interval().midpoint(),
              AbsoluteErrorFrom(initial_osculating.longitude_of_ascending_node,
                                Lt(66 * ArcSecond)));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().midpoint(),
              AbsoluteErrorFrom(*initial_osculating.argument_of_periapsis,
                                Lt(74 * ArcSecond)));

  // Mean element stability.
  EXPECT_THAT(elements.mean_semimajor_axis_interval().measure(),
              Lt(1.1 * Milli(Metre)));
  EXPECT_THAT(elements.mean_eccentricity_interval().measure(),
              Lt(1.1e-10));
  EXPECT_THAT(elements.mean_inclination_interval().measure(),
              Lt(1.4 * Micro(ArcSecond)));
  EXPECT_THAT(elements.mean_longitude_of_ascending_node_interval().measure(),
              Lt(2.3 * ArcMinute));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().measure(),
              Lt(2.4 * ArcMinute));

  Logger logger(
      SOLUTION_DIR / "mathematica" / "unperturbed_elements.generated.wl",
      /*make_unique=*/false);
  logger.Set("unperturbedOsculating",
             elements.osculating_equinoctial_elements(),
             ExpressIn(Metre, Second, Radian));
  logger.Set("unperturbedMean",
             elements.mean_equinoctial_elements(),
             ExpressIn(Metre, Second, Radian));
}

TEST_F(OrbitalElementsTest, J2Perturbation) {
  // The satellite is under the influence of an Earth with a zonal geopotential
  // of degree 2 and no third bodies.

  SolarSystem<ICRS> solar_system(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2451545_000000000.proto.txt");
  std::vector<std::string> const names = solar_system.names();
  for (auto const& name : names) {
    if (name != "Earth") {
      solar_system.RemoveMassiveBody(name);
    }
  }
  solar_system.LimitOblatenessToDegree("Earth", 2);
  solar_system.LimitOblatenessToZonal("Earth");
  auto const ephemeris = solar_system.MakeEphemeris(
      /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<ICRS>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<
              QuinlanTremaine1990Order12,
              Ephemeris<ICRS>::NewtonianMotionEquation>(),
          /*step=*/1 * JulianYear));
  auto const& oblate_earth = dynamic_cast<OblateBody<ICRS> const&>(
      *solar_system.massive_body(*ephemeris, "Earth"));

  Time const mission_duration = 10 * Day;

  KeplerianElements<GCRS> initial_osculating;
  initial_osculating.semimajor_axis = 7000 * Kilo(Metre);
  initial_osculating.eccentricity = 1e-6;
  initial_osculating.inclination = 10 * Milli(ArcSecond);
  initial_osculating.longitude_of_ascending_node = 10 * Degree;
  initial_osculating.argument_of_periapsis = 20 * Degree;
  initial_osculating.mean_anomaly = 30 * Degree;
  auto const status_or_elements = OrbitalElements::ForTrajectory(
      *EarthCentredTrajectory(initial_osculating,
                              J2000,
                              J2000 + mission_duration,
                              *ephemeris),
      oblate_earth,
      MasslessBody{},
      /*fill_osculating_equinoctial_elements=*/true);
  ASSERT_THAT(status_or_elements, IsOk());
  OrbitalElements const& elements = status_or_elements.value();
  EXPECT_THAT(
      elements.anomalistic_period(),
      DifferenceFrom(*initial_osculating.period, IsNear(-7.8_(1) * Second)));
  EXPECT_THAT(
      elements.nodal_period(),
      DifferenceFrom(*initial_osculating.period, IsNear(-23_(1) * Second)));
  EXPECT_THAT(
      elements.sidereal_period(),
      DifferenceFrom(*initial_osculating.period, IsNear(-16_(1) * Second)));

  // The notation for the computation of the theoretical precessions follows
  // [Cap12], section 7.1.1.
  double const η = elements.mean_semimajor_axis_interval().midpoint() /
                   oblate_earth.reference_radius();
  Angle const i = elements.mean_inclination_interval().midpoint();
  AngularFrequency const K0 = 3.0 / 2.0 * oblate_earth.j2() *
                              Sqrt(oblate_earth.gravitational_parameter() /
                                   Pow<3>(oblate_earth.reference_radius())) *
                              Radian;
  // See (7.6).
  AngularFrequency const theoretical_Ωʹ = -K0 * Pow<-7>(Sqrt(η)) * Cos(i);
  // See (7.13).
  AngularFrequency const theoretical_ωʹ =
      0.5 * K0 * Pow<-7>(Sqrt(η)) * (5 * Pow<2>(Cos(i)) - 1);

  EXPECT_THAT(theoretical_Ωʹ, IsNear(-7.2_(1) * Degree / Day));
  EXPECT_THAT(theoretical_ωʹ, IsNear(14_(1) * Degree / Day));

  EXPECT_THAT(RelativeError(theoretical_Ωʹ, elements.nodal_precession()),
              IsNear(0.0029_(1)));

  // Mean element values.  Since Ω and ω precess rapidly, the midpoint of the
  // range of values is of no interest.
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              AbsoluteErrorFrom(*initial_osculating.semimajor_axis,
                                IsNear(25_(1) * Metre)));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.0013_(1)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              AbsoluteErrorFrom(initial_osculating.inclination,
                                Lt(2.0 * Micro(ArcSecond))));

  // Mean element stability: Ω and ω precess as expected, the other elements are
  // stable.
  EXPECT_THAT(elements.mean_semimajor_axis_interval().measure(),
              IsNear(70_(1) * Milli(Metre)));
  EXPECT_THAT(elements.mean_eccentricity_interval().measure(),
              IsNear(3.8e-9_(1)));
  EXPECT_THAT(elements.mean_inclination_interval().measure(),
              Lt(1.1 * Micro(ArcSecond)));
  EXPECT_THAT(
      RelativeError(
          -theoretical_Ωʹ * mission_duration,
          elements.mean_longitude_of_ascending_node_interval().measure()),
      IsNear(0.004_(1)));
  EXPECT_THAT(
      RelativeError(theoretical_ωʹ * mission_duration,
                    elements.mean_argument_of_periapsis_interval().measure()),
      IsNear(0.0029_(1)));

  Logger logger(
      SOLUTION_DIR / "mathematica" / "j2_perturbed_elements.generated.wl",
      /*make_unique=*/false);
  logger.Set("j2PerturbedOsculating",
             elements.osculating_equinoctial_elements(),
             ExpressIn(Metre, Second, Radian));
  logger.Set("j2PerturbedMean",
             elements.mean_equinoctial_elements(),
             ExpressIn(Metre, Second, Radian));
}
#endif

TEST_F(OrbitalElementsTest, RealPerturbation) {
  SolarSystem<ICRS> solar_system(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2451545_000000000.proto.txt");
  auto const ephemeris = solar_system.MakeEphemeris(
      /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<ICRS>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<
              QuinlanTremaine1990Order12,
              Ephemeris<ICRS>::NewtonianMotionEquation>(),
          /*step=*/10 * Minute));
  MassiveBody const& earth = *solar_system.massive_body(*ephemeris, "Earth");

  Time const mission_duration = 10 * Day;

  KeplerianElements<GCRS> initial_osculating;
  initial_osculating.semimajor_axis = 7000 * Kilo(Metre);
  initial_osculating.eccentricity = 1e-6;
  initial_osculating.inclination = 10 * Milli(ArcSecond);
  initial_osculating.longitude_of_ascending_node = 10 * Degree;
  initial_osculating.argument_of_periapsis = 20 * Degree;
  initial_osculating.mean_anomaly = 30 * Degree;
  auto const status_or_elements = OrbitalElements::ForTrajectory(
      *EarthCentredTrajectory(
          initial_osculating, J2000, J2000 + mission_duration, *ephemeris),
      earth,
      MasslessBody{},
      /*fill_osculating_equinoctial_elements=*/true);
  ASSERT_THAT(status_or_elements, IsOk());
  OrbitalElements const& elements = status_or_elements.value();
  EXPECT_THAT(
      elements.anomalistic_period(),
      DifferenceFrom(*initial_osculating.period, IsNear(-8.0_(1) * Second)));
  EXPECT_THAT(
      elements.nodal_period(),
      DifferenceFrom(*initial_osculating.period, IsNear(-14_(1) * Second)));
  EXPECT_THAT(
      elements.sidereal_period(),
      DifferenceFrom(*initial_osculating.period, IsNear(-16_(1) * Second)));

  // This value is meaningless, see below.
  EXPECT_THAT(elements.nodal_precession(), IsNear(2.0_(1) * Degree / Day));

  // Mean element values.
  EXPECT_THAT(elements.mean_semimajor_axis_interval().midpoint(),
              AbsoluteErrorFrom(*initial_osculating.semimajor_axis,
                                IsNear(105_(1) * Metre)));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.0014_(1)));
  EXPECT_THAT(elements.mean_inclination_interval().midpoint(),
              AbsoluteErrorFrom(initial_osculating.inclination,
                                IsNear(6.0_(1) * ArcSecond)));

  // Mean element stability: Ω and ω exhibit a daily oscillation (likely due to
  // the tesseral terms of the geopotential) as the very low inclination means
  // that these elements are singular.
  // The other elements are stable.
  // A closer analysis would show that the longitude of periapsis ϖ = Ω + ω
  // exhibits a precession that is largely free of oscillations, and that, if
  // its oscillations are filtered, the argument of periapsis ω precesses as
  // expected; the longitude of the ascending node Ω exhibits no obvious
  // precession even if its daily oscillation is filtered out.
  EXPECT_THAT(elements.mean_semimajor_axis_interval().measure(),
              IsNear(20_(1) * Metre));
  EXPECT_THAT(elements.mean_eccentricity_interval().measure(),
              IsNear(1.0e-4_(1)));
  EXPECT_THAT(elements.mean_inclination_interval().measure(),
              IsNear(11_(1) * ArcSecond));
  EXPECT_THAT(elements.mean_longitude_of_ascending_node_interval().measure(),
              IsNear(136_(1) * Degree));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().measure(),
              IsNear(154_(1) * Degree));

  Logger logger(
      SOLUTION_DIR / "mathematica" / "fully_perturbed_elements.generated.wl",
      /*make_unique=*/false);
  logger.Set("fullyPerturbedOsculating",
             elements.osculating_equinoctial_elements(),
             ExpressIn(Metre, Second, Radian));
  logger.Set("fullyPerturbedMean",
             elements.mean_equinoctial_elements(),
             ExpressIn(Metre, Second, Radian));
}

TEST_F(OrbitalElementsTest, Escape) {
  SolarSystem<ICRS> solar_system(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2451545_000000000.proto.txt");
  auto const ephemeris = solar_system.MakeEphemeris(
      /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<ICRS>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<
              QuinlanTremaine1990Order12,
              Ephemeris<ICRS>::NewtonianMotionEquation>(),
          /*step=*/10 * Minute));
  MassiveBody const& earth = *solar_system.massive_body(*ephemeris, "Earth");

  Time const mission_duration = 10 * Day;

  KeplerianElements<GCRS> initial_osculating;
  initial_osculating.periapsis_distance = 7000 * Kilo(Metre);
  initial_osculating.eccentricity = 1.2;
  initial_osculating.inclination = 10 * Milli(ArcSecond);
  initial_osculating.longitude_of_ascending_node = 10 * Degree;
  initial_osculating.argument_of_periapsis = 20 * Degree;
  initial_osculating.hyperbolic_mean_anomaly = 30 * Degree;
  EXPECT_THAT(
      OrbitalElements::ForTrajectory(
          *EarthCentredTrajectory(
              initial_osculating, J2000, J2000 + mission_duration, *ephemeris),
          earth,
          MasslessBody{},
          /*fill_osculating_equinoctial_elements=*/true)
          .status(),
      StatusIs(absl::StatusCode::kOutOfRange));
}

TEST_F(OrbitalElementsTest, Years) {
  SolarSystem<ICRS> solar_system(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2451545_000000000.proto.txt");
  auto const ephemeris = solar_system.MakeEphemeris(
      /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<ICRS>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<
              QuinlanTremaine1990Order12,
              Ephemeris<ICRS>::NewtonianMotionEquation>(),
          /*step=*/10 * Minute));
  MassiveBody const& sun = *solar_system.massive_body(*ephemeris, "Sun");
  MassiveBody const& earth = *solar_system.massive_body(*ephemeris, "Earth");

  LOG(ERROR) << "Prolonging...";

  ASSERT_THAT(ephemeris->Prolong("2050-01-01T12:00:00"_TT), IsOk());

  LOG(ERROR) << "Analysing...";

  auto const status_or_elements =
      OrbitalElements::OrbitalElements::ForTrajectory(
          *ephemeris->trajectory(&sun),
          BodyCentredNonRotatingReferenceFrame<ICRS, GCRS>(ephemeris.get(),
                                                           &earth),
          /*primary=*/earth,
          /*secondary=*/sun);
  LOG(ERROR) << "Done.";
  ASSERT_THAT(status_or_elements, IsOk());
  OrbitalElements const& elements = status_or_elements.value();
  // The value given by the Astronomical Almanac is 365.259'636.
  EXPECT_THAT(elements.anomalistic_period(), IsNear(365.259'09_(1) * Day));
  // This should be 365.242190 if the node were with respect to the precessing
  // equator, but we do not have axial precession, so it is just the sidereal
  // period up to the precession of the ecliptic.
  EXPECT_THAT(elements.nodal_period(), IsNear(365.256'350_(1) * Day));
  // https://hpiers.obspm.fr/eop-pc/models/constants.html gives 365.256'363'004.
  EXPECT_THAT(elements.sidereal_period(), IsNear(365.256'836_(1) * Day));
}

#endif

}  // namespace astronomy
}  // namespace principia
