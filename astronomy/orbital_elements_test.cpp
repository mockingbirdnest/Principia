#include "astronomy/orbital_elements.hpp"

#include <limits>
#include <string>
#include <vector>

#include "astronomy/frames.hpp"
#include "base/file.hpp"
#include "base/not_null.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "mathematica/mathematica.hpp"
#include "physics/body_centred_non_rotating_dynamic_frame.hpp"
#include "physics/ephemeris.hpp"
#include "physics/kepler_orbit.hpp"
#include "physics/solar_system.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {

using astronomy::J2000;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::IsNear;
using testing_utilities::RelativeError;
using testing_utilities::operator""_⑴;

namespace mathematica {
namespace internal_mathematica {

std::string ToMathematica(
    astronomy::OrbitalElements::EquinoctialElements const& elements) {
  return ToMathematica(std::make_tuple((elements.t - J2000) / Second,
                                       elements.a / Metre,
                                       elements.h,
                                       elements.k,
                                       elements.λ / Radian,
                                       elements.p,
                                       elements.q,
                                       elements.pʹ,
                                       elements.qʹ));
}
}  // namespace internal_mathematica
}  // namespace mathematica

namespace astronomy {

using base::make_not_null_unique;
using base::not_null;
using base::OFStream;
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
using physics::OblateBody;
using physics::RotatingBody;
using physics::SolarSystem;
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::Cos;
using quantities::Pow;
using quantities::Sqrt;
using quantities::Time;
using quantities::astronomy::JulianYear;
using quantities::si::ArcMinute;
using quantities::si::ArcSecond;
using quantities::si::Day;
using quantities::si::Degree;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Micro;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AbsoluteError;
using testing_utilities::AlmostEquals;
using testing_utilities::IsNear;
using testing_utilities::IsOk;

class OrbitalElementsTest : public ::testing::Test {
 protected:
  OrbitalElementsTest() {}

  // Completes |initial_osculating_elements| and returns a GCRS trajectory
  // obtained by flowing the corresponding initial conditions in |ephemeris|.
  static not_null<std::unique_ptr<DiscreteTrajectory<GCRS>>>
  EarthCentredTrajectory(
      KeplerianElements<GCRS>& initial_osculating_elements,
      Instant const& initial_time,
      Instant const& final_time,
      Ephemeris<ICRS>& ephemeris) {
    MassiveBody const& earth = FindEarthOrDie(ephemeris);
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
            /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
            /*length_integration_tolerance=*/1 * Milli(Metre),
            /*speed_integration_tolerance=*/1 * Milli(Metre) / Second
        },
        /*max_ephemeris_steps=*/std::numeric_limits<std::int64_t>::max());
    auto result = make_not_null_unique<DiscreteTrajectory<GCRS>>();
    for (auto it = icrs_trajectory.Begin(); it != icrs_trajectory.End(); ++it) {
      result->Append(
          it.time(),
          gcrs.ToThisFrameAtTime(it.time())(it.degrees_of_freedom()));
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
          SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                             Position<ICRS>>(),
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
      MasslessBody{});
  ASSERT_THAT(status_or_elements, IsOk());
  OrbitalElements const& elements = status_or_elements.ValueOrDie();
  EXPECT_THAT(
      AbsoluteError(*initial_osculating.period, elements.anomalistic_period()),
      IsNear(243_⑴ * Micro(Second)));
  EXPECT_THAT(
      AbsoluteError(*initial_osculating.period, elements.nodal_period()),
      IsNear(3.7_⑴ * Milli(Second)));
  EXPECT_THAT(
      AbsoluteError(*initial_osculating.period, elements.sidereal_period()),
      IsNear(1.8_⑴ * Micro(Second)));

  EXPECT_THAT(elements.nodal_precession(), IsNear(1.2_⑴ * Degree / JulianYear));

  // Mean element values.
  EXPECT_THAT(AbsoluteError(*initial_osculating.semimajor_axis,
                            elements.mean_semimajor_axis_interval().midpoint()),
              IsNear(330_⑴ * Micro(Metre)));
  EXPECT_THAT(AbsoluteError(*initial_osculating.eccentricity,
                            elements.mean_eccentricity_interval().midpoint()),
              IsNear(4.5e-11_⑴));
  EXPECT_THAT(AbsoluteError(initial_osculating.inclination,
                            elements.mean_inclination_interval().midpoint()),
              IsNear(0.19_⑴ * Micro(ArcSecond)));
  EXPECT_THAT(
      AbsoluteError(
          initial_osculating.longitude_of_ascending_node,
          elements.mean_longitude_of_ascending_node_interval().midpoint()),
      IsNear(57_⑴ * ArcSecond));
  EXPECT_THAT(
      AbsoluteError(*initial_osculating.argument_of_periapsis,
                    elements.mean_argument_of_periapsis_interval().midpoint()),
      IsNear(49_⑴ * ArcSecond));

  // Mean element stability.
  EXPECT_THAT(elements.mean_semimajor_axis_interval().measure(),
              IsNear(0.8_⑴ * Milli(Metre)));
  EXPECT_THAT(elements.mean_eccentricity_interval().measure(),
              IsNear(1.0e-10_⑴));
  EXPECT_THAT(elements.mean_inclination_interval().measure(),
              IsNear(0.66_⑴ * Micro(ArcSecond)));
  EXPECT_THAT(elements.mean_longitude_of_ascending_node_interval().measure(),
              IsNear(2.1_⑴ * ArcMinute));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().measure(),
              IsNear(1.7_⑴ * ArcMinute));

  OFStream f(SOLUTION_DIR / "mathematica" /
             "unperturbed_elements.generated.wl");
  f << mathematica::Assign("unperturbedOsculating",
                           elements.osculating_equinoctial_elements());
  f << mathematica::Assign("unperturbedMean",
                           elements.mean_equinoctial_elements());
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
          SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                             Position<ICRS>>(),
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
      MasslessBody{});
  ASSERT_THAT(status_or_elements, IsOk());
  OrbitalElements const& elements = status_or_elements.ValueOrDie();
  EXPECT_THAT(elements.anomalistic_period() - *initial_osculating.period,
              IsNear(-7.8_⑴ * Second));
  EXPECT_THAT(elements.nodal_period() - *initial_osculating.period,
              IsNear(-23_⑴ * Second));
  EXPECT_THAT(elements.sidereal_period() - *initial_osculating.period,
              IsNear(-16_⑴ * Second));

  // The notation for the computation of the theoretical precessions follows
  // Capderou (2012), Satellites : de Kepler au GPS, section 7.1.1.
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

  EXPECT_THAT(theoretical_Ωʹ, IsNear(-7.2_⑴ * Degree / Day));
  EXPECT_THAT(theoretical_ωʹ, IsNear(14_⑴ * Degree / Day));

  EXPECT_THAT(RelativeError(theoretical_Ωʹ, elements.nodal_precession()),
              IsNear(0.0028_⑴));

  // Mean element values.  Since Ω and ω precess rapidly, the midpoint of the
  // range of values is of no interest.
  EXPECT_THAT(AbsoluteError(*initial_osculating.semimajor_axis,
                            elements.mean_semimajor_axis_interval().midpoint()),
              IsNear(25_⑴ * Metre));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.0013_⑴));
  EXPECT_THAT(AbsoluteError(initial_osculating.inclination,
                            elements.mean_inclination_interval().midpoint()),
              IsNear(1.4_⑴ * Micro(ArcSecond)));

  // Mean element stability: Ω and ω precess as expected, the other elements are
  // stable.
  EXPECT_THAT(elements.mean_semimajor_axis_interval().measure(),
              IsNear(70_⑴ * Milli(Metre)));
  EXPECT_THAT(elements.mean_eccentricity_interval().measure(),
              IsNear(3.8e-9_⑴));
  EXPECT_THAT(elements.mean_inclination_interval().measure(),
              IsNear(0.8_⑴ * Micro(ArcSecond)));
  EXPECT_THAT(
      RelativeError(
          -theoretical_Ωʹ * mission_duration,
          elements.mean_longitude_of_ascending_node_interval().measure()),
      IsNear(0.0039_⑴));
  EXPECT_THAT(
      RelativeError(theoretical_ωʹ * mission_duration,
                    elements.mean_argument_of_periapsis_interval().measure()),
      IsNear(0.0029_⑴));

  OFStream f(SOLUTION_DIR / "mathematica" /
             "j2_perturbed_elements.generated.wl");
  f << mathematica::Assign("j2PerturbedOsculating",
                           elements.osculating_equinoctial_elements());
  f << mathematica::Assign("j2PerturbedMean",
                           elements.mean_equinoctial_elements());
}

TEST_F(OrbitalElementsTest, RealPerturbation) {
  SolarSystem<ICRS> solar_system(
      SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
      SOLUTION_DIR / "astronomy" /
          "sol_initial_state_jd_2451545_000000000.proto.txt");
  auto const ephemeris = solar_system.MakeEphemeris(
      /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                               /*geopotential_tolerance=*/0x1p-24},
      Ephemeris<ICRS>::FixedStepParameters(
          SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                             Position<ICRS>>(),
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
      MasslessBody{});
  ASSERT_THAT(status_or_elements, IsOk());
  OrbitalElements const& elements = status_or_elements.ValueOrDie();
  EXPECT_THAT(elements.anomalistic_period() - *initial_osculating.period,
              IsNear(-8.0_⑴ * Second));
  EXPECT_THAT(elements.nodal_period() - *initial_osculating.period,
              IsNear(-14_⑴ * Second));
  EXPECT_THAT(elements.sidereal_period() - *initial_osculating.period,
              IsNear(-16_⑴ * Second));

  // This value is meaningless, see below.
  EXPECT_THAT(elements.nodal_precession(), IsNear(2.0_⑴ * Degree / Day));

  // Mean element values.
  EXPECT_THAT(AbsoluteError(*initial_osculating.semimajor_axis,
                            elements.mean_semimajor_axis_interval().midpoint()),
              IsNear(105_⑴ * Metre));
  EXPECT_THAT(elements.mean_eccentricity_interval().midpoint(),
              IsNear(0.0014_⑴));
  EXPECT_THAT(AbsoluteError(initial_osculating.inclination,
                            elements.mean_inclination_interval().midpoint()),
              IsNear(6.0_⑴ * ArcSecond));

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
              IsNear(20_⑴ * Metre));
  EXPECT_THAT(elements.mean_eccentricity_interval().measure(),
              IsNear(9.2e-5_⑴));
  EXPECT_THAT(elements.mean_inclination_interval().measure(),
              IsNear(11_⑴ * ArcSecond));
  EXPECT_THAT(elements.mean_longitude_of_ascending_node_interval().measure(),
              IsNear(136_⑴ * Degree));
  EXPECT_THAT(elements.mean_argument_of_periapsis_interval().measure(),
              IsNear(154_⑴ * Degree));

  OFStream f(SOLUTION_DIR / "mathematica" /
             "fully_perturbed_elements.generated.wl");
  f << mathematica::Assign("fullyPerturbedOsculating",
                           elements.osculating_equinoctial_elements());
  f << mathematica::Assign("fullyPerturbedMean",
                           elements.mean_equinoctial_elements());
}

#endif

}  // namespace astronomy
}  // namespace principia
