#include "orbital_elements.hpp"

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
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {

namespace mathematica {
namespace internal_mathematica {

using astronomy::J2000;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;

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
        oblate_earth_(dynamic_cast<OblateBody<ICRS> const&>(
            *solar_system_.massive_body(*oblate_earth_ephemeris_, "Earth"))),
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
  static not_null<std::unique_ptr<DiscreteTrajectory<GCRS>>>
  EarthCentredTrajectory(
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
            /*max_steps=*/std::numeric_limits<std::int64_t>::max(),
            /*length_integration_tolerance=*/1 * Milli(Metre),
            /*speed_integration_tolerance=*/1 * Milli(Metre) / Second
        },
        /*max_ephemeris_steps=*/std::numeric_limits<std::int64_t>::max(),
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
  OblateBody<ICRS> const& oblate_earth_;

  // A spherical Earth with no other bodies.
  not_null<std::unique_ptr<Ephemeris<ICRS>>> spherical_earth_ephemeris_;
  MassiveBody const& spherical_earth_;
};


TEST_F(OrbitalElementsTest, KeplerOrbit) {
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
                              J2000 + 10 * Day,
                              *spherical_earth_ephemeris_),
      spherical_earth_,
      MasslessBody{});
  ASSERT_THAT(status_or_elements, IsOk());
  OrbitalElements const& elements = status_or_elements.ValueOrDie();
  EXPECT_THAT(
      AbsoluteError(*initial_osculating.period, elements.anomalistic_period()),
      IsNear(510 * Micro(Second)));
  EXPECT_THAT(
      AbsoluteError(*initial_osculating.period, elements.nodal_period()),
      IsNear(3'100 * Micro(Second)));
  EXPECT_THAT(
      AbsoluteError(*initial_osculating.period, elements.sidereal_period()),
      IsNear(0.92 * Micro(Second)));

  EXPECT_THAT(elements.nodal_precession(), IsNear(1.0 * Degree / JulianYear));

  // Mean element values.
  EXPECT_THAT(AbsoluteError(*initial_osculating.semimajor_axis,
                            elements.mean_semimajor_axis_range().midpoint()),
              IsNear(100 * Micro(Metre)));
  EXPECT_THAT(AbsoluteError(*initial_osculating.eccentricity,
                            elements.mean_eccentricity_range().midpoint()),
              IsNear(1.5e-11));
  EXPECT_THAT(AbsoluteError(initial_osculating.inclination,
                            elements.mean_inclination_range().midpoint()),
              IsNear(0.56 * Micro(ArcSecond)));
  EXPECT_THAT(AbsoluteError(
                  initial_osculating.longitude_of_ascending_node,
                  elements.mean_longitude_of_ascending_node_range().midpoint()),
              IsNear(54 * ArcSecond));
  EXPECT_THAT(
      AbsoluteError(*initial_osculating.argument_of_periapsis,
                    elements.mean_argument_of_periapsis_range().midpoint()),
      IsNear(61 * ArcSecond));

  // Mean element stability.
  EXPECT_THAT(elements.mean_semimajor_axis_range().measure(),
              IsNear(1.0 * Milli(Metre)));
  EXPECT_THAT(elements.mean_eccentricity_range().measure(), IsNear(1.0e-10));
  EXPECT_THAT(elements.mean_inclination_range().measure(),
              IsNear(0.61 * Micro(ArcSecond)));
  EXPECT_THAT(elements.mean_longitude_of_ascending_node_range().measure(),
              IsNear(1.9 * ArcMinute));
  EXPECT_THAT(elements.mean_argument_of_periapsis_range().measure(),
              IsNear(2.2 * ArcMinute));

  OFStream f(SOLUTION_DIR / "mathematica" /
             "unperturbed_elements.generated.wl");
  f << mathematica::Assign("unperturbedOsculating",
                           elements.osculating_equinoctial_elements());
  f << mathematica::Assign("unperturbedMean",
                           elements.mean_equinoctial_elements());
}

TEST_F(OrbitalElementsTest, J2Perturbation) {
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
                              *oblate_earth_ephemeris_),
      oblate_earth_,
      MasslessBody{});
  ASSERT_THAT(status_or_elements, IsOk());
  OrbitalElements const& elements = status_or_elements.ValueOrDie();
  EXPECT_THAT(elements.anomalistic_period() - *initial_osculating.period,
              IsNear(-7.8 * Second));
  EXPECT_THAT(elements.nodal_period() - *initial_osculating.period,
              IsNear(-23 * Second));
  EXPECT_THAT(elements.sidereal_period() - *initial_osculating.period,
              IsNear(-16 * Second));

  // The notation for the computation of the theoretical precessions follows
  // Capderou (2012), Satellites : de Kepler au GPS, section 7.1.1.
  double const η = elements.mean_semimajor_axis_range().midpoint() /
                   oblate_earth_.reference_radius();
  Angle const i = elements.mean_inclination_range().midpoint();
  AngularFrequency const K0 = 3.0 / 2.0 * oblate_earth_.j2() *
                              Sqrt(oblate_earth_.gravitational_parameter() /
                                   Pow<3>(oblate_earth_.reference_radius())) *
                              Radian;
  // See (7.6).
  AngularFrequency const theoretical_Ωʹ = -K0 * Pow<-7>(Sqrt(η)) * Cos(i);
  // See (7.13).
  AngularFrequency const theoretical_ωʹ =
      0.5 * K0 * Pow<-7>(Sqrt(η)) * (5 * Pow<2>(Cos(i)) - 1);

  EXPECT_THAT(theoretical_Ωʹ, IsNear(-7.2 * Degree / Day));
  EXPECT_THAT(theoretical_ωʹ, IsNear(14 * Degree / Day));

  EXPECT_THAT(elements.nodal_precession(), IsNear(theoretical_Ωʹ));

  // Mean element values.  Since Ω and ω precess rapidly, the midpoint of the
  // range of values is of no interest.
  EXPECT_THAT(AbsoluteError(*initial_osculating.semimajor_axis,
                            elements.mean_semimajor_axis_range().midpoint()),
              IsNear(25 * Metre));
  EXPECT_THAT(elements.mean_eccentricity_range().midpoint(), IsNear(0.0013));
  EXPECT_THAT(AbsoluteError(initial_osculating.inclination,
                            elements.mean_inclination_range().midpoint()),
              IsNear(1.9 * Micro(ArcSecond)));

  // Mean element stability: Ω and ω precess as expected, the other elements are
  // stable.
  EXPECT_THAT(elements.mean_semimajor_axis_range().measure(),
              IsNear(70 * Milli(Metre)));
  EXPECT_THAT(elements.mean_eccentricity_range().measure(), IsNear(3.7e-9));
  EXPECT_THAT(elements.mean_inclination_range().measure(),
              IsNear(1.1 * Micro(ArcSecond)));
  EXPECT_THAT(elements.mean_longitude_of_ascending_node_range().measure(),
              IsNear(-theoretical_Ωʹ * mission_duration));
  EXPECT_THAT(elements.mean_argument_of_periapsis_range().measure(),
              IsNear(theoretical_ωʹ * mission_duration));

  OFStream f(SOLUTION_DIR / "mathematica" /
             "j2_perturbed_elements.generated.wl");
  f << mathematica::Assign("j2PerturbedOsculating",
                           elements.osculating_equinoctial_elements());
  f << mathematica::Assign("j2PerturbedMean",
                           elements.mean_equinoctial_elements());
}

TEST_F(OrbitalElementsTest, RealPerturbation) {
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
                              *real_ephemeris_),
      real_earth_,
      MasslessBody{});
  ASSERT_THAT(status_or_elements, IsOk());
  OrbitalElements const& elements = status_or_elements.ValueOrDie();
  EXPECT_THAT(elements.anomalistic_period() - *initial_osculating.period,
              IsNear(-7.8 * Second));
  EXPECT_THAT(elements.nodal_period() - *initial_osculating.period,
              IsNear(-14 * Second));
  EXPECT_THAT(elements.sidereal_period() - *initial_osculating.period,
              IsNear(-16 * Second));

  // This value is meaningless, see below.
  EXPECT_THAT(elements.nodal_precession(), IsNear(2.0 * Degree / Day));

  // Mean element values.
  EXPECT_THAT(AbsoluteError(*initial_osculating.semimajor_axis,
                            elements.mean_semimajor_axis_range().midpoint()),
              IsNear(100 * Metre));
  EXPECT_THAT(elements.mean_eccentricity_range().midpoint(), IsNear(0.0014));
  EXPECT_THAT(AbsoluteError(initial_osculating.inclination,
                            elements.mean_inclination_range().midpoint()),
              IsNear(6.0 * ArcSecond));

  // Mean element stability: Ω and ω exhibit a daily oscillation (likely due to
  // the tesseral terms of the geopotential) as the very low inclination means
  // that these elements are singular.
  // The other elements are stable.
  // A closer analysis would show that the longitude of periapsis ϖ = Ω + ω
  // exhibits a precession that is largely free of oscillations, and that, if
  // its oscillations are filtered, the argument of periapsis ω precesses as
  // expected; the longitude of the ascending node Ω exhibits no obvious
  // precession even if its daily oscillation is filtered out.
  EXPECT_THAT(elements.mean_semimajor_axis_range().measure(),
              IsNear(20 * Metre));
  EXPECT_THAT(elements.mean_eccentricity_range().measure(), IsNear(9.2e-5));
  EXPECT_THAT(elements.mean_inclination_range().measure(),
              IsNear(11 * ArcSecond));
  EXPECT_THAT(elements.mean_longitude_of_ascending_node_range().measure(),
              IsNear(140 * Degree));
  EXPECT_THAT(elements.mean_argument_of_periapsis_range().measure(),
              IsNear(150 * Degree));

  OFStream f(SOLUTION_DIR / "mathematica" /
             "fully_perturbed_elements.generated.wl");
  f << mathematica::Assign("fullyPerturbedOsculating",
                           elements.osculating_equinoctial_elements());
  f << mathematica::Assign("fullyPerturbedMean",
                           elements.mean_equinoctial_elements());
}

}  // namespace astronomy
}  // namespace principia
