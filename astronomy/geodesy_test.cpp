
#include "astronomy/frames.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/body_surface_dynamic_frame.hpp"
#include "physics/solar_system.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/is_near.hpp"

namespace principia {
namespace astronomy {

using base::dynamic_cast_not_null;
using base::not_null;
using geometry::AngleBetween;
using geometry::Displacement;
using geometry::Instant;
using geometry::Position;
using geometry::Velocity;
using integrators::EmbeddedExplicitRungeKuttaNyströmIntegrator;
using integrators::SymmetricLinearMultistepIntegrator;
using integrators::methods::DormandالمكاوىPrince1986RKN434FM;
using integrators::methods::Quinlan1999Order8A;
using integrators::methods::QuinlanTremaine1990Order12;
using physics::BodySurfaceDynamicFrame;
using physics::ContinuousTrajectory;
using physics::DegreesOfFreedom;
using physics::DiscreteTrajectory;
using physics::Ephemeris;
using physics::KeplerianElements;
using physics::KeplerOrbit;
using physics::MasslessBody;
using physics::OblateBody;
using physics::SolarSystem;
using quantities::si::ArcMinute;
using quantities::si::ArcSecond;
using quantities::si::Deci;
using quantities::si::Degree;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Second;
using testing_utilities::AbsoluteError;
using testing_utilities::IsNear;
using ::testing::Eq;

class GeodesyTest : public ::testing::Test {
 protected:
  GeodesyTest()
      : solar_system_2000_(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2451545_000000000.proto.txt"),
        ephemeris_(solar_system_2000_.MakeEphemeris(
            /*fitting_tolerance=*/5 * Milli(Metre),
            Ephemeris<ICRS>::FixedStepParameters(
                SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                                   Position<ICRS>>(),
                /*step=*/10 * Minute))),
        earth_(dynamic_cast_not_null<OblateBody<ICRS> const*>(
            solar_system_2000_.massive_body(*ephemeris_, "Earth"))),
        earth_trajectory_(*ephemeris_->trajectory(earth_)),
        itrs_(ephemeris_.get(), earth_) {}

  SolarSystem<ICRS> solar_system_2000_;
  not_null<std::unique_ptr<Ephemeris<ICRS>>> const ephemeris_;
  not_null<OblateBody<ICRS> const*> const earth_;
  ContinuousTrajectory<ICRS> const& earth_trajectory_;
  // TODO(egg): find a bound on the error induced by this approximation to the
  // Earth Rotation Angle.
  BodySurfaceDynamicFrame<ICRS, ITRS> itrs_;
};

#if !defined(_DEBUG)

TEST_F(GeodesyTest, LAGEOS2) {
  MasslessBody lageos2;
  // Initial state and expected state from the ILRS official primary analysis
  // product ILRSA, see
  // https://ilrs.cddis.eosdis.nasa.gov/data_and_products/products/index.html;
  // see also the definition of the SP3 format
  // ftp://igs.org/pub/data/format/sp3c.txt.

  // ilrsa.orb.lageos2.160319.v35.sp3, headers and first record, from
  // ftp://cddis.gsfc.nasa.gov/pub/slr/products/orbits/lageos2/160319/.
  // #cV2016  3 13  0  0  0.00000000    5040   SLR SLR08 FIT COMB
  // ## 1888      0.00000000   120.00000000 57460 0.0000000000000
  // +    1   L52  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  // [... Lines 4 through 12 omitted                         ...]
  // %c L  cc UTC ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
  // %c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
  // %f  0.0000000  0.000000000  0.00000000000  0.000000000000000
  // [... Lines 15 through 18 omitted                        ...]
  // %/* ilrsa.orb.lageos2.160319.v35.sp3 Reference TRF: SLRF2008
  // %/* Input orbits: ASI v35, BKG v35, DGFI v35, ESA v35,
  // %/* GFZ v35, GRGS v35, JCET v35, NSGF v35,
  // %/* Combination details in README_CC.ilrsa
  // *  2016  3 13  0  0  0.00000000
  // PL52   2505.232029 -10564.815741  -5129.314404 999999.999999
  // VL52  34323.584344 -10455.947225  38998.988146 999999.999999

  constexpr Instant initial_time = "2016-03-13T00:00:00,000"_UTC;
  DegreesOfFreedom<ITRS> const initial_dof = {
      ITRS::origin + Displacement<ITRS>({2505.232029 * Kilo(Metre),
                                         -10564.815741 * Kilo(Metre),
                                         -5129.314404 * Kilo(Metre)}),
      Velocity<ITRS>({34323.584344 * Deci(Metre) / Second,
                      -10455.947225 * Deci(Metre) / Second,
                      38998.988146 * Deci(Metre) / Second})};

  // ilrsa.orb.lageos2.180804.v70.sp3, headers and first record, from
  // ftp://cddis.gsfc.nasa.gov/pub/slr/products/orbits/lageos2/180804/.
  // #cV2018  7 29  0  0  0.00000000    5040   SLR SLR08 FIT COMB
  // ## 2012      0.00000000   120.00000000 58328 0.0000000000000
  // +    1   L52  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  // [... Lines 4 through 12 omitted                         ...]
  // %c L  cc UTC ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
  // %c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
  // %f  0.0000000  0.000000000  0.00000000000  0.000000000000000
  // [... Lines 15 through 18 omitted                        ...]
  // %/* ilrsa.orb.lageos2.180804.v70.sp3 Reference TRF: SLRF2008
  // %/* Input orbits: ASI v70, BKG v70, DGFI v70, ESA v70,
  // %/* GFZ v70, JCET v70, NSGF v70,  
  // %/* Combination details in README_CC.ilrsa
  // *  2018  7 29  0  0  0.00000000
  // PL52 -11150.750217   5070.184012   1340.324930 999999.999999
  // VL52 -15231.027828 -21132.111357 -44478.560714 999999.999999

  constexpr Instant final_time = "2018-07-29T00:00:00,000"_UTC;
  DegreesOfFreedom<ITRS> const expected_final_dof = {
      ITRS::origin + Displacement<ITRS>({-11150.750217 * Kilo(Metre),
                                         5070.184012 * Kilo(Metre),
                                         1340.324930 * Kilo(Metre)}),
      Velocity<ITRS>({-15231.027828 * Deci(Metre) / Second,
                      -21132.111357 * Deci(Metre) / Second,
                      -44478.560714 * Deci(Metre) / Second})};

  ephemeris_->Prolong(final_time);

  DiscreteTrajectory<ICRS> lageos2_trajectory;
  lageos2_trajectory.Append(
      initial_time, itrs_.FromThisFrameAtTime(initial_time)(initial_dof));
  CHECK_OK(ephemeris_->FlowWithAdaptiveStep(
      &lageos2_trajectory,
      Ephemeris<ICRS>::NoIntrinsicAcceleration,
      final_time,
      Ephemeris<ICRS>::AdaptiveStepParameters(
          EmbeddedExplicitRungeKuttaNyströmIntegrator<
              DormandالمكاوىPrince1986RKN434FM,
              Position<ICRS>>(),
          std::numeric_limits<std::int64_t>::max(),
          /*length_integration_tolerance=*/1 * Milli(Metre),
          /*speed_integration_tolerance=*/1 * Milli(Metre) / Second),
      /*max_ephemeris_steps=*/std::numeric_limits<std::int64_t>::max(),
      /*last_point_only=*/true));
  EXPECT_THAT(lageos2_trajectory.last().time(), Eq(final_time));
  auto const actual_final_dof = itrs_.ToThisFrameAtTime(final_time)(
      lageos2_trajectory.last().degrees_of_freedom());

  // Absolute error in position.
  EXPECT_THAT(
      AbsoluteError(actual_final_dof.position(), expected_final_dof.position()),
      IsNear(570 * Kilo(Metre)));
  // Angular error at the geocentre.
  EXPECT_THAT(AngleBetween(actual_final_dof.position() - ITRS::origin,
                           expected_final_dof.position() - ITRS::origin),
              IsNear(2 * Degree + 41 * ArcMinute));
  // Radial error at the geocentre.
  EXPECT_THAT(
      AbsoluteError((actual_final_dof.position() - ITRS::origin).Norm(),
                    (expected_final_dof.position() - ITRS::origin).Norm()),
      IsNear(270 * Metre));

  // Errors in orbital elements.
  KeplerianElements<ICRS> const expected_elements =
      KeplerOrbit<ICRS>(
          *earth_,
          lageos2,
          itrs_.FromThisFrameAtTime(final_time)(expected_final_dof) -
              earth_trajectory_.EvaluateDegreesOfFreedom(final_time),
          final_time).elements_at_epoch();
  KeplerianElements<ICRS> const actual_elements =
      KeplerOrbit<ICRS>(
          *earth_,
          lageos2,
          itrs_.FromThisFrameAtTime(final_time)(actual_final_dof) -
              earth_trajectory_.EvaluateDegreesOfFreedom(final_time),
          final_time).elements_at_epoch();

  EXPECT_THAT(AbsoluteError(*actual_elements.periapsis_distance,
                            *expected_elements.periapsis_distance),
              IsNear(2 * Kilo(Metre)));
  EXPECT_THAT(AbsoluteError(*actual_elements.apoapsis_distance,
                            *expected_elements.apoapsis_distance),
              IsNear(2 * Kilo(Metre)));
  EXPECT_THAT(AbsoluteError(actual_elements.longitude_of_ascending_node,
                            expected_elements.longitude_of_ascending_node),
              IsNear(4 * ArcMinute + 21 * ArcSecond));
  EXPECT_THAT(AbsoluteError(actual_elements.inclination,
                            expected_elements.inclination),
              IsNear(1.2 * ArcSecond));
  EXPECT_THAT(AbsoluteError(*actual_elements.argument_of_periapsis,
                            *expected_elements.argument_of_periapsis),
              IsNear(15 * ArcMinute + 57 * ArcSecond));
  EXPECT_THAT(AbsoluteError(*actual_elements.mean_anomaly,
                            *expected_elements.mean_anomaly),
              IsNear(2 * Degree + 31 * ArcMinute));
}

#endif

}  // namespace astronomy
}  // namespace principia
