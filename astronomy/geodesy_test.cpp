
#include <limits>

#include "astronomy/frames.hpp"
#include "base/bundle.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/body_surface_dynamic_frame.hpp"
#include "physics/solar_system.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/is_near.hpp"

namespace principia {
namespace astronomy {

using base::Bundle;
using base::dynamic_cast_not_null;
using base::not_null;
using base::Status;
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
using quantities::si::Centi;
using quantities::si::Deci;
using quantities::si::Degree;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Minute;
using quantities::si::Second;
using testing_utilities::AbsoluteError;
using testing_utilities::IsNear;
using ::testing::AnyOf;
using ::testing::Eq;

class GeodesyTest : public ::testing::Test {
 protected:
  GeodesyTest()
      : solar_system_2010_(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2455200_500000000.proto.txt"),
        ephemeris_(solar_system_2010_.MakeEphemeris(
            /*accuracy_parameters=*/{/*fitting_tolerance=*/5 * Milli(Metre),
                                     /*geopotential_tolerance=*/0x1p-24},
            Ephemeris<ICRS>::FixedStepParameters(
                SymmetricLinearMultistepIntegrator<QuinlanTremaine1990Order12,
                                                   Position<ICRS>>(),
                /*step=*/10 * Minute))),
        earth_(dynamic_cast_not_null<OblateBody<ICRS> const*>(
            solar_system_2010_.massive_body(*ephemeris_, "Earth"))),
        earth_trajectory_(*ephemeris_->trajectory(earth_)),
        itrs_(ephemeris_.get(), earth_) {}

  SolarSystem<ICRS> solar_system_2010_;
  not_null<std::unique_ptr<Ephemeris<ICRS>>> const ephemeris_;
  not_null<OblateBody<ICRS> const*> const earth_;
  ContinuousTrajectory<ICRS> const& earth_trajectory_;
  // NOTE(egg): using the WGCCRE 2009 elements instead of the proper ERA as
  // defined by the IERS induces the following errors between 2000 and 2020:
  // - a systematic error of 2.8 arcminutes;
  // - a systematic drift of 2.6 milliarcseconds per year;
  // - short-period errors at the milliarcsecond level;
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
  DegreesOfFreedom<ITRS> const initial_dof_ilrsa = {
      ITRS::origin + Displacement<ITRS>({  2505.232029 * Kilo(Metre),
                                         -10564.815741 * Kilo(Metre),
                                          -5129.314404 * Kilo(Metre)}),
      Velocity<ITRS>({ 34323.584344 * Deci(Metre) / Second,
                      -10455.947225 * Deci(Metre) / Second,
                       38998.988146 * Deci(Metre) / Second})};

  // ilrsb.orb.lageos2.160319.v35.sp3, headers and first record, from
  // ftp://cddis.gsfc.nasa.gov/pub/slr/products/orbits/lageos2/160319/.
  // #cV2016  3 13  0  0  0.00000000    5041   SLR ITRF97 FIT JCET
  // ## 1888      0.00000000   120.00000000 57460 0.0000000000000
  // +   1    L52  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  // [... Lines 4 through 12 omitted                         ...]
  // %c L  cc UTC ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
  // %c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
  // %f  0.0000000  0.000000000  0.00000000000  0.000000000000000
  // [... Lines 15 through 18 omitted                        ...]
  // %/* ilrsb.orb.lageos2.160319.v35.sp3 Reference TRF: SLRF2008
  // %/* Input orbits:  ASI v35 GRGS v35 NSGF v35 ESA v35
  // %/* GFZ v35 DGFI v35 JCET v35
  // %/* Combination details in README_CC.ilrsb
  // * 2016  3 13  0  0  0.00000000
  // PL52   2505.232038 -10564.815750  -5129.314387
  // VL52  34323.584276 -10455.947218  38998.988200

  DegreesOfFreedom<ITRS> const initial_dof_ilrsb = {
      ITRS::origin + Displacement<ITRS>({  2505.232038 * Kilo(Metre),
                                         -10564.815750 * Kilo(Metre),
                                          -5129.314387 * Kilo(Metre)}),
      Velocity<ITRS>({ 34323.584276 * Deci(Metre) / Second,
                      -10455.947218 * Deci(Metre) / Second,
                       38998.988200 * Deci(Metre) / Second})};

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

  DiscreteTrajectory<ICRS> primary_lageos2_trajectory;
  primary_lageos2_trajectory.Append(
      initial_time, itrs_.FromThisFrameAtTime(initial_time)(initial_dof_ilrsa));
  DiscreteTrajectory<ICRS> secondary_lageos2_trajectory;
  secondary_lageos2_trajectory.Append(
      initial_time, itrs_.FromThisFrameAtTime(initial_time)(initial_dof_ilrsb));
  auto flow_lageos2 =
      [this,
       final_time](DiscreteTrajectory<ICRS>& lageos2_trajectory) -> Status {
        return ephemeris_->FlowWithAdaptiveStep(
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
            /*last_point_only=*/true);
  };
  Bundle bundle;
  bundle.Add([&flow_lageos2, &primary_lageos2_trajectory]() {
    return flow_lageos2(primary_lageos2_trajectory);
  });
  bundle.Add([&flow_lageos2, &secondary_lageos2_trajectory]() {
    return flow_lageos2(secondary_lageos2_trajectory);
  });
  bundle.Join();
  EXPECT_THAT(primary_lageos2_trajectory.last().time(), Eq(final_time));
  EXPECT_THAT(secondary_lageos2_trajectory.last().time(), Eq(final_time));

  auto const primary_actual_final_dof = itrs_.ToThisFrameAtTime(final_time)(
      primary_lageos2_trajectory.last().degrees_of_freedom());
  auto const secondary_actual_final_dof = itrs_.ToThisFrameAtTime(final_time)(
      secondary_lageos2_trajectory.last().degrees_of_freedom());

  // Absolute error in position.
  EXPECT_THAT(AbsoluteError(primary_actual_final_dof.position(),
                            expected_final_dof.position()),
              IsNear(191 * Kilo(Metre)));
  // Angular error at the geocentre.
  EXPECT_THAT(AngleBetween(primary_actual_final_dof.position() - ITRS::origin,
                           expected_final_dof.position() - ITRS::origin),
              IsNear(53 * ArcMinute));
  // Radial error at the geocentre.
  EXPECT_THAT(
      AbsoluteError((primary_actual_final_dof.position() - ITRS::origin).Norm(),
                    (expected_final_dof.position() - ITRS::origin).Norm()),
      IsNear(894 * Metre));

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
          itrs_.FromThisFrameAtTime(final_time)(primary_actual_final_dof) -
              earth_trajectory_.EvaluateDegreesOfFreedom(final_time),
          final_time).elements_at_epoch();

  EXPECT_THAT(AbsoluteError(*actual_elements.periapsis_distance,
                            *expected_elements.periapsis_distance),
              IsNear(42 * Metre));
  EXPECT_THAT(AbsoluteError(*actual_elements.apoapsis_distance,
                            *expected_elements.apoapsis_distance),
              IsNear(9.5 * Metre));
  EXPECT_THAT(AbsoluteError(actual_elements.longitude_of_ascending_node,
                            expected_elements.longitude_of_ascending_node),
              IsNear(17 * ArcSecond));
  EXPECT_THAT(AbsoluteError(actual_elements.inclination,
                            expected_elements.inclination),
              IsNear(1.1 * ArcSecond));
  EXPECT_THAT(AbsoluteError(*actual_elements.argument_of_periapsis,
                            *expected_elements.argument_of_periapsis),
              IsNear(3 * ArcMinute + 37 * ArcSecond));
  EXPECT_THAT(AbsoluteError(*actual_elements.mean_anomaly,
                            *expected_elements.mean_anomaly),
              IsNear(58 * ArcMinute));

  // Error arising from uncertainty in the initial state, estimated as the
  // difference between the primary and secondary ILRS products.
  // Absolute error in position.
  EXPECT_THAT(AbsoluteError(secondary_actual_final_dof.position(),
                            primary_actual_final_dof.position()),
              AnyOf(IsNear(237 * Metre),    // Linux.
                    IsNear(28 * Metre),     // No FMA.
                    IsNear(8.9 * Metre)));  // FMA.
  // Angular error at the geocentre.
  EXPECT_THAT(AngleBetween(secondary_actual_final_dof.position() - ITRS::origin,
                           primary_actual_final_dof.position() - ITRS::origin),
              AnyOf(IsNear(4.0 * ArcSecond),     // Linux.
                    IsNear(0.47 * ArcSecond),    // No FMA.
                    IsNear(0.15 * ArcSecond)));  // FMA.
  // Radial error at the geocentre.
  EXPECT_THAT(AbsoluteError(
                  (secondary_actual_final_dof.position() - ITRS::origin).Norm(),
                  (primary_actual_final_dof.position() - ITRS::origin).Norm()),
              AnyOf(IsNear(99 * Centi(Metre)),     // Linux.
                    IsNear(11 * Centi(Metre)),     // No FMA.
                    IsNear(3.3 * Centi(Metre))));  // FMA.
}

#endif

}  // namespace astronomy
}  // namespace principia
