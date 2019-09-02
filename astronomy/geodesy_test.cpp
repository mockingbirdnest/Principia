
#include <limits>

#include "astronomy/frames.hpp"
#include "astronomy/standard_product_3.hpp"
#include "base/bundle.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/body_surface_dynamic_frame.hpp"
#include "physics/solar_system.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/numerics.hpp"

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
using testing_utilities::operator""_⑴;
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

  SolarSystem<ICRS> const solar_system_2010_;
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

TEST_F(GeodesyTest, DISABLED_LAGEOS2) {
  MasslessBody lageos2;

  StandardProduct3 initial_ilrsa(SOLUTION_DIR / "astronomy" /
                                     "standard_product_3" /
                                     "ilrsa.orb.lageos2.160319.v35.sp3",
                                 StandardProduct3::Dialect::ILRSA);
  StandardProduct3 initial_ilrsb(SOLUTION_DIR / "astronomy" /
                                     "standard_product_3" /
                                     "ilrsb.orb.lageos2.160319.v35.sp3",
                                 StandardProduct3::Dialect::ILRSB);
  StandardProduct3 final_ilrsa(SOLUTION_DIR / "astronomy" /
                                   "standard_product_3" /
                                   "ilrsa.orb.lageos2.180804.v70.sp3",
                               StandardProduct3::Dialect::ILRSA);

  StandardProduct3::SatelliteIdentifier const lageos2_id{
      StandardProduct3::SatelliteGroup::General, 52};

  CHECK_EQ(initial_ilrsa.orbit(lageos2_id).front()->front().time,
           initial_ilrsb.orbit(lageos2_id).front()->front().time);

  Instant const initial_time =
      initial_ilrsa.orbit(lageos2_id).front()->front().time;
  DegreesOfFreedom<ITRS> const initial_dof_ilrsa =
      initial_ilrsa.orbit(lageos2_id).front()->front().degrees_of_freedom;

  DegreesOfFreedom<ITRS> const initial_dof_ilrsb =
      initial_ilrsb.orbit(lageos2_id).front()->front().degrees_of_freedom;

  Instant const final_time =
      final_ilrsa.orbit(lageos2_id).front()->front().time;
  DegreesOfFreedom<ITRS> const expected_final_dof =
      final_ilrsa.orbit(lageos2_id).front()->front().degrees_of_freedom;

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
            /*max_ephemeris_steps=*/std::numeric_limits<std::int64_t>::max());
  };
  Bundle bundle;
  bundle.Add([&flow_lageos2, &primary_lageos2_trajectory]() {
    return flow_lageos2(primary_lageos2_trajectory);
  });
  bundle.Add([&flow_lageos2, &secondary_lageos2_trajectory]() {
    return flow_lageos2(secondary_lageos2_trajectory);
  });
  bundle.Join();
  EXPECT_THAT(primary_lageos2_trajectory.back().time, Eq(final_time));
  EXPECT_THAT(secondary_lageos2_trajectory.back().time, Eq(final_time));

  auto const primary_actual_final_dof = itrs_.ToThisFrameAtTime(final_time)(
      primary_lageos2_trajectory.back().degrees_of_freedom);
  auto const secondary_actual_final_dof = itrs_.ToThisFrameAtTime(final_time)(
      secondary_lageos2_trajectory.back().degrees_of_freedom);

  // Absolute error in position.
  EXPECT_THAT(AbsoluteError(primary_actual_final_dof.position(),
                            expected_final_dof.position()),
              IsNear(191_⑴ * Kilo(Metre)));
  // Angular error at the geocentre.
  EXPECT_THAT(AngleBetween(primary_actual_final_dof.position() - ITRS::origin,
                           expected_final_dof.position() - ITRS::origin),
              IsNear(53_⑴ * ArcMinute));
  // Radial error at the geocentre.
  EXPECT_THAT(
      AbsoluteError((primary_actual_final_dof.position() - ITRS::origin).Norm(),
                    (expected_final_dof.position() - ITRS::origin).Norm()),
      IsNear(0.89_⑴ * Kilo(Metre)));

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
              IsNear(42_⑴ * Metre));
  EXPECT_THAT(AbsoluteError(*actual_elements.apoapsis_distance,
                            *expected_elements.apoapsis_distance),
              IsNear(9.5_⑴ * Metre));
  EXPECT_THAT(AbsoluteError(actual_elements.longitude_of_ascending_node,
                            expected_elements.longitude_of_ascending_node),
              IsNear(17_⑴ * ArcSecond));
  EXPECT_THAT(AbsoluteError(actual_elements.inclination,
                            expected_elements.inclination),
              IsNear(1.1_⑴ * ArcSecond));
  EXPECT_THAT(AbsoluteError(*actual_elements.argument_of_periapsis,
                            *expected_elements.argument_of_periapsis),
              IsNear(217_⑴ * ArcSecond));
  EXPECT_THAT(AbsoluteError(*actual_elements.mean_anomaly,
                            *expected_elements.mean_anomaly),
              IsNear(58_⑴ * ArcMinute));

#if 0
  // Error arising from uncertainty in the initial state, estimated as the
  // difference between the primary and secondary ILRS products.
  // Absolute error in position.
  EXPECT_THAT(AbsoluteError(secondary_actual_final_dof.position(),
                            primary_actual_final_dof.position()),
              AnyOf(IsNear(237_⑴ * Metre),    // Linux.
                    IsNear(28_⑴ * Metre),     // No FMA.
                    IsNear(98_⑴ * Metre),     // FMA.
                    IsNear(220_⑴ * Metre)));  // VS 2019.
  // Angular error at the geocentre.
  EXPECT_THAT(AngleBetween(secondary_actual_final_dof.position() - ITRS::origin,
                           primary_actual_final_dof.position() - ITRS::origin),
                AnyOf(IsNear(4.0_⑴ * ArcSecond),    // Linux.
                      IsNear(0.47_⑴ * ArcSecond),   // No FMA.
                      IsNear(1.6_⑴ * ArcSecond),    // FMA.
                      IsNear(3.7_⑴ * ArcSecond)));  // VS 2019.
  // Radial error at the geocentre.
  EXPECT_THAT(AbsoluteError(
                  (secondary_actual_final_dof.position() - ITRS::origin).Norm(),
                  (primary_actual_final_dof.position() - ITRS::origin).Norm()),
              AnyOf(IsNear(99_⑴ * Centi(Metre)),    // Linux.
                    IsNear(11_⑴ * Centi(Metre)),    // No FMA.
                    IsNear(43_⑴ * Centi(Metre)),    // FMA.
                    IsNear(93_⑴ * Centi(Metre))));  // VS 2019.
#endif
}

#endif

}  // namespace astronomy
}  // namespace principia
