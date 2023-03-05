#include "physics/body_surface_dynamic_frame.hpp"

#include <memory>

#include "astronomy/frames.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/named_quantities.hpp"
#include "geometry/rotation.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "integrators/methods.hpp"
#include "integrators/symplectic_runge_kutta_nyström_integrator.hpp"
#include "physics/ephemeris.hpp"
#include "physics/solar_system.hpp"
#include "quantities/constants.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "serialization/physics.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace physics {
namespace internal_body_surface_dynamic_frame {

using astronomy::ICRS;
using ::testing::IsNull;
using ::testing::Lt;
using ::testing::Not;
using ::testing::Return;
using ::testing::_;
using namespace principia::base::_not_null;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::geometry::_named_quantities;
using namespace principia::integrators::_methods;
using namespace principia::integrators::
    _symplectic_runge_kutta_nyström_integrator;
using namespace principia::quantities::_elementary_functions;
using namespace principia::quantities::_named_quantities;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::testing_utilities::_almost_equals;
using namespace principia::testing_utilities::_componentwise;
using namespace principia::testing_utilities::_numerics;
using namespace principia::testing_utilities::_vanishes_before;

namespace {

char constexpr big[] = "Big";
char constexpr small[] = "Small";

}  // namespace

class BodySurfaceDynamicFrameTest : public ::testing::Test {
 protected:
  // The rotating frame centred on the big body and directed to the small one.
  using BigSmallFrame = Frame<serialization::Frame::TestTag,
                              Arbitrary,
                              Handedness::Right,
                              serialization::Frame::TEST>;

  BodySurfaceDynamicFrameTest()
      : period_(10 * π * sqrt(5.0 / 7.0) * Second),
        solar_system_(SOLUTION_DIR / "astronomy" /
                          "test_gravity_model_two_bodies.proto.txt",
                      SOLUTION_DIR / "astronomy" /
                          "test_initial_state_two_bodies_circular.proto.txt"),
        t0_(solar_system_.epoch()),
        ephemeris_(solar_system_.MakeEphemeris(
            /*accuracy_parameters=*/{/*fitting_tolerance=*/1 * Milli(Metre),
                                     /*geopotential_tolerance=*/0x1p-24},
            Ephemeris<ICRS>::FixedStepParameters(
                SymplecticRungeKuttaNyströmIntegrator<
                    McLachlanAtela1992Order4Optimal,
                    Ephemeris<ICRS>::NewtonianMotionEquation>(),
                /*step=*/10 * Milli(Second)))),
        big_(dynamic_cast_not_null<RotatingBody<ICRS> const*>(
            solar_system_.massive_body(*ephemeris_, big))),
        big_initial_state_(solar_system_.degrees_of_freedom(big)),
        big_gravitational_parameter_(
            solar_system_.gravitational_parameter(big)),
        small_initial_state_(solar_system_.degrees_of_freedom(small)),
        small_gravitational_parameter_(
            solar_system_.gravitational_parameter(small)),
        // A body that rotates at the same speed as the one in
        // BodyCentredBodyDirectionDynamicFrameTest, so it produces the same
        // fictitious forces.
        centre_(MassiveBody::Parameters(1 * Kilogram),
                RotatingBody<ICRS>::Parameters(
                    /*mean_radius=*/1 * Metre,
                    /*reference_angle=*/0 * Radian,
                    /*reference_instant=*/t0_,
                    /*angular_frequency=*/10 * Radian / Second,
                    /*ascension_of_pole=*/0 * Radian,
                    /*declination_of_pole=*/π / 2 * Radian)),
        massive_centre_(&centre_) {
    EXPECT_OK(ephemeris_->Prolong(t0_ + 2 * period_));
    big_frame_ = std::make_unique<BodySurfaceDynamicFrame<ICRS, BigSmallFrame>>(
        ephemeris_.get(), big_);
  }

  Time const period_;
  SolarSystem<ICRS> solar_system_;
  Instant const t0_;
  std::unique_ptr<Ephemeris<ICRS>> const ephemeris_;
  RotatingBody<ICRS> const* const big_;
  DegreesOfFreedom<ICRS> big_initial_state_;
  GravitationalParameter big_gravitational_parameter_;
  DegreesOfFreedom<ICRS> small_initial_state_;
  GravitationalParameter small_gravitational_parameter_;
  RotatingBody<ICRS> const centre_;
  not_null<MassiveBody const*> const massive_centre_;

  std::unique_ptr<BodySurfaceDynamicFrame<ICRS, BigSmallFrame>> big_frame_;
};


TEST_F(BodySurfaceDynamicFrameTest, ToBigSmallFrameAtTime) {
  int const steps = 100;

  for (Instant t = t0_; t < t0_ + 1 * period_; t += period_ / steps) {
    auto const to_big_frame_at_t = big_frame_->ToThisFrameAtTime(t);

    // Check that the big body is at the origin and doesn't move.  Check that
    // the small body is at a fixed position in the sky.
    DegreesOfFreedom<ICRS> const big_in_inertial_frame_at_t =
        solar_system_.trajectory(*ephemeris_, big).
            EvaluateDegreesOfFreedom(t);
    DegreesOfFreedom<ICRS> const small_in_inertial_frame_at_t =
        solar_system_.trajectory(*ephemeris_, small).
            EvaluateDegreesOfFreedom(t);

    DegreesOfFreedom<BigSmallFrame> const big_in_big_small_at_t =
        to_big_frame_at_t(big_in_inertial_frame_at_t);
    DegreesOfFreedom<BigSmallFrame> const small_in_big_small_at_t =
        to_big_frame_at_t(small_in_inertial_frame_at_t);
    EXPECT_THAT(AbsoluteError(big_in_big_small_at_t.position(),
                              BigSmallFrame::origin),
                Lt(1.0e-6 * Metre));
    EXPECT_THAT(AbsoluteError(big_in_big_small_at_t.velocity(),
                              BigSmallFrame::unmoving),
                Lt(1.0e-4 * Metre / Second));
    EXPECT_THAT(AbsoluteError(small_in_big_small_at_t.position(),
                              Displacement<BigSmallFrame>({
                                  0 * Kilo(Metre),
                                  5.0 * Kilo(Metre),
                                  0 * Kilo(Metre)}) + BigSmallFrame::origin),
                Lt(2.7e-4 * Metre));
    EXPECT_THAT(AbsoluteError(small_in_big_small_at_t.velocity(),
                              BigSmallFrame::unmoving),
                Lt(4.0e-3 * Metre / Second));
  }
}

TEST_F(BodySurfaceDynamicFrameTest, Inverse) {
  int const steps = 100;
  for (Instant t = t0_; t < t0_ + 1 * period_; t += period_ / steps) {
    auto const from_big_frame_at_t =
        big_frame_->FromThisFrameAtTime(t);
    auto const to_big_frame_at_t = big_frame_->ToThisFrameAtTime(t);
    auto const small_initial_state_transformed_and_back =
        from_big_frame_at_t(to_big_frame_at_t(
            small_initial_state_));
    EXPECT_THAT(
        AbsoluteError(small_initial_state_transformed_and_back.position(),
                      small_initial_state_.position()),
        Lt(1.0e-11 * Metre));
    EXPECT_THAT(
        AbsoluteError(small_initial_state_transformed_and_back.velocity(),
                      small_initial_state_.velocity()),
        Lt(1.0e-11 * Metre / Second));
  }
}

TEST_F(BodySurfaceDynamicFrameTest, GeometricAcceleration) {
  Instant const t = t0_ + period_;
  DegreesOfFreedom<BigSmallFrame> const point_dof =
      {Displacement<BigSmallFrame>({10 * Metre, 20 * Metre, 30 * Metre}) +
           BigSmallFrame::origin,
       Velocity<BigSmallFrame>({3 * Metre / Second,
                                2 * Metre / Second,
                                1 * Metre / Second})};
  // We trust the functions to compute the values correctly, but this test
  // ensures that we don't get NaNs.
  EXPECT_THAT(big_frame_->GeometricAcceleration(t, point_dof),
              AlmostEquals(Vector<Acceleration, BigSmallFrame>({
                  -9.54504983899937710e5 * Metre / Pow<2>(Second),
                  -1.90900569196062256e6 * Metre / Pow<2>(Second),
                  -2.86351379198155506e6 * Metre / Pow<2>(Second)}), 0, 2));
}

TEST_F(BodySurfaceDynamicFrameTest, Serialization) {
  serialization::DynamicFrame message;
  big_frame_->WriteToMessage(&message);

  EXPECT_TRUE(message.HasExtension(
      serialization::BodySurfaceDynamicFrame::extension));
  auto const extension = message.GetExtension(
      serialization::BodySurfaceDynamicFrame::extension);
  EXPECT_TRUE(extension.has_centre());
  EXPECT_EQ(0, extension.centre());

  auto const read_big_frame =
      DynamicFrame<ICRS, BigSmallFrame>::ReadFromMessage(message,
                                                         ephemeris_.get());
  EXPECT_THAT(read_big_frame, Not(IsNull()));

  Instant const t = t0_ + period_;
  DegreesOfFreedom<BigSmallFrame> const point_dof =
      {Displacement<BigSmallFrame>({10 * Metre, 20 * Metre, 30 * Metre}) +
           BigSmallFrame::origin,
       Velocity<BigSmallFrame>({3 * Metre / Second,
                                2 * Metre / Second,
                                1 * Metre / Second})};
  EXPECT_EQ(big_frame_->GeometricAcceleration(t, point_dof),
            read_big_frame->GeometricAcceleration(t, point_dof));
}

}  // namespace internal_body_surface_dynamic_frame
}  // namespace physics
}  // namespace principia
