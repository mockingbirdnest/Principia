
#include "physics/body_surface_frame_field.hpp"

#include "astronomy/epoch.hpp"
#include "astronomy/frames.hpp"
#include "geometry/frame.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/mock_continuous_trajectory.hpp"
#include "physics/mock_ephemeris.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {
namespace physics {
namespace internal_body_surface_dynamic_frame {

using astronomy::ICRFJ2000Equator;
using astronomy::J2000;
using geometry::Displacement;
using geometry::InnerProduct;
using geometry::Frame;
using geometry::Rotation;
using geometry::Vector;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::VanishesBefore;
using ::testing::Return;
using ::testing::_;

class BodySurfaceFrameFieldTest : public ::testing::Test {
 protected:
  using TestFrame = Frame<serialization::Frame::TestTag,
                              serialization::Frame::TEST, /*inertial=*/false>;

  BodySurfaceFrameFieldTest()
      : body_(MassiveBody::Parameters(1 * Kilogram),
              RotatingBody<ICRFJ2000Equator>::Parameters(
                  /*mean_radius=*/1 * Metre,
                  /*reference_angle=*/0 * Radian,
                  /*reference_instant=*/J2000,
                  /*angular_frequency=*/10 * Radian / Second,
                  /*ascension_of_pole=*/π / 4 * Radian,
                  /*declination_of_pole=*/π / 4 * Radian)) {
    // The polar axis is {1/2, 1/2, 1/Sqrt[2]}.
    EXPECT_CALL(ephemeris_, trajectory(_)).WillOnce(Return(&trajectory_));
    EXPECT_CALL(trajectory_, EvaluatePosition(J2000, nullptr))
        .WillOnce(Return(ICRFJ2000Equator::origin));
    field_ =
        std::make_unique<BodySurfaceFrameField<ICRFJ2000Equator, TestFrame>>(
            ephemeris_, J2000, &body_);
  }

  MockEphemeris<ICRFJ2000Equator> const ephemeris_;
  RotatingBody<ICRFJ2000Equator> const body_;
  MockContinuousTrajectory<ICRFJ2000Equator> trajectory_;
  std::unique_ptr<BodySurfaceFrameField<ICRFJ2000Equator, TestFrame>> field_;
};

TEST_F(BodySurfaceFrameFieldTest, FromThisFrame) {
  Displacement<ICRFJ2000Equator> const displacement(
      {2 * Metre, 3 * Metre, 4 * Metre});
  Rotation<TestFrame, ICRFJ2000Equator> const rotation =
      field_->FromThisFrame(ICRFJ2000Equator::origin + displacement);
  {
    auto const actual = rotation(Vector<double, TestFrame>({1, 0, 0}));
    auto const expected = Vector<double, ICRFJ2000Equator>(
        {Sqrt((4531 + 1624 * Sqrt(2)) / 8149),
         -2 * Sqrt((419 - 116 * Sqrt(2)) / 8149),
         -Sqrt(((2 * (971 - 580 * Sqrt(2))) / 8149))});
    EXPECT_THAT(actual, AlmostEquals(expected, 43));
    EXPECT_THAT(InnerProduct(actual, displacement),
                VanishesBefore(displacement.Norm(), 11));
  }
  {
    auto const actual = rotation(Vector<double, TestFrame>({0, 1, 0}));
    auto const expected = Vector<double, ICRFJ2000Equator>(
        {-Sqrt((86 - 56 * Sqrt(2)) / 281),
         -2 * Sqrt(2 * (17 + 2 * Sqrt(2)) / 281),
         Sqrt((59 + 40 * Sqrt(2)) / 281)});
    EXPECT_THAT(actual, AlmostEquals(expected, 74));
    EXPECT_THAT(InnerProduct(actual, displacement),
                VanishesBefore(displacement.Norm(), 17));
    EXPECT_THAT(InnerProduct(actual, body_.polar_axis()),
                VanishesBefore(1, 25));
  }
  {
    auto const actual = rotation(Vector<double, TestFrame>({0, 0, 1}));
    auto const expected = Vector<double, ICRFJ2000Equator>(
        {-2 / Sqrt(29), -3 / Sqrt(29), -4 / Sqrt(29)});
    EXPECT_THAT(actual, AlmostEquals(expected, 76));
  }
}

}  // namespace internal_body_surface_dynamic_frame
}  // namespace physics
}  // namespace principia
