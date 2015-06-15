#include "physics/transformz.hpp"

#include <limits>

#include "geometry/frame.hpp"
#include "base/not_null.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/degrees_of_freedom.hpp"
#include "physics/massive_body.hpp"
#include "physics/massless_body.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {

using base::make_not_null_unique;
using geometry::Frame;
using geometry::InnerProduct;
using geometry::Normalize;
using geometry::Rotation;
using quantities::Abs;
using quantities::Length;
using quantities::Mass;
using quantities::Speed;
using quantities::Time;
using si::Kilogram;
using si::Metre;
using si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using testing_utilities::VanishesBefore;
using ::testing::Eq;
using ::testing::Lt;

namespace physics {

namespace {
const int kNumberOfPoints = 33;
const Time kStep = 1 * Second;
const Length kLowTolerance = 0.001 * Metre;
const Length kHighTolerance = 0.01 * Metre;
}  // namespace

class TransformzTest : public testing::Test {
 protected:
  using From = Frame<serialization::Frame::TestTag,
                     serialization::Frame::FROM, true>;
  using Through = Frame<serialization::Frame::TestTag,
                        serialization::Frame::THROUGH, false>;
  using To = Frame<serialization::Frame::TestTag,
                   serialization::Frame::TO, true>;

  static void SetUpTestCase() {
    // This will need to change if we ever have continuous trajectories with
    // more than 8 divisions.
    CHECK_EQ(1, kNumberOfPoints % 8);
  }

  TransformzTest()
      : body1_(MassiveBody(1 * Kilogram)),
        body2_(MassiveBody(3 * Kilogram)),
        body1_from_(kStep, kLowTolerance, kHighTolerance),
        body1_to_(kStep, kLowTolerance, kHighTolerance),
        body2_from_(kStep, kLowTolerance, kHighTolerance),
        body2_to_(kStep, kLowTolerance, kHighTolerance),
        satellite_from_(make_not_null_unique<Trajectory<From>>(&satellite_)),
        satellite_through_(
            make_not_null_unique<Trajectory<Through>>(&satellite_)),
        satellite_to_(make_not_null_unique<Trajectory<To>>(&satellite_)) {
    // The various bodies move have both a position and a velocity that
    // increases linearly with time.  This is not a situation that's physically
    // possible, but we don't care, all we want is to make sure that the
    // transforms are properly performed: |Transforms| doesn't know anything
    // about physics.  Also, the trajectories were chosen so that we are not in
    // any "special case" with respect to the positions or the velocities.
    for (int i = 1; i <= kNumberOfPoints; ++i) {
      body1_from_.Append(Instant(i * Second),
                         DegreesOfFreedom<From>(
                             Position<From>(
                                 Displacement<From>({1 * i * Metre,
                                                     2 * i * Metre,
                                                     3 * i * Metre})),
                             Velocity<From>({4 * i * Metre / Second,
                                             8 * i * Metre / Second,
                                             16 * i * Metre / Second})));
      body1_to_.Append(Instant(i * Second),
                       DegreesOfFreedom<To>(
                           Position<To>(
                               Displacement<To>({3 * i * Metre,
                                                 1 * i * Metre,
                                                 2 * i * Metre})),
                           Velocity<To>({16 * i * Metre / Second,
                                         4 * i * Metre / Second,
                                         8 * i * Metre / Second})));
      body2_from_.Append(Instant(i * Second),
                         DegreesOfFreedom<From>(
                             Position<From>(
                                 Displacement<From>({-1 * i * Metre,
                                                     -2 * i * Metre,
                                                     3 * i * Metre})),
                             Velocity<From>({-4 * i * Metre / Second,
                                             8 * i * Metre / Second,
                                             -16 * i * Metre / Second})));
      body2_to_.Append(Instant(i * Second),
                       DegreesOfFreedom<To>(
                           Position<To>(
                               Displacement<To>({3 * i * Metre,
                                                -1 * i * Metre,
                                                -2 * i * Metre})),
                           Velocity<To>({-16 * i * Metre / Second,
                                         4 * i * Metre / Second,
                                         8 * i * Metre / Second})));
      satellite_from_->Append(Instant(i * Second),
                              DegreesOfFreedom<From>(
                                  Position<From>(
                                      Displacement<From>({10 * i * Metre,
                                                          -20 * i * Metre,
                                                          30 * i * Metre})),
                                  Velocity<From>({40 * i * Metre / Second,
                                                  -80 * i * Metre / Second,
                                                  160 * i * Metre / Second})));
    }
  }

  MassiveBody body1_;
  MassiveBody body2_;
  MasslessBody satellite_;
  ContinuousTrajectory<From> body1_from_;
  ContinuousTrajectory<To> body1_to_;
  ContinuousTrajectory<From> body2_from_;
  ContinuousTrajectory<To> body2_to_;
  not_null<std::unique_ptr<Trajectory<From>>> satellite_from_;
  not_null<std::unique_ptr<Trajectory<Through>>> satellite_through_;
  not_null<std::unique_ptr<Trajectory<To>>> satellite_to_;
};

// This transform is simple enough that we can compute its effect by hand.  This
// test verifies that we get the expected result both in |Through| and in |To|.
TEST_F(TransformzTest, BodyCentredNonRotating) {
  auto const transforms =
      Transformz<From, Through, To>::BodyCentredNonRotating(
          body1_, body1_from_, body1_to_);

  int i = 1;
  for (auto it = transforms->first(*satellite_from_);
       !it.at_end();
       ++it, ++i) {
    DegreesOfFreedom<Through> const degrees_of_freedom =
        it.degrees_of_freedom();
    EXPECT_THAT(degrees_of_freedom,
                Componentwise(
                    AlmostEquals(Through::origin + Displacement<Through>(
                        {9 * i * Metre,
                         -22 * i * Metre,
                         27 * i * Metre}),
                        0, 384),
                    AlmostEquals(Velocity<Through>(
                        {36 * i * Metre / Second,
                         -88 * i * Metre / Second,
                         144 * i * Metre / Second}),
                        0, 720))) << i;
    satellite_through_->Append(Instant(i * Second), degrees_of_freedom);
  }

  i = 1;
  for (auto it = transforms->second(Instant(kNumberOfPoints * Second),
                                    *satellite_through_);
       !it.at_end();
       ++it, ++i) {
    DegreesOfFreedom<To> const degrees_of_freedom =
        it.degrees_of_freedom();
    EXPECT_THAT(degrees_of_freedom,
                Componentwise(
                    AlmostEquals(To::origin + Displacement<To>(
                        {(99 + 9 * i) * Metre,
                         (33 -22 * i) * Metre,
                         (66 + 27 * i) * Metre}),
                        2, 1036),
                    AlmostEquals(Velocity<To>(
                        {36 * i * Metre / Second,
                         -88 * i * Metre / Second,
                         144 * i * Metre / Second}),
                        0, 720))) << i;
  }
}

// Check that the computations we do match those done using Mathematica.
TEST_F(TransformzTest, SatelliteBarycentricRotating) {
  auto const transforms =
      Transformz<From, Through, To>::BarycentricRotating(
          body1_, body1_from_, body1_to_,
          body2_, body2_from_, body2_to_);
  Trajectory<Through> satellite_through(&satellite_);

  int i = 1;
  for (auto it = transforms->first(*satellite_from_);
       !it.at_end();
       ++it, ++i) {
    DegreesOfFreedom<Through> const degrees_of_freedom =
        it.degrees_of_freedom();
    EXPECT_THAT(degrees_of_freedom,
                Componentwise(
                    AlmostEquals(Through::origin + Displacement<Through>(
                        {-5.5 * sqrt(5.0) * i * Metre,
                         62.0 * sqrt(5.0 / 21.0) * i * Metre,
                         53.0 / sqrt(21.0) * i * Metre}),
                        2, 4607),
                    AlmostEquals(Velocity<Through>(
                        {(362.0 / sqrt(5.0)) * i * Metre / Second,
                         (2776.0 / sqrt(105.0)) * i * Metre / Second,
                         176.0 / sqrt(21.0) * i * Metre / Second}),
                        8, 51340))) << i;
    satellite_through.Append(Instant(i * Second), degrees_of_freedom);
  }

  i = 1;
  for (auto it = transforms->second(Instant(kNumberOfPoints * Second),
                                    satellite_through);
       !it.at_end();
       ++it, ++i) {
    DegreesOfFreedom<To> const degrees_of_freedom =
        it.degrees_of_freedom();
    EXPECT_THAT(degrees_of_freedom,
                Componentwise(
                    AlmostEquals(To::origin + Displacement<To>(
                        {(99.0 + (62.0 * sqrt(5.0 / 21.0)) * i) * Metre,
                         (-16.5 + (-5.5 + 106.0 / sqrt(105.0)) * i) * Metre,
                         (-33.0 + (-11.0 - 53.0 / sqrt(105.0)) * i) * Metre}),
                        129, 28312),
                    AlmostEquals(Velocity<To>(
                        {2776.0 / sqrt(105.0) * i * Metre / Second,
                         (72.4 + 352.0 / sqrt(105.0)) * i * Metre / Second,
                         (144.8 - 176.0 / sqrt(105.0)) * i * Metre / Second}),
                        349, 22568))) << i;
  }
}

// Check that the bodies remain on the invariant line and at the right distance
// from each other and from the barycentre, and that the barycentre is the
// centre of the coordinates.
TEST_F(TransformzTest, BodiesBarycentricRotating) {
  auto const transforms =
      Transformz<From, Through, To>::BarycentricRotating(
          body1_, body1_from_, body1_to_,
          body2_, body2_from_, body2_to_);

  // Compute discrete trajectories from the continuous ones since we can only
  // transform discrete trajectories.
  Trajectory<From> discrete_body1_from(&body1_);
  Trajectory<From> discrete_body2_from(&body2_);
  for (int i = 1; i <= kNumberOfPoints; ++i) {
    discrete_body1_from.Append(
        Instant(i * Second),
        body1_from_.EvaluateDegreesOfFreedom(Instant(i * Second),
                                             nullptr /*hint*/));
    discrete_body2_from.Append(
        Instant(i * Second),
        body2_from_.EvaluateDegreesOfFreedom(Instant(i * Second),
                                             nullptr /*hint*/));
  }

  int i = 1;
  for (auto it1 = transforms->first(discrete_body1_from),
            it2 = transforms->first(discrete_body2_from);
       !it1.at_end() && !it2.at_end();
       ++it1, ++it2, ++i) {
    Length const l = i * Metre;
    Speed const s = i * Metre / Second;
    DegreesOfFreedom<Through> const degrees_of_freedom1 =
        it1.degrees_of_freedom();
    DegreesOfFreedom<Through> const degrees_of_freedom2 =
        it2.degrees_of_freedom();

    EXPECT_THAT(degrees_of_freedom1.position() - Through::origin,
                Componentwise(AlmostEquals(1.5 * sqrt(5.0) * l, 0, 2403),
                              VanishesBefore(l, 0, 12),
                              VanishesBefore(l, 0, 4)));
    EXPECT_THAT(degrees_of_freedom2.position() - Through::origin,
                Componentwise(AlmostEquals(-0.5 * sqrt(5.0) * l, 1, 1602),
                              VanishesBefore(l, 0, 6),
                              VanishesBefore(l, 0, 2)));
    EXPECT_THAT(degrees_of_freedom1.velocity(),
                Componentwise(AlmostEquals(6.0 / sqrt(5.0) * s, 6, 40219),
                              VanishesBefore(s, 0, 80),
                              VanishesBefore(s, 0, 10)));
    EXPECT_THAT(degrees_of_freedom2.velocity(),
                Componentwise(AlmostEquals(-2.0 / sqrt(5.0) * s, 4, 53633),
                              VanishesBefore(s, 0, 20),
                              VanishesBefore(s, 0, 13)));

    DegreesOfFreedom<Through> const barycentre_degrees_of_freedom =
        Barycentre<Through, Mass>({degrees_of_freedom1, degrees_of_freedom2},
                                  {body1_.mass(), body2_.mass()});
    EXPECT_THAT(barycentre_degrees_of_freedom.position() - Through::origin,
                Componentwise(VanishesBefore(l, 0, 2),
                              VanishesBefore(l, 0, 4),
                              VanishesBefore(l, 0, 2)));
    EXPECT_THAT(barycentre_degrees_of_freedom.velocity(),
                Componentwise(VanishesBefore(s, 0, 35),
                              VanishesBefore(s, 0, 18),
                              VanishesBefore(s, 0, 9)));

    Length const length = (degrees_of_freedom1.position() -
                           degrees_of_freedom2.position()).Norm();
    EXPECT_THAT(length,
                AlmostEquals(2.0 * sqrt(5.0) * i * Metre, 0, 1602));
  }
}

TEST_F(TransformzTest, CoordinateFrame) {
  Instant const t(kNumberOfPoints * Second);
  {
    auto const transforms =
        Transformz<From, Through, To>::BodyCentredNonRotating(
            body1_, body1_from_, body1_to_);
    auto const identity = Rotation<To, To>::Identity();
    EXPECT_EQ(identity.quaternion(),
              transforms->coordinate_frame(t)(To::origin).quaternion());
  }
  {
    auto const transforms =
        Transformz<From, Through, To>::BarycentricRotating(
            body1_, body1_from_, body1_to_,
            body2_, body2_from_, body2_to_);

    Vector<double, To> const x({1, 0, 0});
    Vector<double, To> const y({0, 1, 0});
    Vector<double, To> const z({0, 0, 1});
    Vector<double, To> const body1_body2 =
        Normalize(
            body1_to_.EvaluatePosition(t, nullptr /*hint*/) -
            body2_to_.EvaluatePosition(t, nullptr /*hint*/));
    EXPECT_THAT(
        transforms->coordinate_frame(t)(To::origin)(x),
        Componentwise(VanishesBefore(1, 666),
                      AlmostEquals(body1_body2.coordinates().y, 1),
                      AlmostEquals(body1_body2.coordinates().z, 1)));
    EXPECT_GT(
        InnerProduct(
            transforms->coordinate_frame(t)(To::origin)(y),
            Normalize(body1_to_.EvaluateVelocity(t, nullptr /*hint*/))),
        0);
    EXPECT_THAT(
        InnerProduct(
                transforms->coordinate_frame(t)(To::origin)(z),
                Normalize(
                    body1_to_.EvaluateVelocity(t, nullptr /*hint*/))),
        VanishesBefore(1, 0));
  }
}

}  // namespace physics
}  // namespace principia
