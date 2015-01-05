#include "physics/transforms.hpp"

#include <limits>

#include "base/not_null.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/degrees_of_freedom.hpp"
#include "physics/frame.hpp"
#include "physics/massive_body.hpp"
#include "physics/massless_body.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/vanishes_before.hpp"

using principia::base::check_not_null;
using principia::quantities::Length;
using principia::quantities::Mass;
using principia::quantities::SIUnit;
using principia::quantities::Speed;
using principia::quantities::Time;
using principia::testing_utilities::AlmostEquals;
using principia::testing_utilities::Componentwise;
using principia::testing_utilities::VanishesBefore;
using testing::Eq;
using testing::Lt;

namespace principia {
namespace physics {

namespace {
const int kNumberOfPoints = 20;
}  // namespace

class TransformsTest : public testing::Test {
 protected:
  enum class Tag {
    kFrom,
    kThrough,
    kTo,
  };
  using From = Frame<Tag, Tag::kFrom, true>;
  using Through = Frame<Tag, Tag::kThrough, false>;
  using To = Frame<Tag, Tag::kTo, true>;

  TransformsTest()
      : body1_(MassiveBody(1 * SIUnit<Mass>())),
        body2_(MassiveBody(3 * SIUnit<Mass>())) {
    satellite_from_ = std::make_unique<Trajectory<From>>(
        check_not_null(&satellite_));
    body1_from_ = std::make_unique<Trajectory<From>>(check_not_null(&body1_));
    body2_from_ = std::make_unique<Trajectory<From>>(check_not_null(&body2_));
    body1_to_ = std::make_unique<Trajectory<To>>(check_not_null(&body1_));
    body2_to_ = std::make_unique<Trajectory<To>>(check_not_null(&body2_));
    body1_from_fn_ =
        [this]() -> Trajectory<From> const& { return *this->body1_from_; };
    body2_from_fn_ =
        [this]() -> Trajectory<From> const& { return *this->body2_from_; };
    body1_to_fn_ =
        [this]() -> Trajectory<To> const& { return *this->body1_to_; };
    body2_to_fn_ =
        [this]() -> Trajectory<To> const& { return *this->body2_to_; };

    // The various bodies move have both a position and a velocity that
    // increases linearly with time.  This is not a situation that's physically
    // possible, but we don't care, all we want is to make sure that the
    // transforms are properly performed: |Transforms| doesn't know anything
    // about physics.  Also, the trajectories were chosen so that we are not in
    // any "special case" with respect to the positions or the velocities.
    for (int i = 1; i <= kNumberOfPoints; ++i) {
      body1_from_->Append(
          Instant(i * SIUnit<Time>()),
                  DegreesOfFreedom<From>(
                      Position<From>(
                          Displacement<From>({1 * i * SIUnit<Length>(),
                                              2 * i * SIUnit<Length>(),
                                              3 * i * SIUnit<Length>()})),
                      Velocity<From>({4 * i * SIUnit<Speed>(),
                                      8 * i * SIUnit<Speed>(),
                                      16 * i * SIUnit<Speed>()})));
      body2_from_->Append(
          Instant(i * SIUnit<Time>()),
                  DegreesOfFreedom<From>(
                      Position<From>(
                          Displacement<From>({-1 * i * SIUnit<Length>(),
                                              -2 * i * SIUnit<Length>(),
                                              3 * i * SIUnit<Length>()})),
                      Velocity<From>({-4 * i * SIUnit<Speed>(),
                                      8 * i * SIUnit<Speed>(),
                                      -16 * i * SIUnit<Speed>()})));
      satellite_from_->Append(
          Instant(i * SIUnit<Time>()),
                  DegreesOfFreedom<From>(
                      Position<From>(
                          Displacement<From>({10 * i * SIUnit<Length>(),
                                              -20 * i * SIUnit<Length>(),
                                              30 * i * SIUnit<Length>()})),
                      Velocity<From>({40 * i * SIUnit<Speed>(),
                                      -80 * i * SIUnit<Speed>(),
                                      160 * i * SIUnit<Speed>()})));
    }
  }

  MassiveBody body1_;
  MassiveBody body2_;
  MasslessBody satellite_;
  std::unique_ptr<Trajectory<From>> body1_from_;
  std::unique_ptr<Trajectory<From>> body2_from_;
  std::unique_ptr<Trajectory<To>> body1_to_;
  std::unique_ptr<Trajectory<To>> body2_to_;
  Transforms<From, Through, To>::LazyTrajectory<From> body1_from_fn_;
  Transforms<From, Through, To>::LazyTrajectory<From> body2_from_fn_;
  Transforms<From, Through, To>::LazyTrajectory<To> body1_to_fn_;
  Transforms<From, Through, To>::LazyTrajectory<To> body2_to_fn_;
  std::unique_ptr<Trajectory<From>> satellite_from_;

  std::unique_ptr<Transforms<From, Through, To>> transforms_;
};

// This transform is simple enough that we can compute its effect by hand.  This
// test verifies that we get the expected result both in |Through| and in |To|.
TEST_F(TransformsTest, BodyCentredNonRotating) {
  transforms_ = Transforms<From, Through, To>::BodyCentredNonRotating(
                    body1_from_fn_, body1_to_fn_);
  Trajectory<Through> body1_through(check_not_null(&body1_));

  int i = 1;
  for (auto it = transforms_->first(satellite_from_.get());
       !it.at_end();
       ++it, ++i) {
    DegreesOfFreedom<Through> const degrees_of_freedom =
        it.degrees_of_freedom();
    EXPECT_THAT(degrees_of_freedom,
                Componentwise(
                    Eq(Through::origin +
                       Displacement<Through>({9 * i * SIUnit<Length>(),
                                              -22 * i * SIUnit<Length>(),
                                              27 * i * SIUnit<Length>()})),
                    Eq(Velocity<Through>({36 * i * SIUnit<Speed>(),
                                          -88 * i * SIUnit<Speed>(),
                                          144 * i * SIUnit<Speed>()})))) << i;
    body1_through.Append(Instant(i * SIUnit<Time>()), degrees_of_freedom);
  }

  i = 1;
  for (auto it = transforms_->second(&body1_through);
       !it.at_end();
       ++it, ++i) {
    body1_to_->Append(
        Instant(i * SIUnit<Time>()),
                DegreesOfFreedom<To>(
                    Position<To>(
                        Displacement<To>({3 * i * SIUnit<Length>(),
                                          1 * i * SIUnit<Length>(),
                                          2 * i * SIUnit<Length>()})),
                    Velocity<To>({16 * i * SIUnit<Speed>(),
                                  4 * i * SIUnit<Speed>(),
                                  8 * i * SIUnit<Speed>()})));

    DegreesOfFreedom<To> const degrees_of_freedom =
        it.degrees_of_freedom();
    EXPECT_THAT(degrees_of_freedom,
                Componentwise(
                    Eq(To::origin +
                        Displacement<To>({12 * i * SIUnit<Length>(),
                                          -21 * i * SIUnit<Length>(),
                                          29 * i * SIUnit<Length>()})),
                    Eq(Velocity<To>({36 * i * SIUnit<Speed>(),
                                      -88 * i * SIUnit<Speed>(),
                                      144 * i * SIUnit<Speed>()})))) << i;
  }
}

// Check that the computations we do match those done using Mathematica.
TEST_F(TransformsTest, SatelliteBarycentricRotating) {
  transforms_ = Transforms<From, Through, To>::BarycentricRotating(
                    body1_from_fn_, body1_to_fn_,
                    body2_from_fn_, body2_to_fn_);
  Trajectory<Through> satellite_through(check_not_null(&satellite_));

  int i = 1;
  for (auto it = transforms_->first(satellite_from_.get());
       !it.at_end();
       ++it, ++i) {
    DegreesOfFreedom<Through> const degrees_of_freedom =
        it.degrees_of_freedom();
    EXPECT_THAT(degrees_of_freedom.position() - Position<Through>(),
                AlmostEquals(Displacement<Through>(
                    {-5.5 * sqrt(5.0) * i * SIUnit<Length>(),
                     62.0 * sqrt(5.0 / 21.0) * i * SIUnit<Length>(),
                     53.0 / sqrt(21.0) * i * SIUnit<Length>()}),
                    1, 8)) << i;
    EXPECT_THAT(degrees_of_freedom.velocity(),
                AlmostEquals(Velocity<Through>(
                    {(362.0 / sqrt(5.0)) * i * SIUnit<Speed>(),
                     (2776.0 / sqrt(105.0)) * i * SIUnit<Speed>(),
                     176.0 / sqrt(21.0) * i * SIUnit<Speed>()}),
                    1, 10)) << i;
    satellite_through.Append(Instant(i * SIUnit<Time>()), degrees_of_freedom);
  }

  i = 1;
  for (auto it = transforms_->second(&satellite_through);
       !it.at_end();
       ++it, ++i) {
      body1_to_->Append(
          Instant(i * SIUnit<Time>()),
                  DegreesOfFreedom<To>(
                      Position<To>(
                          Displacement<To>({3 * i * SIUnit<Length>(),
                                            1 * i * SIUnit<Length>(),
                                            2 * i * SIUnit<Length>()})),
                      Velocity<To>({16 * i * SIUnit<Speed>(),
                                    4 * i * SIUnit<Speed>(),
                                    8 * i * SIUnit<Speed>()})));
      body2_to_->Append(
          Instant(i * SIUnit<Time>()),
                  DegreesOfFreedom<To>(
                      Position<To>(
                          Displacement<To>({3 * i * SIUnit<Length>(),
                                           -1 * i * SIUnit<Length>(),
                                           -2 * i * SIUnit<Length>()})),
                      Velocity<To>({-16 * i * SIUnit<Speed>(),
                                    4 * i * SIUnit<Speed>(),
                                    8 * i * SIUnit<Speed>()})));

    DegreesOfFreedom<To> const degrees_of_freedom =
        it.degrees_of_freedom();
    EXPECT_THAT(degrees_of_freedom.position() - To::origin,
                AlmostEquals(Displacement<To>(
                    {(3.0 + 62.0 * sqrt(5.0 / 21.0)) * i * SIUnit<Length>(),
                     (-6.0 + 106.0 / sqrt(105.0)) * i * SIUnit<Length>(),
                     (-12.0 - 53.0 / sqrt(105.0)) * i * SIUnit<Length>()}),
                    3, 21))
        << i;
    EXPECT_THAT(degrees_of_freedom.velocity(),
                AlmostEquals(Velocity<To>(
                    {2776.0 / sqrt(105.0) * i * SIUnit<Speed>(),
                     (72.4 + 352.0 / sqrt(105.0)) * i * SIUnit<Speed>(),
                     (144.8 - 176.0 / sqrt(105.0)) * i * SIUnit<Speed>()}),
                    1, 8)) << i;
  }
}

// Check that the bodies remain on the invariant line and at the right distance
// from each other and from the barycentre, and that the barycentre is the
// centre of the coordinates.
TEST_F(TransformsTest, BodiesBarycentricRotating) {
  transforms_ = Transforms<From, Through, To>::BarycentricRotating(
                    body1_from_fn_, body1_to_fn_,
                    body2_from_fn_, body2_to_fn_);
  Trajectory<Through> body1_through(check_not_null(&body1_));
  Trajectory<Through> body2_through(check_not_null(&body2_));

  int i = 1;
  for (auto it1 = transforms_->first(body1_from_.get()),
            it2 = transforms_->first(body2_from_.get());
       !it1.at_end() && !it2.at_end();
       ++it1, ++it2, ++i) {
    Length const l = i * SIUnit<Length>();
    Speed const s = i * SIUnit<Speed>();
    DegreesOfFreedom<Through> const degrees_of_freedom1 =
        it1.degrees_of_freedom();
    DegreesOfFreedom<Through> const degrees_of_freedom2 =
        it2.degrees_of_freedom();

    EXPECT_THAT(degrees_of_freedom1.position() - Through::origin,
                Componentwise(AlmostEquals(1.5 * sqrt(5.0) * l, 0, 1),
                              VanishesBefore(l, 0, 4),
                              VanishesBefore(l, 0, 2)));
    EXPECT_THAT(degrees_of_freedom2.position() - Through::origin,
                Componentwise(AlmostEquals(-0.5 * sqrt(5.0) * l, 0, 2),
                              VanishesBefore(l, 0, 2),
                              VanishesBefore(l, 0, 1)));
    EXPECT_THAT(degrees_of_freedom1.velocity(),
                Componentwise(AlmostEquals(6.0 / sqrt(5.0) * s, 0, 7),
                              VanishesBefore(s, 0, 34),
                              VanishesBefore(s, 0, 2)));
    EXPECT_THAT(degrees_of_freedom2.velocity(),
                Componentwise(AlmostEquals(-2.0 / sqrt(5.0) * s, 0, 14),
                              VanishesBefore(s, 0, 16),
                              VanishesBefore(s, 0, 1)));

    DegreesOfFreedom<Through> const barycentre_degrees_of_freedom =
        Barycentre<Through, Mass>({degrees_of_freedom1, degrees_of_freedom2},
                                  {body1_.mass(), body2_.mass()});
    EXPECT_THAT(barycentre_degrees_of_freedom.position() - Through::origin,
                Componentwise(VanishesBefore(l, 0, 1),
                              VanishesBefore(l, 0, 2),
                              VanishesBefore(l, 0, 1)));
    EXPECT_THAT(barycentre_degrees_of_freedom.velocity(),
                Componentwise(VanishesBefore(s, 0, 14),
                              VanishesBefore(s, 0, 7),
                              VanishesBefore(s, 0, 1)));

    Length const length = (degrees_of_freedom1.position() -
                           degrees_of_freedom2.position()).Norm();
    EXPECT_THAT(length,
                AlmostEquals(2.0 * sqrt(5.0) * i * SIUnit<Length>(), 0, 2));
  }
}

}  // namespace physics
}  // namespace principia
