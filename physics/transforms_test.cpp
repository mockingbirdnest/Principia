#include "physics/transforms.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/degrees_of_freedom.hpp"
#include "physics/frame.hpp"
#include "physics/massive_body.hpp"
#include "physics/massless_body.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/numerics.hpp"

using principia::quantities::Length;
using principia::quantities::Mass;
using principia::quantities::SIUnit;
using principia::quantities::Speed;
using principia::quantities::Time;
using principia::testing_utilities::AbsoluteError;
using principia::testing_utilities::AlmostEquals;
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

  void SetUp() override {
    body1_ = std::make_unique<MassiveBody>(1 * SIUnit<Mass>());
    body2_ = std::make_unique<MassiveBody>(3 * SIUnit<Mass>());
    satellite1_from_ = std::make_unique<Trajectory<From>>(satellite1_);
    satellite2_from_ = std::make_unique<Trajectory<From>>(satellite2_);
    satellite3_from_ = std::make_unique<Trajectory<From>>(satellite3_);
    body1_from_ = std::make_unique<Trajectory<From>>(*body1_);
    body2_from_ = std::make_unique<Trajectory<From>>(*body2_);
    body1_to_ = std::make_unique<Trajectory<To>>(*body1_);
    body2_to_ = std::make_unique<Trajectory<To>>(*body2_);

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
                      Velocity<From>({4 * i * SIUnit<Speed>(),
                                      8 * i * SIUnit<Speed>(),
                                      -16 * i * SIUnit<Speed>()})));
      satellite1_from_->Append(
          Instant(i * SIUnit<Time>()),
                  DegreesOfFreedom<From>(
                      Position<From>(
                          Displacement<From>({10 * i * SIUnit<Length>(),
                                              -20 * i * SIUnit<Length>(),
                                              30 * i * SIUnit<Length>()})),
                      Velocity<From>({40 * i * SIUnit<Speed>(),
                                      -80 * i * SIUnit<Speed>(),
                                      160 * i * SIUnit<Speed>()})));
      // At the midpoint between body1_ and body2_.
      satellite2_from_->Append(
          Instant(i * SIUnit<Time>()),
                  DegreesOfFreedom<From>(
                      Position<From>(
                          Displacement<From>({0 * i * SIUnit<Length>(),
                                              0 * i * SIUnit<Length>(),
                                              6 * i * SIUnit<Length>()})),
                      Velocity<From>({8 * i * SIUnit<Speed>(),
                                      16 * i * SIUnit<Speed>(),
                                      0 * i * SIUnit<Speed>()})));
    }
  }

  std::unique_ptr<MassiveBody> body1_;
  std::unique_ptr<MassiveBody> body2_;
  MasslessBody satellite1_;
  MasslessBody satellite2_;
  MasslessBody satellite3_;
  std::unique_ptr<Trajectory<From>> body1_from_;
  std::unique_ptr<Trajectory<From>> body2_from_;
  std::unique_ptr<Trajectory<To>> body1_to_;
  std::unique_ptr<Trajectory<To>> body2_to_;
  std::unique_ptr<Trajectory<From>> satellite1_from_;
  std::unique_ptr<Trajectory<From>> satellite2_from_;
  std::unique_ptr<Trajectory<From>> satellite3_from_;

  std::unique_ptr<Transforms<From, Through, To>> transforms_;
};

TEST_F(TransformsTest, BodyCentredNonRotating) {
  transforms_ = Transforms<From, Through, To>::BodyCentredNonRotating(
                    *body1_from_, *body1_to_);
  Trajectory<Through> body1_through(*body1_);

  int i = 1;
  for (auto it = transforms_->first(satellite1_from_.get());
       !it.at_end();
       ++it, ++i) {
    DegreesOfFreedom<Through> const degrees_of_freedom =
        it.degrees_of_freedom();
    EXPECT_THAT(degrees_of_freedom.position - Position<Through>(),
                Eq(Displacement<Through>({9 * i * SIUnit<Length>(),
                                          -22 * i * SIUnit<Length>(),
                                          27 * i * SIUnit<Length>()}))) << i;
    EXPECT_THAT(degrees_of_freedom.velocity,
                Eq(Velocity<Through>({36 * i * SIUnit<Speed>(),
                                      -88 * i * SIUnit<Speed>(),
                                      144 * i * SIUnit<Speed>()}))) << i;
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
    EXPECT_THAT(degrees_of_freedom.position - Position<To>(),
                Eq(Displacement<To>({12 * i * SIUnit<Length>(),
                                     -21 * i * SIUnit<Length>(),
                                     29 * i * SIUnit<Length>()}))) << i;
    EXPECT_THAT(degrees_of_freedom.velocity,
                Eq(Velocity<To>({36 * i * SIUnit<Speed>(),
                                 -88 * i * SIUnit<Speed>(),
                                 144 * i * SIUnit<Speed>()}))) << i;
  }
}

TEST_F(TransformsTest, SatelliteBarycentricRotating) {
  transforms_ = Transforms<From, Through, To>::BarycentricRotating(
                    *body1_from_, *body1_to_,
                    *body2_from_, *body2_to_);
  Trajectory<Through> satellite1_through(satellite1_);

  int i = 1;
  for (auto it = transforms_->first(satellite1_from_.get());
       !it.at_end();
       ++it, ++i) {
    DegreesOfFreedom<Through> const degrees_of_freedom =
        it.degrees_of_freedom();
    EXPECT_THAT(degrees_of_freedom.position - Position<Through>(),
                AlmostEquals(Displacement<Through>(
                                 {-5.5 * sqrt(5.0) * i * SIUnit<Length>(),
                                  27 * i * SIUnit<Length>(),
                                  8 * sqrt(5.0) * i * SIUnit<Length>()}))) << i;
    EXPECT_THAT(degrees_of_freedom.velocity,
                AlmostEquals(Velocity<Through>(
                                 {-28 * sqrt(5.0) * i * SIUnit<Speed>(),
                                  168 * i * SIUnit<Speed>(),
                                  32 * sqrt(5.0) * i * SIUnit<Speed>()}))) << i;
    satellite1_through.Append(Instant(i * SIUnit<Time>()), degrees_of_freedom);
  }

  i = 1;
  for (auto it = transforms_->second(&satellite1_through);
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
    EXPECT_THAT(degrees_of_freedom.position - Position<To>(),
                AlmostEquals(Displacement<To>(
                                 {30 * i * SIUnit<Length>(),
                                  10 * i * SIUnit<Length>(),
                                  -20 * i * SIUnit<Length>()}))) << i;
    EXPECT_THAT(degrees_of_freedom.velocity,
                AlmostEquals(Velocity<To>(
                                 {168 * i * SIUnit<Speed>(),
                                  36 * i * SIUnit<Speed>(),
                                  -88 * i * SIUnit<Speed>()}))) << i;
  }
}

TEST_F(TransformsTest, BodiesBarycentricRotating) {
  // TODO(phl): We need two transforms to avoid getting confused by the cache.
  auto transforms1 = Transforms<From, Through, To>::BarycentricRotating(
                         *body1_from_, *body1_to_,
                         *body2_from_, *body2_to_);
  auto transforms2 = Transforms<From, Through, To>::BarycentricRotating(
                         *body1_from_, *body1_to_,
                         *body2_from_, *body2_to_);
  Trajectory<Through> body1_through(*body1_);
  Trajectory<Through> body2_through(*body2_);

  int i = 1;
  for (auto it1 = transforms1->first(body1_from_.get()),
            it2 = transforms2->first(body2_from_.get());
       !it1.at_end() && !it2.at_end();
       ++it1, ++it2, ++i) {
    DegreesOfFreedom<Through> const degrees_of_freedom1 =
        it1.degrees_of_freedom();
    DegreesOfFreedom<Through> const degrees_of_freedom2 =
        it2.degrees_of_freedom();

    EXPECT_THAT(
        AbsoluteError(degrees_of_freedom1.position - Position<Through>(),
                      Displacement<Through>(
                          {1.5 * sqrt(5.0) * i * SIUnit<Length>(),
                           0 * SIUnit<Length>(),
                           0 * SIUnit<Length>()})),
        Lt(40 * DBL_EPSILON * SIUnit<Length>()));
    EXPECT_THAT(
        AbsoluteError(degrees_of_freedom2.position - Position<Through>(),
                      Displacement<Through>(
                          {-0.5 * sqrt(5.0) * i * SIUnit<Length>(),
                           0 * SIUnit<Length>(),
                           0 * SIUnit<Length>()})),
        Lt(50 * DBL_EPSILON * SIUnit<Length>()));
    // TODO(phl): Add a test that the velocities are along X once we have fixed
    // them.

    DegreesOfFreedom<Through> const barycenter_degrees_of_freedom =
        Barycentre<Through, Mass>({degrees_of_freedom1, degrees_of_freedom2},
                                  {body1_->mass(), body2_->mass()});
    EXPECT_THAT(AbsoluteError(barycenter_degrees_of_freedom.position -
                                  Position<Through>(),
                              Displacement<Through>()),
                Lt(40 * DBL_EPSILON * SIUnit<Length>())) << i;
    EXPECT_THAT(AbsoluteError(barycenter_degrees_of_freedom.velocity,
                              Velocity<Through>()),
                Lt(140 * DBL_EPSILON * SIUnit<Speed>())) << i;

    Length const length = (degrees_of_freedom1.position -
                           degrees_of_freedom2.position).Norm();
    EXPECT_THAT(length, AlmostEquals(2.0 * sqrt(5.0) * i * SIUnit<Length>()));
  }
}

}  // namespace physics
}  // namespace principia
