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

using principia::quantities::Length;
using principia::quantities::Mass;
using principia::quantities::SIUnit;
using principia::quantities::Speed;
using principia::quantities::Time;
using principia::testing_utilities::AlmostEquals;
using testing::Eq;

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
    body_from_ = std::make_unique<Trajectory<From>>(body_);
    body1_from_ = std::make_unique<Trajectory<From>>(*body1_);
    body2_from_ = std::make_unique<Trajectory<From>>(*body2_);
    body1_to_ = std::make_unique<Trajectory<To>>(*body1_);
    body2_to_ = std::make_unique<Trajectory<To>>(*body2_);

    for (int i = 1; i <= kNumberOfPoints; ++i) {
      body_from_->Append(
          Instant(i * SIUnit<Time>()),
                  DegreesOfFreedom<From>(
                      Position<From>(
                          Displacement<From>({10 * i * SIUnit<Length>(),
                                              -20 * i * SIUnit<Length>(),
                                              30 * i * SIUnit<Length>()})),
                      Velocity<From>({40 * i * SIUnit<Speed>(),
                                      -80 * i * SIUnit<Speed>(),
                                      160 * i * SIUnit<Speed>()})));
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
    }
  }

  MasslessBody body_;
  std::unique_ptr<MassiveBody> body1_;
  std::unique_ptr<MassiveBody> body2_;
  std::unique_ptr<Trajectory<From>> body_from_;
  std::unique_ptr<Trajectory<From>> body1_from_;
  std::unique_ptr<Trajectory<From>> body2_from_;
  std::unique_ptr<Trajectory<To>> body1_to_;
  std::unique_ptr<Trajectory<To>> body2_to_;

  std::unique_ptr<Transforms<From, Through, To>> transforms_;
};

TEST_F(TransformsTest, BodyCentredNonRotating) {
  transforms_ = Transforms<From, Through, To>::BodyCentredNonRotating(
                    *body1_from_, *body1_to_);
  Trajectory<Through> body1_through(*body1_);

  int i = 1;
  for (auto it = transforms_->first(body_from_.get());
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

TEST_F(TransformsTest, BarycentricRotating) {
  transforms_ = Transforms<From, Through, To>::BarycentricRotating(
                    *body1_from_, *body1_to_,
                    *body2_from_, *body2_to_);
  Trajectory<Through> body1_through(*body1_);

  int i = 1;
  for (auto it = transforms_->first(body_from_.get());
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

}  // namespace physics
}  // namespace principia
