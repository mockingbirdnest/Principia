#include "physics/transforms.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/degrees_of_freedom.hpp"
#include "physics/frame.hpp"
#include "physics/massive_body.hpp"
#include "physics/massless_body.hpp"
#include "quantities/named_quantities.hpp"
#include "quantities/quantities.hpp"

using principia::quantities::Length;
using principia::quantities::Mass;
using principia::quantities::SIUnit;
using principia::quantities::Speed;
using principia::quantities::Time;
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
    body2_ = std::make_unique<MassiveBody>(2 * SIUnit<Mass>());
    body_from_ = std::make_unique<Trajectory<From>>(body_);
    body1_from_ = std::make_unique<Trajectory<From>>(*body1_);
    body2_from_ = std::make_unique<Trajectory<From>>(*body2_);
    body1_to_ = std::make_unique<Trajectory<To>>(*body1_);
    body2_to_ = std::make_unique<Trajectory<To>>(*body2_);

    for (int i = 0; i < kNumberOfPoints; ++i) {
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
      body1_to_->Append(
          Instant(i * SIUnit<Time>()),
                  DegreesOfFreedom<To>(
                      Position<To>(
                          Displacement<To>({3 * i * SIUnit<Length>(),
                                            1 * i * SIUnit<Length>(),
                                            2 * i * SIUnit<Length>()})),
                      Velocity<To>({16 * i * SIUnit<Speed>(),
                                    8 * i * SIUnit<Speed>(),
                                    4 * i * SIUnit<Speed>()})));
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
  int i = 0;
  for (auto it = transforms_->first(body_from_.get());
       !it.at_end();
       ++it, ++i) {
    DegreesOfFreedom<Through> const degrees_of_freedom =
        it.degrees_of_freedom();
    EXPECT_THAT(degrees_of_freedom.position,
                Eq(Position<Through>(
                       Displacement<Through>({9 * i * SIUnit<Length>(),
                                              -22 * i * SIUnit<Length>(),
                                              27 * i * SIUnit<Length>()}))));
    EXPECT_THAT(degrees_of_freedom.velocity,
                Eq(Velocity<Through>({36 * i * SIUnit<Speed>(),
                                      -88 * i * SIUnit<Speed>(),
                                      144 * i * SIUnit<Speed>()})));
  }
}

}  // namespace physics
}  // namespace principia
