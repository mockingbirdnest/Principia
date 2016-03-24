
#include "ksp_plugin/vessel.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/mock_ephemeris.hpp"

namespace principia {

using physics::MockEphemeris;
using quantities::si::Kilogram;
using quantities::si::Metre;
using quantities::si::Second;

namespace ksp_plugin {

// TODO(egg): We would want to use a real ephemeris to properly exercise the
// limit cases.
class VesselTest : public testing::Test {
 protected:
  VesselTest()
      : body_(1 * Kilogram),
        parent_(&body_),
        adaptive_parameters_(
            DormandElMikkawyPrince1986RKN434FM<Position<Barycentric>>(),
            /*max_steps=*/1,
            /*length_integration_tolerance=*/1 * Metre,
            /*speed_integration_tolerance=*/1 * Metre / Second),
        fixed_parameters_(
            McLachlanAtela1992Order5Optimal<Position<Barycentric>>(),
            /*step=*/1 * Second),
        vessel_(make_not_null_unique<Vessel>(&parent_,
                                             &ephemeris_,
                                             adaptive_parameters_,
                                             fixed_parameters_)) {}

  MassiveBody body_;
  Celestial parent_;
  MockEphemeris<Barycentric> ephemeris_;
  Ephemeris<Barycentric>::AdaptiveStepParameters const adaptive_parameters_;
  Ephemeris<Barycentric>::FixedStepParameters const fixed_parameters_;
  not_null<std::unique_ptr<Vessel>> vessel_;
  DegreesOfFreedom<Barycentric> d1_ =
      {Barycentric::origin +
           Displacement<Barycentric>({1 * Metre, 2 * Metre, 3 * Metre}),
       Velocity<Barycentric>({4 * Metre / Second,
                              5 * Metre / Second,
                              6 * Metre / Second})};
  DegreesOfFreedom<Barycentric> d2_ =
      {Barycentric::origin +
           Displacement<Barycentric>({11 * Metre, 12 * Metre, 13 * Metre}),
       Velocity<Barycentric>({14 * Metre / Second,
                              15 * Metre / Second,
                              16 * Metre / Second})};
  Instant const t1_ = kUniversalTimeEpoch;
  Instant const t2_ = kUniversalTimeEpoch + 42.3 * Second;
};

using VesselDeathTest = VesselTest;

TEST_F(VesselDeathTest, Uninitialized) {
  EXPECT_DEATH({
    vessel_->history();
  }, "is_initialized");
  EXPECT_DEATH({
    vessel_->prolongation();
  }, "is_initialized");
}

TEST_F(VesselTest, Initialization) {
  EXPECT_FALSE(vessel_->is_initialized());
  vessel_->CreateHistoryAndForkProlongation(t2_, d2_);
  EXPECT_TRUE(vessel_->is_initialized());
  auto const& prolongation = vessel_->prolongation();
  EXPECT_EQ(t2_, prolongation.last().time());
  auto const& history = vessel_->history();
  EXPECT_EQ(t2_, history.last().time());
  EXPECT_FALSE(vessel_->has_flight_plan());
  EXPECT_FALSE(vessel_->has_prediction());
}

TEST_F(VesselTest, Dirty) {
  EXPECT_FALSE(vessel_->is_dirty());
  vessel_->set_dirty();
  EXPECT_TRUE(vessel_->is_dirty());
}

TEST_F(VesselTest, Parent) {
  MassiveBody body(2 * Kilogram);
  Celestial celestial(&body);

  EXPECT_EQ(&parent_, vessel_->parent());
  vessel_->set_parent(&celestial);
  EXPECT_EQ(&celestial, vessel_->parent());
}

TEST_F(VesselTest, AdvanceTime) {
  vessel_->CreateHistoryAndForkProlongation(t1_, d1_);
  vessel_->AdvanceTimeNotInBubble(t2_);
  EXPECT_EQ(kUniversalTimeEpoch + 42 * Second,
            vessel_->history().last().time());
  EXPECT_EQ(t2_, vessel_->prolongation().last().time());
}

TEST_F(VesselDeathTest, SerializationError) {
  EXPECT_DEATH({
    serialization::Vessel message;
    vessel_->WriteToMessage(&message);
  }, "is_initialized");
  EXPECT_DEATH({
    serialization::Vessel message;
    Vessel::ReadFromMessage(message, &ephemeris_, &parent_);
  }, "message.has_history");
}

TEST_F(VesselTest, SerializationSuccess) {
  serialization::Vessel message;
  EXPECT_FALSE(message.has_history());
  vessel_->CreateHistoryAndForkProlongation(t2_, d2_);
  vessel_->WriteToMessage(&message);
  EXPECT_TRUE(message.has_history());
  vessel_ = Vessel::ReadFromMessage(message, &ephemeris_, &parent_);
  EXPECT_TRUE(vessel_->is_initialized());
}

}  // namespace ksp_plugin
}  // namespace principia
