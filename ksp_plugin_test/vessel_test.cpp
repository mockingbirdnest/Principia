
#include "ksp_plugin/vessel.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using principia::si::Kilogram;
using principia::si::Metre;
using principia::si::Second;
using testing::IsNull;
using testing::NotNull;

namespace principia {
namespace ksp_plugin {

class VesselTest : public testing::Test {
 protected:
  VesselTest()
      : parent_(make_not_null_unique<MassiveBody>(1 * Kilogram)),
        vessel_(make_not_null_unique<Vessel>(&parent_)) {}

  Celestial parent_;
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
  Instant const t2_ = kUniversalTimeEpoch + 42 * Second;
};

using VesselDeathTest = VesselTest;

TEST_F(VesselDeathTest, Uninitialized) {
  EXPECT_DEATH({vessel_->history();}, "is_synchronized");
  EXPECT_DEATH({vessel_->mutable_history();}, "is_synchronized");
  EXPECT_DEATH({vessel_->prolongation();}, "is_initialized");
  EXPECT_DEATH({vessel_->mutable_prolongation();}, "is_initialized");
}

TEST_F(VesselDeathTest, Unsynchronized) {
  EXPECT_DEATH({
    vessel_->CreateProlongation(t1_, d1_);
    vessel_->history();
  }, "is_synchronized");
  EXPECT_DEATH({
    vessel_->CreateProlongation(t1_, d1_);
    vessel_->mutable_history();
  }, "is_synchronized");
}

TEST_F(VesselTest, InitializationAndSynchronization) {
  // TODO(egg): mutable_history() should probably return a |not_null| and check
  // that |is_initialized()|.  This is even more confusing for the non-mutable
  // versions, which return a |const&|.
  EXPECT_FALSE(vessel_->is_initialized());
  EXPECT_FALSE(vessel_->is_synchronized());
  vessel_->CreateProlongation(t1_, d1_);
  EXPECT_TRUE(vessel_->is_initialized());
  EXPECT_FALSE(vessel_->is_synchronized());
  EXPECT_THAT(vessel_->mutable_prolongation(), NotNull());
  vessel_->CreateHistoryAndForkProlongation(t2_, d2_);
  EXPECT_TRUE(vessel_->is_initialized());
  EXPECT_TRUE(vessel_->is_synchronized());
  EXPECT_THAT(vessel_->mutable_history(), NotNull());
  EXPECT_THAT(vessel_->mutable_prolongation(), NotNull());
}

TEST_F(VesselDeathTest, SerializationError) {
  EXPECT_DEATH({
    serialization::Vessel message;
    vessel_->WriteToMessage(&message);
  }, "is_initialized");
  EXPECT_DEATH({
    serialization::Vessel message;
    Vessel::ReadFromMessage(message, &parent_);
  }, "message does not represent an initialized Vessel");
}

TEST_F(VesselTest, SerializationSuccess) {
  serialization::Vessel message;
  EXPECT_FALSE(message.has_owned_prolongation());
  EXPECT_FALSE(message.has_history_and_prolongation());
  vessel_->CreateProlongation(t1_, d1_);
  vessel_->WriteToMessage(&message);
  EXPECT_TRUE(message.has_owned_prolongation());
  EXPECT_FALSE(message.has_history_and_prolongation());
  vessel_ = Vessel::ReadFromMessage(message, &parent_);
  EXPECT_TRUE(vessel_->is_initialized());
  EXPECT_FALSE(vessel_->is_synchronized());
  vessel_->CreateHistoryAndForkProlongation(t2_, d2_);
  message.Clear();
  vessel_->WriteToMessage(&message);
  EXPECT_FALSE(message.has_owned_prolongation());
  EXPECT_TRUE(message.has_history_and_prolongation());
  vessel_ = Vessel::ReadFromMessage(message, &parent_);
  EXPECT_TRUE(vessel_->is_initialized());
  EXPECT_TRUE(vessel_->is_synchronized());
}

}  // namespace ksp_plugin
}  // namespace principia
