
#include "ksp_plugin/celestial.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using principia::si::Kilogram;
using principia::si::Metre;
using principia::si::Second;

namespace principia {
namespace ksp_plugin {

class CelestialTest : public testing::Test {
 protected:
  CelestialTest()
      : celestial_(
            make_not_null_unique<Celestial>(
                make_not_null_unique<MassiveBody>(1 * Kilogram))) {}

  not_null<std::unique_ptr<Celestial>> celestial_;
  DegreesOfFreedom<Barycentric> d1_ =
      {Barycentric::origin +
           Displacement<Barycentric>({1 * Metre, 2 * Metre, 3 * Metre}),
       Velocity<Barycentric>({4 * Metre / Second,
                              5 * Metre / Second,
                              6 * Metre / Second})};
  Instant const t1_ = kUniversalTimeEpoch + 42 * Second;
};

using CelestialDeathTest = CelestialTest;

TEST_F(CelestialDeathTest, Uninitialized) {
  EXPECT_DEATH({celestial_->history();}, "is_initialized");
  EXPECT_DEATH({celestial_->mutable_history();}, "is_initialized");
  EXPECT_DEATH({celestial_->prolongation();}, "is_initialized");
  EXPECT_DEATH({celestial_->mutable_prolongation();}, "is_initialized");
}

TEST_F(CelestialTest, Initialization) {
  EXPECT_FALSE(celestial_->is_initialized());
  celestial_->CreateHistoryAndForkProlongation(t1_, d1_);
  EXPECT_TRUE(celestial_->is_initialized());
}

TEST_F(CelestialDeathTest, SerializationError) {
  EXPECT_DEATH({
    serialization::Celestial message;
    celestial_->WriteToMessage(&message);
  }, "is_initialized");
}

TEST_F(CelestialTest, SerializationSuccess) {
  serialization::Celestial message;
  celestial_->CreateHistoryAndForkProlongation(t1_, d1_);
  celestial_->WriteToMessage(&message);
  EXPECT_TRUE(message.has_history_and_prolongation());
  celestial_ = Celestial::ReadFromMessage(message);
  EXPECT_TRUE(celestial_->is_initialized());
  EXPECT_EQ(d1_, celestial_->prolongation().last().degrees_of_freedom());
  EXPECT_EQ(t1_, celestial_->prolongation().last().time());
}

}  // namespace ksp_plugin
}  // namespace principia
