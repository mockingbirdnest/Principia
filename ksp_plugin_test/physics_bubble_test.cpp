#include "ksp_plugin/physics_bubble.hpp"

#include "gtest/gtest.h"
#include "ksp_plugin/celestial.hpp"
#include "physics/body.hpp"
#include "quantities/quantities.hpp"

namespace principia {
namespace ksp_plugin {

class PhysicsBubbleTest : public testing::Test {
 protected:
  PhysicsBubbleTest()
    : celestial_(std::make_unique<MassiveBody>(10 * SIUnit<Mass>())),
      vessel1_(&celestial_),
      vessel2_(&celestial_) {}

  PhysicsBubble bubble_;
  Celestial celestial_;
  Vessel vessel1_;
  Vessel vessel2_;
};

using PhysicsBubbleDeathTest = PhysicsBubbleTest;

TEST_F(PhysicsBubbleDeathTest, EmptyError) {
  EXPECT_DEATH({
    bubble_.vessels();
  }, "Empty bubble");
  EXPECT_DEATH({
    bubble_.displacements_from_centre_of_mass(&vessel1_);
  }, "Empty bubble");
  EXPECT_DEATH({
    bubble_.velocities_from_centre_of_mass(&vessel1_);
  }, "Empty bubble");
  EXPECT_DEATH({
    bubble_.centre_of_mass_trajectory();
  }, "Empty bubble");
  EXPECT_DEATH({
    bubble_.mutable_centre_of_mass_trajectory();
  }, "Empty bubble");
}

TEST_F(PhysicsBubbleTest, EmptySuccess) {
  EXPECT_TRUE(bubble_.empty());
  EXPECT_EQ(0, bubble_.size());
  EXPECT_EQ(0, bubble_.number_of_vessels());
  EXPECT_FALSE(bubble_.contains(&vessel1_));
}

}  // namespace ksp_plugin
}  // namespace principia
