#include "ksp_plugin/physics_bubble.hpp"

#include "gtest/gtest.h"
#include "ksp_plugin/celestial.hpp"
#include "physics/body.hpp"
#include "quantities/quantities.hpp"

using principia::quantities::SIUnit;
using principia::quantities::Time;

namespace principia {
namespace ksp_plugin {

class PhysicsBubbleTest : public testing::Test {
 protected:
  PhysicsBubbleTest()
    : celestial_(std::make_unique<MassiveBody>(10 * SIUnit<Mass>())),
      vessel1_(&celestial_),
      vessel2_(&celestial_),
      t1_(1 * SIUnit<Time>()),
      t2_(2 * SIUnit<Time>()) {}

  PhysicsBubble bubble_;
  PhysicsBubble::PlanetariumRotation rotation_;
  Celestial celestial_;
  Vessel vessel1_;
  Vessel vessel2_;
  Instant t1_;
  Instant t2_;
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
  EXPECT_DEATH({
    bubble_.DisplacementCorrection(rotation_, celestial_, Position<World>());
  }, "Empty bubble");
  EXPECT_DEATH({
    bubble_.VelocityCorrection(rotation_, celestial_);
  }, "Empty bubble");
}

TEST_F(PhysicsBubbleTest, EmptySuccess) {
  EXPECT_TRUE(bubble_.empty());
  EXPECT_EQ(0, bubble_.size());
  EXPECT_EQ(0, bubble_.number_of_vessels());
  EXPECT_FALSE(bubble_.contains(&vessel1_));
  // Check that the following doesn't fail.  It does mostly nothing.
  bubble_.Prepare(rotation_, t1_, t2_);
}

}  // namespace ksp_plugin
}  // namespace principia
