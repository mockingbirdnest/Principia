#include "ksp_plugin/physics_bubble.hpp"

#include "gtest/gtest.h"
#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/part.hpp"
#include "ksp_plugin/vessel.hpp"
#include "physics/body.hpp"
#include "quantities/quantities.hpp"

using principia::quantities::Acceleration;
using principia::quantities::Speed;
using principia::quantities::SIUnit;
using principia::quantities::Time;

namespace principia {
namespace ksp_plugin {

class PhysicsBubbleTest : public testing::Test {
 protected:
  PhysicsBubbleTest()
    : p1a_(std::make_unique<Part<World>>(
               DegreesOfFreedom<World>(
                   Position<World>(Displacement<World>(
                       {14 * SIUnit<Length>(),
                        15 * SIUnit<Length>(),
                        16 * SIUnit<Length>()})),
                   Velocity<World>({114 * SIUnit<Speed>(),
                                    115 * SIUnit<Speed>(),
                                    116 * SIUnit<Speed>()})),
               11 * SIUnit<Mass>(),
               Vector<Acceleration, World>(
                   {114 * SIUnit<Acceleration>(),
                    115 * SIUnit<Acceleration>(),
                    116 * SIUnit<Acceleration>()}))),
      p1b_(std::make_unique<Part<World>>(
               DegreesOfFreedom<World>(
                   Position<World>(Displacement<World>(
                       {14 * SIUnit<Length>(),
                        15 * SIUnit<Length>(),
                        16 * SIUnit<Length>()})),
                   Velocity<World>({124 * SIUnit<Speed>(),
                                    125 * SIUnit<Speed>(),
                                    126 * SIUnit<Speed>()})),
               12 * SIUnit<Mass>(),
               Vector<Acceleration, World>(
                   {124 * SIUnit<Acceleration>(),
                    125 * SIUnit<Acceleration>(),
                    126 * SIUnit<Acceleration>()}))),
      p2a_(std::make_unique<Part<World>>(
               DegreesOfFreedom<World>(
                   Position<World>(Displacement<World>(
                       {24 * SIUnit<Length>(),
                        25 * SIUnit<Length>(),
                        26 * SIUnit<Length>()})),
                   Velocity<World>({214 * SIUnit<Speed>(),
                                    215 * SIUnit<Speed>(),
                                    216 * SIUnit<Speed>()})),
               21 * SIUnit<Mass>(),
               Vector<Acceleration, World>(
                   {214 * SIUnit<Acceleration>(),
                    215 * SIUnit<Acceleration>(),
                    216 * SIUnit<Acceleration>()}))),
      p2b_(std::make_unique<Part<World>>(
               DegreesOfFreedom<World>(
                   Position<World>(Displacement<World>(
                       {24 * SIUnit<Length>(),
                        25 * SIUnit<Length>(),
                        26 * SIUnit<Length>()})),
                   Velocity<World>({224 * SIUnit<Speed>(),
                                    225 * SIUnit<Speed>(),
                                    226 * SIUnit<Speed>()})),
                22 * SIUnit<Mass>(),
                Vector<Acceleration, World>(
                    {224 * SIUnit<Acceleration>(),
                     225 * SIUnit<Acceleration>(),
                     226 * SIUnit<Acceleration>()}))),
      p2c_(std::make_unique<Part<World>>(
               DegreesOfFreedom<World>(
                   Position<World>(Displacement<World>(
                       {24 * SIUnit<Length>(),
                        25 * SIUnit<Length>(),
                        26 * SIUnit<Length>()})),
                   Velocity<World>({234 * SIUnit<Speed>(),
                                    235 * SIUnit<Speed>(),
                                    236 * SIUnit<Speed>()})),
               23 * SIUnit<Mass>(),
               Vector<Acceleration, World>(
                   {234 * SIUnit<Acceleration>(),
                    235 * SIUnit<Acceleration>(),
                    236 * SIUnit<Acceleration>()}))),
      celestial_(std::make_unique<MassiveBody>(100 * SIUnit<Mass>())),
      vessel1_(&celestial_),
      vessel2_(&celestial_),
      t1_(1 * SIUnit<Time>()),
      t2_(2 * SIUnit<Time>()) {}

  PhysicsBubble bubble_;
  PhysicsBubble::PlanetariumRotation rotation_;
  std::unique_ptr<Part<World>> p1a_;
  std::unique_ptr<Part<World>> p1b_;
  std::unique_ptr<Part<World>> p2a_;
  std::unique_ptr<Part<World>> p2b_;
  std::unique_ptr<Part<World>> p2c_;
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

TEST_F(PhysicsBubbleTest, OneVessel) {
  std::vector<IdAndOwnedPart> parts;
  parts.push_back({11, std::move(p1a_)});
  parts.push_back({12, std::move(p1b_)});
  bubble_.AddVesselToNext(&vessel1_, std::move(parts));
  EXPECT_TRUE(bubble_.empty());
  bubble_.Prepare(rotation_, t1_, t2_);
}

}  // namespace ksp_plugin
}  // namespace principia
