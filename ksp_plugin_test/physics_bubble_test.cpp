#include "ksp_plugin/physics_bubble.hpp"

#include <map>
#include <string>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/grassmann.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ksp_plugin/celestial.hpp"
#include "ksp_plugin/frames.hpp"
#include "ksp_plugin/part.hpp"
#include "ksp_plugin/vessel.hpp"
#include "physics/body.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/componentwise.hpp"
#include "testing_utilities/vanishes_before.hpp"

namespace principia {

using base::make_not_null_unique;
using geometry::Bivector;
using geometry::Rotation;
using quantities::Acceleration;
using quantities::Speed;
using quantities::SIUnit;
using quantities::Time;
using quantities::si::Degree;
using quantities::si::Second;
using testing_utilities::AlmostEquals;
using testing_utilities::Componentwise;
using testing_utilities::VanishesBefore;
using ::testing::ElementsAre;

namespace ksp_plugin {

// All the numerical values in the tests below were computed exactly based on
// the formulae in the code and the data set up at construction.
class PhysicsBubbleTest : public testing::Test {
 protected:
  PhysicsBubbleTest()
    : rotation_(
          Rotation<Barycentric, WorldSun>(
              90 * Degree,
              Bivector<double, Barycentric>({0, 0, 1})).Forget()),
      celestial_dof_(Position<Barycentric>(Displacement<Barycentric>(
                         {-4 * SIUnit<Length>(),
                          -5 * SIUnit<Length>(),
                          -6 * SIUnit<Length>()})),
                     Velocity<Barycentric>({-4 * SIUnit<Speed>(),
                                            -5 * SIUnit<Speed>(),
                                            -6 * SIUnit<Speed>()})),
      dof1_(Position<Barycentric>(Displacement<Barycentric>(
                {104 * SIUnit<Length>(),
                 105 * SIUnit<Length>(),
                 106 * SIUnit<Length>()})),
            Velocity<Barycentric>({104 * SIUnit<Speed>(),
                                   105 * SIUnit<Speed>(),
                                   106 * SIUnit<Speed>()})),
      dof2_(Position<Barycentric>(Displacement<Barycentric>(
                {204 * SIUnit<Length>(),
                 205 * SIUnit<Length>(),
                 206 * SIUnit<Length>()})),
            Velocity<Barycentric>({204 * SIUnit<Speed>(),
                                   205 * SIUnit<Speed>(),
                                   206 * SIUnit<Speed>()})),
      body_(100 * SIUnit<Mass>()),
      celestial_(&body_),
      celestial_trajectory_(1 * Second, 1 * Metre),
      celestial_world_position_(Position<World>(Displacement<World>(
                                    {99 * SIUnit<Length>(),
                                     98 * SIUnit<Length>(),
                                     97 * SIUnit<Length>()}))),
      vessel1_(&celestial_),
      vessel2_(&celestial_),
      t1_(1 * SIUnit<Time>()),
      t2_(1.5 * SIUnit<Time>()),
      t3_(2 * SIUnit<Time>()) {
    celestial_.set_trajectory(&celestial_trajectory_);
    for (int i = 0; i < 9; ++i) {
      celestial_trajectory_.Append(t1_ + i * Second,
                                   {celestial_dof_.position() +
                                        i * Second * celestial_dof_.velocity(),
                                    celestial_dof_.velocity()});
    }
    vessel1_.CreateProlongation(t1_, dof1_);
    vessel2_.CreateProlongation(t1_, dof2_);

#if 0
    // Useful for debugging.
    google::SetStderrLogging(google::INFO);
    FLAGS_v = 1;
    FLAGS_logbuflevel = google::INFO - 1;
#endif
  }

  void CreateParts() {
    p1a_ = std::make_unique<Part<World>>(
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
                    116 * SIUnit<Acceleration>()}));
    p1b_ = std::make_unique<Part<World>>(
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
                    126 * SIUnit<Acceleration>()}));
    p2a_ = std::make_unique<Part<World>>(
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
                    216 * SIUnit<Acceleration>()}));
    p2b_ = std::make_unique<Part<World>>(
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
                    226 * SIUnit<Acceleration>()}));
    p2c_ = std::make_unique<Part<World>>(
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
                    236 * SIUnit<Acceleration>()}));
  }

  void CheckOneVesselDegreesOfFreedom(PhysicsBubble const& bubble) {
    // Since we have only one vessel, it is at the centre of mass of the bubble.
    EXPECT_THAT(bubble.from_centre_of_mass(&vessel1_),
                Componentwise(
                    Componentwise(VanishesBefore(1 * SIUnit<Length>(), 0),
                                  VanishesBefore(1 * SIUnit<Length>(), 0),
                                  VanishesBefore(1 * SIUnit<Length>(), 0)),
                    Componentwise(VanishesBefore(1 * SIUnit<Speed>(), 0),
                                  VanishesBefore(1 * SIUnit<Speed>(), 0),
                                  VanishesBefore(1 * SIUnit<Speed>(), 0))));

    // The trajectory of the centre of mass has one point which matches that
    // of the vessel.
    DiscreteTrajectory<Barycentric> const& trajectory =
        bubble.centre_of_mass_trajectory();
    EXPECT_EQ(dof1_, trajectory.last().degrees_of_freedom());
    DiscreteTrajectory<Barycentric>* mutable_trajectory =
        bubble.mutable_centre_of_mass_trajectory();
    EXPECT_EQ(dof1_, mutable_trajectory->last().degrees_of_freedom());

    if (mutable_trajectory->last().time() == t1_) {
      EXPECT_THAT(bubble_.DisplacementCorrection(
                      rotation_, celestial_, celestial_world_position_),
                  AlmostEquals(
                      Displacement<World>({-25 * SIUnit<Length>(),
                                           191 * SIUnit<Length>(),
                                           193 * SIUnit<Length>()}),
                      8));;
    } else if (mutable_trajectory->last().time() == t2_) {
      EXPECT_THAT(bubble_.DisplacementCorrection(
                      rotation_, celestial_, celestial_world_position_),
                  AlmostEquals(
                      Displacement<World>({-27.5 * SIUnit<Length>(),
                                           193 * SIUnit<Length>(),
                                           196 * SIUnit<Length>()}),
                      8));
    } else {
      ADD_FAILURE() << "unexpected last time";
    }
    EXPECT_THAT(bubble.VelocityCorrection(rotation_, celestial_),
                AlmostEquals(Velocity<World>(
                    {(-110.0 - 2742.0 / 23.0) * SIUnit<Speed>(),
                     (108.0 - 2765.0 / 23.0) * SIUnit<Speed>(),
                     (112.0 - 2788.0 / 23.0) * SIUnit<Speed>()}), 16));
  }

  void CheckTwoVesselsDegreesOfFreedom(PhysicsBubble const& bubble) {
    // Check the degrees of freedom of the vessels with respect to the centre of
    // mass.
    Displacement<World> const cdm_position({1906.0 / 89.0 * SIUnit<Length>(),
                                            1995.0 / 89.0 * SIUnit<Length>(),
                                            2084.0 / 89.0 * SIUnit<Length>()});
    Velocity<World> const cdm_velocity({17546.0 / 89.0 * SIUnit<Speed>(),
                                        17635.0 / 89.0 * SIUnit<Speed>(),
                                        17724.0 / 89.0 * SIUnit<Speed>()});
    EXPECT_THAT(bubble.from_centre_of_mass(&vessel1_),
                Componentwise(
                    AlmostEquals(Displacement<Barycentric>(
                        {15 * SIUnit<Length>() -
                            cdm_position.coordinates().y,
                         -14 * SIUnit<Length>() +
                            cdm_position.coordinates().x,
                         16 * SIUnit<Length>() -
                            cdm_position.coordinates().z}), 1),
                    AlmostEquals(Velocity<Barycentric>(
                        {2765.0 / 23.0 * SIUnit<Speed>() -
                            cdm_velocity.coordinates().y,
                         -2742.0 / 23.0 * SIUnit<Speed>() +
                            cdm_velocity.coordinates().x,
                         2788.0 / 23.0 * SIUnit<Speed>() -
                            cdm_velocity.coordinates().z}), 1)));
    EXPECT_THAT(bubble.from_centre_of_mass(&vessel2_),
                Componentwise(
                    AlmostEquals(Displacement<Barycentric>(
                        {25 * SIUnit<Length>() -
                            cdm_position.coordinates().y,
                         -24 * SIUnit<Length>() +
                            cdm_position.coordinates().x,
                         26 * SIUnit<Length>() -
                            cdm_position.coordinates().z}), 2),
                    AlmostEquals(Velocity<Barycentric>(
                        {7435.0 / 33.0 * SIUnit<Speed>() -
                            cdm_velocity.coordinates().y,
                         -7402.0 / 33.0 * SIUnit<Speed>() +
                            cdm_velocity.coordinates().x,
                         7468.0 / 33.0 * SIUnit<Speed>() -
                            cdm_velocity.coordinates().z}), 2)));

    // The trajectory of the centre of mass has only one point which is at the
    // barycentre of the trajectories of the vessels.
    DiscreteTrajectory<Barycentric> const& trajectory =
        bubble.centre_of_mass_trajectory();
    DegreesOfFreedom<Barycentric> const expected_dof =
        physics::Barycentre<Barycentric, double>({dof1_, dof2_}, {23, 66});
    EXPECT_EQ(expected_dof, trajectory.last().degrees_of_freedom());

    EXPECT_THAT(bubble.DisplacementCorrection(
                    rotation_, celestial_, celestial_world_position_),
                AlmostEquals(Displacement<World>(
                    {-7579.0 / 89.0 * SIUnit<Length>(),
                     24934.0 / 89.0 * SIUnit<Length>(),
                     25201 / 89.0  * SIUnit<Length>()}) - cdm_position, 5));
    EXPECT_THAT(bubble.VelocityCorrection(rotation_, celestial_),
                AlmostEquals(Velocity<World>(
                    {-16390.0 / 89.0 * SIUnit<Speed>(),
                     16212.0 / 89.0 * SIUnit<Speed>(),
                     16568.0 / 89.0 * SIUnit<Speed>()}) - cdm_velocity, 32));
  }

  PhysicsBubble bubble_;
  PhysicsBubble::BarycentricToWorldSun rotation_;
  DegreesOfFreedom<Barycentric> celestial_dof_;
  DegreesOfFreedom<Barycentric> dof1_;
  DegreesOfFreedom<Barycentric> dof2_;
  std::unique_ptr<Part<World>> p1a_;
  std::unique_ptr<Part<World>> p1b_;
  std::unique_ptr<Part<World>> p2a_;
  std::unique_ptr<Part<World>> p2b_;
  std::unique_ptr<Part<World>> p2c_;
  MassiveBody body_;
  Celestial celestial_;
  ContinuousTrajectory<Barycentric> celestial_trajectory_;
  Position<World> celestial_world_position_;
  Vessel vessel1_;
  Vessel vessel2_;
  Instant t1_;
  Instant t2_;
  Instant t3_;
};

using PhysicsBubbleDeathTest = PhysicsBubbleTest;

TEST_F(PhysicsBubbleDeathTest, EmptyError) {
  EXPECT_DEATH({
    bubble_.vessels();
  }, "Empty bubble");
  EXPECT_DEATH({
    bubble_.from_centre_of_mass(&vessel1_);
  }, "Empty bubble");
  EXPECT_DEATH({
    bubble_.centre_of_mass_trajectory();
  }, "Empty bubble");
  EXPECT_DEATH({
    bubble_.mutable_centre_of_mass_trajectory();
  }, "Empty bubble");
  EXPECT_DEATH({
    bubble_.DisplacementCorrection(rotation_,
                                   celestial_,
                                   Position<World>());
  }, "Empty bubble");
  EXPECT_DEATH({
    bubble_.VelocityCorrection(rotation_, celestial_);
  }, "Empty bubble");
}

TEST_F(PhysicsBubbleTest, EmptySuccess) {
  EXPECT_TRUE(bubble_.empty());
  EXPECT_EQ(0, bubble_.count());
  EXPECT_EQ(0, bubble_.number_of_vessels());
  EXPECT_FALSE(bubble_.contains(&vessel1_));
  // Check that the following doesn't fail.  It does mostly nothing.
  bubble_.Prepare(rotation_, t1_, t2_);
}

TEST_F(PhysicsBubbleTest, OneVesselOneStep) {
  std::vector<IdAndOwnedPart> parts;
  CreateParts();
  parts.emplace_back(11, std::move(p1a_));
  parts.emplace_back(12, std::move(p1b_));
  bubble_.AddVesselToNext(&vessel1_, std::move(parts));
  EXPECT_TRUE(bubble_.empty());

  bubble_.Prepare(rotation_, t1_, t2_);
  EXPECT_FALSE(bubble_.empty());
  EXPECT_EQ(1, bubble_.number_of_vessels());
  EXPECT_TRUE(bubble_.contains(&vessel1_));
  EXPECT_THAT(bubble_.vessels(), ElementsAre(&vessel1_));

  // The trajectory of the centre of mass has only one point and no
  // acceleration.
  DiscreteTrajectory<Barycentric> const& trajectory =
      bubble_.centre_of_mass_trajectory();
  EXPECT_THAT(trajectory.Times(), ElementsAre(t1_));
  EXPECT_FALSE(trajectory.has_intrinsic_acceleration());
  DiscreteTrajectory<Barycentric>* mutable_trajectory =
      bubble_.mutable_centre_of_mass_trajectory();
  EXPECT_THAT(mutable_trajectory->Times(), ElementsAre(t1_));
  EXPECT_FALSE(mutable_trajectory->has_intrinsic_acceleration());

  // Check the positions and velocities.
  CheckOneVesselDegreesOfFreedom(bubble_);
}

TEST_F(PhysicsBubbleTest, OneVesselTwoSteps) {
  std::vector<IdAndOwnedPart> parts;
  CreateParts();
  parts.emplace_back(11, std::move(p1a_));
  parts.emplace_back(12, std::move(p1b_));
  bubble_.AddVesselToNext(&vessel1_, std::move(parts));
  EXPECT_TRUE(bubble_.empty());

  bubble_.Prepare(rotation_, t1_, t2_);
  EXPECT_FALSE(bubble_.empty());
  // Necessary for the second call to Prepare to compute the acceleration.
  bubble_.VelocityCorrection(rotation_, celestial_);

  CreateParts();
  parts.emplace_back(11, std::move(p1a_));
  parts.emplace_back(12, std::move(p1b_));
  bubble_.AddVesselToNext(&vessel1_, std::move(parts));

  bubble_.Prepare(rotation_, t2_, t3_);
  EXPECT_FALSE(bubble_.empty());
  EXPECT_EQ(1, bubble_.number_of_vessels());
  EXPECT_TRUE(bubble_.contains(&vessel1_));
  EXPECT_THAT(bubble_.vessels(), ElementsAre(&vessel1_));

  // The trajectory now has an intrinsic acceleration.
  DiscreteTrajectory<Barycentric> const& trajectory =
      bubble_.centre_of_mass_trajectory();
  EXPECT_THAT(trajectory.Times(), ElementsAre(t1_));
  EXPECT_TRUE(trajectory.has_intrinsic_acceleration());
  EXPECT_THAT(trajectory.evaluate_intrinsic_acceleration(t2_),
              AlmostEquals(Vector<Acceleration, Barycentric>(
                               {(-2203.0 / 23.0) * SIUnit<Acceleration>(),
                                (-7802.0 / 23.0) * SIUnit<Acceleration>(),
                                (-2364.0 / 23.0) * SIUnit<Acceleration>()}),
                           3));

  // All the other assertions remain as in |OneVesselOneStep| since we didn't
  // restart or shift the bubble.
  CheckOneVesselDegreesOfFreedom(bubble_);
}

TEST_F(PhysicsBubbleTest, OneVesselPartRemoved) {
  std::vector<IdAndOwnedPart> parts;
  CreateParts();
  parts.emplace_back(11, std::move(p1a_));
  parts.emplace_back(12, std::move(p1b_));
  bubble_.AddVesselToNext(&vessel1_, std::move(parts));
  EXPECT_TRUE(bubble_.empty());

  bubble_.Prepare(rotation_, t1_, t2_);
  EXPECT_FALSE(bubble_.empty());
  // Necessary for the second call to Prepare to compute the acceleration.
  bubble_.VelocityCorrection(rotation_, celestial_);

  // For the second step the vessel only has one part.
  CreateParts();
  parts.emplace_back(12, std::move(p1b_));
  bubble_.AddVesselToNext(&vessel1_, std::move(parts));
  EXPECT_FALSE(bubble_.empty());

  bubble_.Prepare(rotation_, t2_, t3_);
  EXPECT_FALSE(bubble_.empty());
  EXPECT_EQ(1, bubble_.number_of_vessels());
  EXPECT_TRUE(bubble_.contains(&vessel1_));
  EXPECT_THAT(bubble_.vessels(), ElementsAre(&vessel1_));

  // The trajectory was reset by Shift.  Its intrinsic acceleration comes only
  // from |p1b_|.  The velocity of the centre of mass was updated to reflect the
  // change of parts.
  DiscreteTrajectory<Barycentric> const& trajectory =
      bubble_.centre_of_mass_trajectory();
  EXPECT_THAT(trajectory.Times(), ElementsAre(t2_));
  EXPECT_EQ(dof1_.position(),
            trajectory.last().degrees_of_freedom().position());
  EXPECT_EQ(dof1_.velocity() +
            Velocity<Barycentric>({(125.0 - 2765.0 / 23.0) * SIUnit<Speed>(),
                                   (-124.0 + 2742.0 / 23.0) * SIUnit<Speed>(),
                                   (126.0 - 2788.0 / 23.0) * SIUnit<Speed>()}),
            trajectory.last().degrees_of_freedom().velocity());
  EXPECT_TRUE(trajectory.has_intrinsic_acceleration());
  EXPECT_THAT(trajectory.evaluate_intrinsic_acceleration(t2_),
              AlmostEquals(Vector<Acceleration, Barycentric>(
                               {(-2313.0 / 23.0) * SIUnit<Acceleration>(),
                                (-7692.0 / 23.0) * SIUnit<Acceleration>(),
                                (-2474.0 / 23.0) * SIUnit<Acceleration>()}),
                           2));

  // The corrections are as in |OneVesselOneStep| except that the velocity is a
  // bit more precise (probably because we do less computations).
  EXPECT_THAT(bubble_.DisplacementCorrection(
                  rotation_, celestial_, celestial_world_position_),
              AlmostEquals(Displacement<World>({-27.5 * SIUnit<Length>(),
                                                193 * SIUnit<Length>(),
                                                196 * SIUnit<Length>()}), 8));
  EXPECT_THAT(bubble_.VelocityCorrection(rotation_, celestial_),
              AlmostEquals(Velocity<World>(
                  {(-110.0 - 2742.0 / 23.0) * SIUnit<Speed>(),
                   (108.0 - 2765.0 / 23.0) * SIUnit<Speed>(),
                   (112.0 - 2788.0 / 23.0) * SIUnit<Speed>()}), 8));
}

TEST_F(PhysicsBubbleTest, OneVesselPartAdded) {
  std::vector<IdAndOwnedPart> parts;
  // For the first step the vessel has only one part.
  CreateParts();
  parts.emplace_back(12, std::move(p1b_));
  bubble_.AddVesselToNext(&vessel1_, std::move(parts));
  EXPECT_TRUE(bubble_.empty());

  bubble_.Prepare(rotation_, t1_, t2_);
  EXPECT_FALSE(bubble_.empty());
  // Necessary for the second call to Prepare to compute the acceleration.
  bubble_.VelocityCorrection(rotation_, celestial_);

  // For the second step the vessel has two parts.
  CreateParts();
  parts.emplace_back(11, std::move(p1a_));
  parts.emplace_back(12, std::move(p1b_));
  bubble_.AddVesselToNext(&vessel1_, std::move(parts));

  bubble_.Prepare(rotation_, t2_, t3_);
  EXPECT_FALSE(bubble_.empty());
  EXPECT_EQ(1, bubble_.number_of_vessels());
  EXPECT_TRUE(bubble_.contains(&vessel1_));
  EXPECT_THAT(bubble_.vessels(), ElementsAre(&vessel1_));

  // The trajectory was reset by Shift.  Its intrinsic acceleration comes only
  // from |p1b_|.  The velocity of the centre of mass was updated to reflect the
  // change of parts, the adjustment is opposite to that of
  // |OneVesselPartRemoved|.
  DiscreteTrajectory<Barycentric> const& trajectory =
      bubble_.centre_of_mass_trajectory();
  EXPECT_THAT(trajectory.Times(), ElementsAre(t2_));
  EXPECT_EQ(dof1_.position(),
            trajectory.last().degrees_of_freedom().position());
  EXPECT_EQ(dof1_.velocity() +
            Velocity<Barycentric>({(-125.0 + 2765.0 / 23.0) * SIUnit<Speed>(),
                                   (124.0 - 2742.0 / 23.0) * SIUnit<Speed>(),
                                   (-126.0 + 2788.0 / 23.0) * SIUnit<Speed>()}),
            trajectory.last().degrees_of_freedom().velocity());
  EXPECT_TRUE(trajectory.has_intrinsic_acceleration());
  EXPECT_THAT(trajectory.evaluate_intrinsic_acceleration(t2_),
              AlmostEquals(Vector<Acceleration, Barycentric>(
                               {-91 * SIUnit<Acceleration>(),
                                -344 * SIUnit<Acceleration>(),
                                -98 * SIUnit<Acceleration>()}), 1));

  EXPECT_THAT(bubble_.DisplacementCorrection(
                  rotation_, celestial_, celestial_world_position_),
              AlmostEquals(Displacement<World>({-27.5 * SIUnit<Length>(),
                                                193 * SIUnit<Length>(),
                                                196 * SIUnit<Length>()}), 8));
  EXPECT_THAT(bubble_.VelocityCorrection(rotation_, celestial_),
              AlmostEquals(Velocity<World>(
                  {-234 * SIUnit<Speed>(),
                   -17 * SIUnit<Speed>(),
                   -14 * SIUnit<Speed>()}), 8));
}

TEST_F(PhysicsBubbleTest, OneVesselNoCommonParts) {
  std::vector<IdAndOwnedPart> parts;
  CreateParts();
  parts.emplace_back(21, std::move(p2a_));
  bubble_.AddVesselToNext(&vessel1_, std::move(parts));
  EXPECT_TRUE(bubble_.empty());

  bubble_.Prepare(rotation_, t1_, t2_);
  EXPECT_FALSE(bubble_.empty());
  // Necessary for the second call to Prepare to compute the acceleration.
  bubble_.VelocityCorrection(rotation_, celestial_);

  CreateParts();
  parts.emplace_back(11, std::move(p1a_));
  parts.emplace_back(12, std::move(p1b_));
  bubble_.AddVesselToNext(&vessel1_, std::move(parts));

  bubble_.Prepare(rotation_, t2_, t3_);
  EXPECT_FALSE(bubble_.empty());
  EXPECT_EQ(1, bubble_.number_of_vessels());
  EXPECT_TRUE(bubble_.contains(&vessel1_));
  EXPECT_THAT(bubble_.vessels(), ElementsAre(&vessel1_));

  // The bubble was restarted.
  DiscreteTrajectory<Barycentric> const& trajectory =
      bubble_.centre_of_mass_trajectory();
  EXPECT_THAT(trajectory.Times(), ElementsAre(t2_));
  EXPECT_FALSE(trajectory.has_intrinsic_acceleration());

  // All the other assertions remain as in |OneVesselOneStep| since we restarted
  // the bubble with parts |p1a_| and |p1b_|.
  CheckOneVesselDegreesOfFreedom(bubble_);
}

TEST_F(PhysicsBubbleTest, TwoVessels) {
  std::vector<IdAndOwnedPart> parts;
  CreateParts();
  parts.emplace_back(11, std::move(p1a_));
  parts.emplace_back(12, std::move(p1b_));
  bubble_.AddVesselToNext(&vessel1_, std::move(parts));
  parts.emplace_back(21, std::move(p2a_));
  parts.emplace_back(22, std::move(p2b_));
  parts.emplace_back(23, std::move(p2c_));
  bubble_.AddVesselToNext(&vessel2_, std::move(parts));
  EXPECT_TRUE(bubble_.empty());

  bubble_.Prepare(rotation_, t1_, t2_);
  EXPECT_FALSE(bubble_.empty());
  EXPECT_EQ(2, bubble_.number_of_vessels());
  EXPECT_TRUE(bubble_.contains(&vessel1_));
  EXPECT_TRUE(bubble_.contains(&vessel2_));
  EXPECT_THAT(bubble_.vessels(), ElementsAre(&vessel1_, &vessel2_));


  // The trajectory of the centre of mass has only one point.
  DiscreteTrajectory<Barycentric> const& trajectory =
      bubble_.centre_of_mass_trajectory();
  EXPECT_THAT(trajectory.Times(), ElementsAre(t1_));
  EXPECT_FALSE(trajectory.has_intrinsic_acceleration());

  // Check the positions and velocities.
  CheckTwoVesselsDegreesOfFreedom(bubble_);

  // The second step.
  CreateParts();
  parts.emplace_back(11, std::move(p1a_));
  parts.emplace_back(12, std::move(p1b_));
  bubble_.AddVesselToNext(&vessel1_, std::move(parts));
  parts.emplace_back(21, std::move(p2a_));
  parts.emplace_back(22, std::move(p2b_));
  parts.emplace_back(23, std::move(p2c_));
  bubble_.AddVesselToNext(&vessel2_, std::move(parts));
  EXPECT_FALSE(bubble_.empty());

  bubble_.Prepare(rotation_, t2_, t3_);
  EXPECT_FALSE(bubble_.empty());
  EXPECT_EQ(2, bubble_.number_of_vessels());
  EXPECT_TRUE(bubble_.contains(&vessel1_));
  EXPECT_TRUE(bubble_.contains(&vessel2_));
  EXPECT_THAT(bubble_.vessels(), ElementsAre(&vessel1_, &vessel2_));

  // Not much has changed, except that the trajectory now has an acceleration.
  EXPECT_THAT(trajectory.Times(), ElementsAre(t1_));
  Vector<Acceleration, World> const acceleration_correction =
      bubble_.VelocityCorrection(rotation_, celestial_) / (t3_ - t2_);
  EXPECT_TRUE(trajectory.has_intrinsic_acceleration());
  EXPECT_THAT(trajectory.evaluate_intrinsic_acceleration(t2_),
            AlmostEquals(Vector<Acceleration, Barycentric>(
                              {-acceleration_correction.coordinates().y -
                                   17635.0 / 89.0 * SIUnit<Acceleration>(),
                               acceleration_correction.coordinates().x +
                                   17546.0 / 89.0 * SIUnit<Acceleration>(),
                               -acceleration_correction.coordinates().z -
                                   17724.0 / 89.0 * SIUnit<Acceleration>()}),
                         2));

  // All the rest is identical since we didn't restart the bubble.
  CheckTwoVesselsDegreesOfFreedom(bubble_);
}

TEST_F(PhysicsBubbleTest, Serialization) {
  // Build a bubble similar to OneVesselOneStep.
  std::vector<IdAndOwnedPart> parts;
  CreateParts();
  parts.emplace_back(11, std::move(p1a_));
  parts.emplace_back(12, std::move(p1b_));
  bubble_.AddVesselToNext(&vessel1_, std::move(parts));
  bubble_.Prepare(rotation_, t1_, t2_);
  DiscreteTrajectory<Barycentric> const& trajectory =
      bubble_.centre_of_mass_trajectory();
  DiscreteTrajectory<Barycentric>* mutable_trajectory =
      bubble_.mutable_centre_of_mass_trajectory();
  bubble_.DisplacementCorrection(rotation_,
                                 celestial_,
                                 celestial_world_position_);
  bubble_.VelocityCorrection(rotation_, celestial_);

  std::map<std::string, not_null<Vessel*>> guid_to_vessel = {{"v1", &vessel1_}};
  std::map<not_null<Vessel const*>, std::string> vessel_to_guid =
      {{&vessel1_, "v1"}};

  serialization::PhysicsBubble message;
  bubble_.WriteToMessage([&vessel_to_guid](not_null<Vessel const*> const v) {
                           return vessel_to_guid.find(v)->second;
                         },
                         &message);
  EXPECT_TRUE(message.has_body());
  EXPECT_TRUE(message.has_current());
  EXPECT_EQ(2, message.current().part_size());
  EXPECT_EQ(11, message.current().part(0).part_id());
  EXPECT_EQ(11, message.current().part(0).part().mass().magnitude());
  EXPECT_EQ(12, message.current().part(1).part_id());
  EXPECT_EQ(12, message.current().part(1).part().mass().magnitude());
  EXPECT_EQ(1, message.current().vessel_size());
  EXPECT_EQ("v1", message.current().vessel(0).guid());
  EXPECT_EQ(2, message.current().vessel(0).part_id_size());
  EXPECT_EQ(11, message.current().vessel(0).part_id(0));
  EXPECT_EQ(12, message.current().vessel(0).part_id(1));
  EXPECT_TRUE(message.current().has_centre_of_mass());
  EXPECT_TRUE(message.current().has_centre_of_mass_trajectory());
  EXPECT_EQ(1, message.current().from_centre_of_mass_size());
  EXPECT_EQ("v1", message.current().from_centre_of_mass(0).guid());
  EXPECT_EQ(0, message.current().from_centre_of_mass(0).degrees_of_freedom().
                   t1().multivector().vector().x().quantity().magnitude());
  EXPECT_EQ(0, message.current().from_centre_of_mass(0).degrees_of_freedom().
                   t1().multivector().vector().y().quantity().magnitude());
  EXPECT_EQ(0, message.current().from_centre_of_mass(0).degrees_of_freedom().
                   t1().multivector().vector().z().quantity().magnitude());
  EXPECT_EQ(0, message.current().from_centre_of_mass(0).degrees_of_freedom().
                   t2().multivector().vector().x().quantity().magnitude());
  EXPECT_EQ(0, message.current().from_centre_of_mass(0).degrees_of_freedom().
                   t2().multivector().vector().y().quantity().magnitude());
  EXPECT_EQ(0, message.current().from_centre_of_mass(0).degrees_of_freedom().
                   t2().multivector().vector().z().quantity().magnitude());
  EXPECT_TRUE(message.current().has_displacement_correction());
  EXPECT_TRUE(message.current().has_velocity_correction());

  auto const bubble = PhysicsBubble::ReadFromMessage(
                          [&guid_to_vessel](std::string const& v) {
                            return guid_to_vessel.find(v)->second;
                          },
                          message);
  CheckOneVesselDegreesOfFreedom(*bubble);
}


}  // namespace ksp_plugin
}  // namespace principia
