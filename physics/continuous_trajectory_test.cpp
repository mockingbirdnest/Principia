#include "physics/continuous_trajectory.hpp"

#include <functional>

#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {

using astronomy::JulianYear;
using geometry::Displacement;
using geometry::Frame;
using geometry::Velocity;
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::Time;
using si::Kilo;
using si::Metre;
using si::Milli;
using si::Radian;
using si::Second;
using testing_utilities::AbsoluteError;
using testing_utilities::AlmostEquals;

namespace physics {

class ContinuousTrajectoryTest : public testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST1, true>;

  void FillTrajectory(
      int const number_of_steps,
      Time const& step,
      std::function<Position<World>(Instant const)> position_function,
      std::function<Velocity<World>(Instant const)> velocity_function) {
    Instant time;
    for (int i = 0; i < number_of_steps; ++i) {
      time += step;
      trajectory_->Append(time,
                          DegreesOfFreedom<World>(position_function(time),
                                                  velocity_function(time)));
    }
  }

  std::unique_ptr<ContinuousTrajectory<World>> trajectory_;
};

TEST_F(ContinuousTrajectoryTest, Polynomial) {
  int const kNumberOfSteps = 20;
  int const kNumberOfSubsteps = 50;
  Time const kStep = 0.01 * Second;
  Instant const t0;

  auto position_function =
      [t0](Instant const t) {
        return World::origin +
            Displacement<World>({(t - t0) * 3 * Metre / Second,
                                (t - t0) * 5 * Metre / Second,
                                (t - t0) * (-2) * Metre / Second});
      };
  auto velocity_function =
      [t0](Instant const t) {
        return Velocity<World>({3 * Metre / Second,
                                5 * Metre / Second,
                                -2 * Metre / Second});
      };

  trajectory_ = std::make_unique<ContinuousTrajectory<World>>(
                    kStep,
                    0.05 * Metre /*low_tolerance*/,
                    0.1 * Metre /*high_tolerance*/);

  EXPECT_TRUE(trajectory_->empty());
  FillTrajectory(kNumberOfSteps, kStep, position_function, velocity_function);
  EXPECT_FALSE(trajectory_->empty());
  EXPECT_EQ(t0 + kStep, trajectory_->t_min());
  EXPECT_EQ(t0 + (((kNumberOfSteps - 1) / 8) * 8 + 1) * kStep,
            trajectory_->t_max());

  ContinuousTrajectory<World>::Hint hint;
  for (Instant time = trajectory_->t_min();
       time <= trajectory_->t_max();
       time += kStep / kNumberOfSubsteps) {
    EXPECT_THAT(trajectory_->EvaluatePosition(time, &hint) - World::origin,
                AlmostEquals(position_function(time) - World::origin, 0, 7));
    EXPECT_THAT(trajectory_->EvaluateVelocity(time, &hint),
                AlmostEquals(velocity_function(time), 1, 3));
    EXPECT_EQ(trajectory_->EvaluateDegreesOfFreedom(time, &hint),
              DegreesOfFreedom<World>(
                  trajectory_->EvaluatePosition(time, &hint),
                  trajectory_->EvaluateVelocity(time, &hint)));
  }
}

// An approximation to the trajectory of Io.
TEST_F(ContinuousTrajectoryTest, Io) {
  int const kNumberOfSteps = 200;
  int const kNumberOfSubsteps = 50;
  Length const kSunJupiterDistance = 778500000 * Kilo(Metre);
  Length const kJupiterIoDistance = 421700 * Kilo(Metre);
  Time const kJupiterPeriod = 11.8618 * JulianYear;
  Time const kIoPeriod = 152853.5047 * Second;
  Time const kStep = 3600 * Second;
  Instant const t0;

  auto position_function =
      [t0, kSunJupiterDistance, kJupiterIoDistance, kJupiterPeriod, kIoPeriod](
          Instant const t) {
        Angle const jupiter_angle = 2 * π * Radian * (t - t0) / kJupiterPeriod;
        Angle const io_angle = 2 * π * Radian * (t - t0) / kIoPeriod;
        return World::origin +
            Displacement<World>({
                kSunJupiterDistance * Cos(jupiter_angle) +
                    kJupiterIoDistance * Cos(io_angle),
                kSunJupiterDistance * Sin(jupiter_angle) +
                    kJupiterIoDistance * Sin(io_angle),
                0 * Metre});
      };
  auto velocity_function =
      [t0, kSunJupiterDistance, kJupiterIoDistance, kJupiterPeriod, kIoPeriod](
          Instant const t) {
        AngularFrequency const jupiter_ω = 2 * π * Radian / kJupiterPeriod;
        AngularFrequency const io_ω = 2 * π * Radian / kIoPeriod;
        Angle const jupiter_angle = jupiter_ω * (t - t0);
        Angle const io_angle = io_ω * (t - t0);
        return Velocity<World>({
            (-jupiter_ω * kSunJupiterDistance * Sin(jupiter_angle) -
                io_ω * kJupiterIoDistance * Sin(io_angle)) / Radian,
            (jupiter_ω * kSunJupiterDistance * Cos(jupiter_angle) +
                io_ω * kJupiterIoDistance * Cos(io_angle)) / Radian,
            0 * Metre / Second});
      };

  trajectory_ = std::make_unique<ContinuousTrajectory<World>>(
                    kStep,
                    1 * Milli(Metre) /*low_tolerance*/,
                    5 * Milli(Metre) /*high_tolerance*/);

  EXPECT_TRUE(trajectory_->empty());
  FillTrajectory(kNumberOfSteps, kStep, position_function, velocity_function);
  EXPECT_FALSE(trajectory_->empty());
  EXPECT_EQ(t0 + kStep, trajectory_->t_min());
  EXPECT_EQ(t0 + (((kNumberOfSteps - 1) / 8) * 8 + 1) * kStep,
            trajectory_->t_max());

  ContinuousTrajectory<World>::Hint hint;
  for (Instant time = trajectory_->t_min();
       time <= trajectory_->t_max();
       time += kStep / kNumberOfSubsteps) {
    Displacement<World> const actual_displacement =
        trajectory_->EvaluatePosition(time, &hint) - World::origin;
    Displacement<World> const expected_displacement =
        position_function(time) - World::origin;
    Velocity<World> const actual_velocity =
        trajectory_->EvaluateVelocity(time, &hint);
    Velocity<World> const expected_velocity = velocity_function(time);
    EXPECT_GT(0.492 * Milli(Metre),
              AbsoluteError(expected_displacement, actual_displacement));
    EXPECT_GT(1.60E-7 * Metre / Second,
              AbsoluteError(expected_velocity, actual_velocity));
  }

  Instant const kForgetTime = t0 + 44444 * Second;
  trajectory_->ForgetBefore(kForgetTime);
  EXPECT_EQ(kForgetTime, trajectory_->t_min());
  EXPECT_EQ(t0 + (((kNumberOfSteps - 1) / 8) * 8 + 1) * kStep,
            trajectory_->t_max());
  for (Instant time = trajectory_->t_min();
       time <= trajectory_->t_max();
       time += kStep / kNumberOfSubsteps) {
    Displacement<World> const actual_displacement =
        trajectory_->EvaluatePosition(time, &hint) - World::origin;
    Displacement<World> const expected_displacement =
        position_function(time) - World::origin;
    Velocity<World> const actual_velocity =
        trajectory_->EvaluateVelocity(time, &hint);
    Velocity<World> const expected_velocity = velocity_function(time);
    EXPECT_GT(0.492 * Milli(Metre),
              AbsoluteError(expected_displacement, actual_displacement));
    EXPECT_GT(1.60E-7 * Metre / Second,
              AbsoluteError(expected_velocity, actual_velocity));
  }
}

}  // namespace physics
}  // namespace principia
