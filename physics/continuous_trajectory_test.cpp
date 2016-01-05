#include "physics/continuous_trajectory.hpp"

#include <deque>
#include <functional>
#include <limits>
#include <vector>

#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "numerics/чебышёв_series.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "serialization/physics.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {

using geometry::Displacement;
using geometry::Frame;
using geometry::Velocity;
using numerics::ЧебышёвSeries;
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::Time;
using quantities::astronomy::JulianYear;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AbsoluteError;
using testing_utilities::AlmostEquals;

namespace physics {

class ContinuousTrajectoryTest : public testing::Test {
 public:
  using World = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST1, true>;

 protected:
  static ЧебышёвSeries<Displacement<World>> SimulatedNewhallApproximation(
      int const degree,
      std::vector<Displacement<World>> const& q,
      std::vector<Velocity<World>> const& v,
      Instant const& t_min,
      Instant const& t_max) {
    Displacement<World> const error_estimate = error_estimates_->front();
    error_estimates_->pop_front();
    return ЧебышёвSeries<Displacement<World>>({error_estimate}, t_min, t_max);
  }

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

  void ComputeBestNewhallApproximation(
      std::deque<Displacement<World>> const& error_estimates) {
    delete error_estimates_;
    error_estimates_ = new std::deque<Displacement<World>>(error_estimates);

    Instant const t0;
    Instant const t = t0 + 1 * Second;
    std::vector<Displacement<World>> const q;
    std::vector<Velocity<World>> const v;
    trajectory_->ComputeBestNewhallApproximation(
        t, q, v, &SimulatedNewhallApproximation);
  }

  int degree() const {
    return trajectory_->degree_;
  }

  Length adjusted_tolerance() const {
    return trajectory_->adjusted_tolerance_;
  }

  bool is_unstable() const {
    return trajectory_->is_unstable_;
  }

  void ResetBestNewhallApproximation() {
    trajectory_->degree_age_ = std::numeric_limits<int>::max();
  }

  static std::deque<Displacement<World>>* error_estimates_;
  std::unique_ptr<ContinuousTrajectory<World>> trajectory_;
};

std::deque<Displacement<ContinuousTrajectoryTest::World>>*
ContinuousTrajectoryTest::error_estimates_ = nullptr;

TEST_F(ContinuousTrajectoryTest, BestNewhallApproximation) {
  Time const kStep = 1 * Second;
  Length const kTolerance = 1 * Metre;

  trajectory_ = std::make_unique<ContinuousTrajectory<World>>(
                    kStep,
                    kTolerance);
  trajectory_->Append(Instant(),
                      DegreesOfFreedom<World>(Position<World>(),
                                              Velocity<World>()));

  // A case where the errors smoothly decrease.
  ComputeBestNewhallApproximation(
      {Displacement<World>({3 * Metre, 4 * Metre, 5 * Metre}),
       Displacement<World>({2 * Metre, 1 * Metre, 2 * Metre}),
       Displacement<World>({0.1 * Metre, 2 * Metre, 0 * Metre}),
       Displacement<World>({0.5 * Metre, 0.5 * Metre, 0.1 * Metre})});
  EXPECT_EQ(6, degree());
  EXPECT_EQ(kTolerance, adjusted_tolerance());
  EXPECT_FALSE(is_unstable());
  ResetBestNewhallApproximation();

  // A case where the errors increase at the end, but after we have reached the
  // desired tolerance.
  ComputeBestNewhallApproximation(
      {Displacement<World>({3 * Metre, 4 * Metre, 5 * Metre}),
       Displacement<World>({2 * Metre, 1 * Metre, 2 * Metre}),
       Displacement<World>({0.1 * Metre, 2 * Metre, 0 * Metre}),
       Displacement<World>({0.5 * Metre, 0.5 * Metre, 0.1 * Metre}),
       Displacement<World>({1 * Metre, 3 * Metre, 1 * Metre})});
  EXPECT_EQ(6, degree());
  EXPECT_EQ(kTolerance, adjusted_tolerance());
  EXPECT_FALSE(is_unstable());
  ResetBestNewhallApproximation();

  // A case where the errors increase before we have reach the desired
  // tolerance...
  ComputeBestNewhallApproximation(
      {Displacement<World>({3 * Metre, 4 * Metre, 5 * Metre}),
       Displacement<World>({2 * Metre, 1 * Metre, 2 * Metre}),
       Displacement<World>({0.1 * Metre, 2 * Metre, 0 * Metre}),
       Displacement<World>({1 * Metre, 3 * Metre, 1 * Metre}),
       Displacement<World>({0.5 * Metre, 0.5 * Metre, 0.1 * Metre})});
  EXPECT_EQ(5, degree());
  EXPECT_EQ(sqrt(4.01) * Metre, adjusted_tolerance());
  EXPECT_TRUE(is_unstable());

  // ... then the error decreases...
  ComputeBestNewhallApproximation(
      {Displacement<World>({0.1 * Metre, 1.5 * Metre, 0 * Metre})});
  EXPECT_EQ(5, degree());
  EXPECT_EQ(sqrt(4.01) * Metre, adjusted_tolerance());
  EXPECT_TRUE(is_unstable());

  // ... then the error increases forcing us to go back to square one...
  ComputeBestNewhallApproximation(
      {Displacement<World>({0.1 * Metre, 2 * Metre, 0.5 * Metre}),
       Displacement<World>({3 * Metre, 4 * Metre, 5 * Metre}),
       Displacement<World>({2 * Metre, 1 * Metre, 2 * Metre}),
       Displacement<World>({1 * Metre, 2 * Metre, 1 * Metre}),
       Displacement<World>({1 * Metre, 1.5 * Metre, 1 * Metre}),
       Displacement<World>({1 * Metre, 1.2 * Metre, 1 * Metre}),
       Displacement<World>({1 * Metre, 1.3 * Metre, 1 * Metre})});
  EXPECT_EQ(7, degree());
  EXPECT_EQ(sqrt(3.44) * Metre, adjusted_tolerance());
  EXPECT_TRUE(is_unstable());

  // ... it does it again but then the computation becomes stable.
  ComputeBestNewhallApproximation(
      {Displacement<World>({1 * Metre, 1.3 * Metre, 1 * Metre}),
       Displacement<World>({3 * Metre, 4 * Metre, 5 * Metre}),
       Displacement<World>({0.1 * Metre, 0.5 * Metre, 0.2 * Metre})});
  EXPECT_EQ(4, degree());
  EXPECT_EQ(kTolerance, adjusted_tolerance());
  EXPECT_FALSE(is_unstable());
  ResetBestNewhallApproximation();

  // Check that the degree is properly lowered when the age of the approximation
  // exceeds the limit.
  // First, the errors force usage of degree 6.
  ComputeBestNewhallApproximation(
      {Displacement<World>({3 * Metre, 3 * Metre, 3 * Metre}),
       Displacement<World>({2 * Metre, 2 * Metre, 2 * Metre}),
       Displacement<World>({1 * Metre, 1 * Metre, 1 * Metre}),
       Displacement<World>({0.1 * Metre, 0.1 * Metre, 0.1 * Metre})});
  EXPECT_EQ(6, degree());
  EXPECT_EQ(kTolerance, adjusted_tolerance());
  EXPECT_FALSE(is_unstable());

  // Then we get low errors for a long time.
  for (int i = 0; i < 99; ++i) {
    ComputeBestNewhallApproximation(
        {Displacement<World>({0.1 * Metre, 0.1 * Metre, 0.1 * Metre})});
  }
  EXPECT_EQ(6, degree());
  EXPECT_EQ(kTolerance, adjusted_tolerance());
  EXPECT_FALSE(is_unstable());

  // Finally we try all the degrees again and discover that degree 5 works.
  ComputeBestNewhallApproximation(
      {Displacement<World>({3 * Metre, 3 * Metre, 3 * Metre}),
       Displacement<World>({2 * Metre, 2 * Metre, 2 * Metre}),
       Displacement<World>({0.2 * Metre, 0.2 * Metre, 0.2 * Metre})});
  EXPECT_EQ(5, degree());
  EXPECT_EQ(kTolerance, adjusted_tolerance());
  EXPECT_FALSE(is_unstable());
  ResetBestNewhallApproximation();
}

// A trajectory defined by a degree-1 polynomial.
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
                    0.1 * Metre /*tolerance*/);

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
                    5 * Milli(Metre) /*tolerance*/);

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
    Position<World> const actual_position =
        trajectory_->EvaluatePosition(time, &hint);
    Position<World> const expected_position = position_function(time);
    Velocity<World> const actual_velocity =
        trajectory_->EvaluateVelocity(time, &hint);
    Velocity<World> const expected_velocity = velocity_function(time);
    EXPECT_GT(0.491 * Milli(Metre),
              AbsoluteError(expected_position, actual_position));
    EXPECT_GT(1.60E-7 * Metre / Second,
              AbsoluteError(expected_velocity, actual_velocity));
  }

  trajectory_->ForgetBefore(trajectory_->t_min() - kStep);

  Instant const kForgetBeforeTime = t0 + 44444 * Second;
  trajectory_->ForgetBefore(kForgetBeforeTime);
  EXPECT_EQ(kForgetBeforeTime, trajectory_->t_min());
  EXPECT_EQ(t0 + (((kNumberOfSteps - 1) / 8) * 8 + 1) * kStep,
            trajectory_->t_max());
  for (Instant time = trajectory_->t_min();
       time <= trajectory_->t_max();
       time += kStep / kNumberOfSubsteps) {
    Position<World> const actual_position =
        trajectory_->EvaluatePosition(time, &hint);
    Position<World> const expected_position = position_function(time);
    Velocity<World> const actual_velocity =
        trajectory_->EvaluateVelocity(time, &hint);
    Velocity<World> const expected_velocity = velocity_function(time);
    EXPECT_GT(0.492 * Milli(Metre),
              AbsoluteError(expected_position, actual_position));
    EXPECT_GT(1.60E-7 * Metre / Second,
              AbsoluteError(expected_velocity, actual_velocity));
  }

  Instant const kForgetAfterTime = t0 + 25 * kStep;
  trajectory_->ForgetAfter(
      kForgetAfterTime,
      trajectory_->EvaluateDegreesOfFreedom(kForgetAfterTime,
                                            nullptr /*hint*/));
  EXPECT_EQ(kForgetBeforeTime, trajectory_->t_min());
  EXPECT_EQ(kForgetAfterTime, trajectory_->t_max());
  for (Instant time = trajectory_->t_min();
       time <= trajectory_->t_max();
       time += kStep / kNumberOfSubsteps) {
    Position<World> const actual_position =
      trajectory_->EvaluatePosition(time, &hint);
    Position<World> const expected_position = position_function(time);
    Velocity<World> const actual_velocity =
      trajectory_->EvaluateVelocity(time, &hint);
    Velocity<World> const expected_velocity = velocity_function(time);
    EXPECT_GT(0.492 * Milli(Metre),
              AbsoluteError(expected_position, actual_position));
    EXPECT_GT(1.60E-7 * Metre / Second,
              AbsoluteError(expected_velocity, actual_velocity));
  }

  trajectory_->ForgetBefore(trajectory_->t_max() + kStep);
}

TEST_F(ContinuousTrajectoryTest, Serialization) {
  int const kNumberOfSteps = 20;
  int const kNumberOfSubsteps = 50;
  Time const kStep = 0.01 * Second;
  Length const kTolerance = 0.1 * Metre;
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
                    kStep, kTolerance);

  EXPECT_TRUE(trajectory_->empty());
  FillTrajectory(kNumberOfSteps, kStep, position_function, velocity_function);

  serialization::ContinuousTrajectory message;
  trajectory_->WriteToMessage(&message);
  EXPECT_EQ(kStep / Second, message.step().magnitude());
  EXPECT_EQ(kTolerance / Metre, message.tolerance().magnitude());
  EXPECT_GE(message.adjusted_tolerance().magnitude(),
            message.tolerance().magnitude());
  EXPECT_TRUE(message.has_is_unstable());
  EXPECT_EQ(3, message.degree());
  EXPECT_GE(100, message.degree_age());
  EXPECT_EQ(2, message.series_size());
  EXPECT_TRUE(message.has_first_time());
  EXPECT_EQ(4, message.last_point_size());

  auto const trajectory = ContinuousTrajectory<World>::ReadFromMessage(message);
  EXPECT_EQ(trajectory->t_min(), trajectory_->t_min());
  EXPECT_EQ(trajectory->t_max(), trajectory_->t_max());
  for (Instant time = trajectory_->t_min();
       time <= trajectory_->t_max();
       time += kStep / kNumberOfSubsteps) {
    EXPECT_EQ(trajectory->EvaluateDegreesOfFreedom(time, nullptr /*hint*/),
              trajectory_->EvaluateDegreesOfFreedom(time, nullptr /*hint*/));
  }
}

}  // namespace physics
}  // namespace principia
