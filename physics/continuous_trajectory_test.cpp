
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
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace physics {
namespace internal_continuous_trajectory {

using geometry::Displacement;
using geometry::Frame;
using geometry::Velocity;
using numerics::ЧебышёвSeries;
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::Cos;
using quantities::Length;
using quantities::Sin;
using quantities::Time;
using quantities::astronomy::JulianYear;
using quantities::si::Kilo;
using quantities::si::Metre;
using quantities::si::Milli;
using quantities::si::Radian;
using quantities::si::Second;
using testing_utilities::AbsoluteError;
using testing_utilities::AlmostEquals;
using testing_utilities::EqualsProto;

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
      std::function<Position<World>(Instant const)> const& position_function,
      std::function<Velocity<World>(Instant const)> const& velocity_function,
      Instant const& time) {
    for (int i = 0; i < number_of_steps; ++i) {
      // We use this way of computing the time (as opposed to consecutive
      // additions) because it results in a bit of jitter in the intervals,
      // which matters for continuity.
      Instant ti = time + (i + 1) * step;
      trajectory_->Append(ti,
                          DegreesOfFreedom<World>(position_function(ti),
                                                  velocity_function(ti)));
    }
  }

  void ComputeBestNewhallApproximation(
      std::deque<Displacement<World>> const& error_estimates) {
    delete error_estimates_;
    error_estimates_ = new std::deque<Displacement<World>>(error_estimates);

    Instant const t = t0_ + 1 * Second;
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
  Instant const t0_;
  std::unique_ptr<ContinuousTrajectory<World>> trajectory_;
};

std::deque<Displacement<ContinuousTrajectoryTest::World>>*
ContinuousTrajectoryTest::error_estimates_ = nullptr;

TEST_F(ContinuousTrajectoryTest, BestNewhallApproximation) {
  Time const step = 1 * Second;
  Length const tolerance = 1 * Metre;

  trajectory_ = std::make_unique<ContinuousTrajectory<World>>(
                    step,
                    tolerance);
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
  EXPECT_EQ(tolerance, adjusted_tolerance());
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
  EXPECT_EQ(tolerance, adjusted_tolerance());
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
  EXPECT_EQ(tolerance, adjusted_tolerance());
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
  EXPECT_EQ(tolerance, adjusted_tolerance());
  EXPECT_FALSE(is_unstable());

  // Then we get low errors for a long time.
  for (int i = 0; i < 99; ++i) {
    ComputeBestNewhallApproximation(
        {Displacement<World>({0.1 * Metre, 0.1 * Metre, 0.1 * Metre})});
  }
  EXPECT_EQ(6, degree());
  EXPECT_EQ(tolerance, adjusted_tolerance());
  EXPECT_FALSE(is_unstable());

  // Finally we try all the degrees again and discover that degree 5 works.
  ComputeBestNewhallApproximation(
      {Displacement<World>({3 * Metre, 3 * Metre, 3 * Metre}),
       Displacement<World>({2 * Metre, 2 * Metre, 2 * Metre}),
       Displacement<World>({0.2 * Metre, 0.2 * Metre, 0.2 * Metre})});
  EXPECT_EQ(5, degree());
  EXPECT_EQ(tolerance, adjusted_tolerance());
  EXPECT_FALSE(is_unstable());
  ResetBestNewhallApproximation();
}

// A trajectory defined by a degree-1 polynomial.
TEST_F(ContinuousTrajectoryTest, Polynomial) {
  int const number_of_steps = 20;
  int const number_of_substeps = 50;
  Time const step = 0.01 * Second;

  auto position_function =
      [this](Instant const t) {
        return World::origin +
            Displacement<World>({(t - t0_) * 3 * Metre / Second,
                                 (t - t0_) * 5 * Metre / Second,
                                 (t - t0_) * (-2) * Metre / Second});
      };
  auto velocity_function =
      [](Instant const t) {
        return Velocity<World>({3 * Metre / Second,
                                5 * Metre / Second,
                                -2 * Metre / Second});
      };

  trajectory_ = std::make_unique<ContinuousTrajectory<World>>(
                    step,
                    /*tolerance=*/0.1 * Metre);

  EXPECT_TRUE(trajectory_->empty());
  FillTrajectory(
      number_of_steps, step, position_function, velocity_function, t0_);
  EXPECT_FALSE(trajectory_->empty());
  EXPECT_EQ(t0_ + step, trajectory_->t_min());
  EXPECT_EQ(t0_ + (((number_of_steps - 1) / 8) * 8 + 1) * step,
            trajectory_->t_max());

  for (Instant time = trajectory_->t_min();
       time <= trajectory_->t_max();
       time += step / number_of_substeps) {
    EXPECT_THAT(trajectory_->EvaluatePosition(time) - World::origin,
                AlmostEquals(position_function(time) - World::origin, 0, 11));
    EXPECT_THAT(trajectory_->EvaluateVelocity(time),
                AlmostEquals(velocity_function(time), 1, 3));
    EXPECT_EQ(trajectory_->EvaluateDegreesOfFreedom(time),
              DegreesOfFreedom<World>(trajectory_->EvaluatePosition(time),
                                      trajectory_->EvaluateVelocity(time)));
  }
}

// An approximation to the trajectory of Io.
TEST_F(ContinuousTrajectoryTest, Io) {
  int const number_of_steps = 200;
  int const number_of_substeps = 50;
  Length const sun_jupiter_distance = 778500000 * Kilo(Metre);
  Length const jupiter_io_distance = 421700 * Kilo(Metre);
  Time const jupiter_period = 11.8618 * JulianYear;
  Time const io_period = 152853.5047 * Second;
  Time const step = 3600 * Second;

  auto position_function =
      [this,
       sun_jupiter_distance,
       jupiter_io_distance,
       jupiter_period,
       io_period](Instant const t) {
        Angle const jupiter_angle = 2 * π * Radian * (t - t0_) / jupiter_period;
        Angle const io_angle = 2 * π * Radian * (t - t0_) / io_period;
        return World::origin +
            Displacement<World>({
                sun_jupiter_distance * Cos(jupiter_angle) +
                    jupiter_io_distance * Cos(io_angle),
                sun_jupiter_distance * Sin(jupiter_angle) +
                    jupiter_io_distance * Sin(io_angle),
                0 * Metre});
      };
  auto velocity_function = [this,
                            sun_jupiter_distance,
                            jupiter_io_distance,
                            jupiter_period,
                            io_period](Instant const t) {
    AngularFrequency const jupiter_ω = 2 * π * Radian / jupiter_period;
    AngularFrequency const io_ω = 2 * π * Radian / io_period;
    Angle const jupiter_angle = jupiter_ω *(t - t0_);
    Angle const io_angle = io_ω *(t - t0_);
    return Velocity<World>({
        (-jupiter_ω * sun_jupiter_distance * Sin(jupiter_angle) -
            io_ω * jupiter_io_distance * Sin(io_angle)) / Radian,
        (jupiter_ω * sun_jupiter_distance * Cos(jupiter_angle) +
            io_ω * jupiter_io_distance * Cos(io_angle)) / Radian,
        0 * Metre / Second});
  };

  trajectory_ = std::make_unique<ContinuousTrajectory<World>>(
                    step,
                    /*tolerance=*/5 * Milli(Metre));

  EXPECT_TRUE(trajectory_->empty());
  FillTrajectory(
      number_of_steps, step, position_function, velocity_function, t0_);
  EXPECT_FALSE(trajectory_->empty());
  EXPECT_EQ(t0_ + step, trajectory_->t_min());
  EXPECT_EQ(t0_ + (((number_of_steps - 1) / 8) * 8 + 1) * step,
            trajectory_->t_max());

  for (Instant time = trajectory_->t_min();
       time <= trajectory_->t_max();
       time += step / number_of_substeps) {
    Position<World> const actual_position = trajectory_->EvaluatePosition(time);
    Position<World> const expected_position = position_function(time);
    Velocity<World> const actual_velocity =
        trajectory_->EvaluateVelocity(time);
    Velocity<World> const expected_velocity = velocity_function(time);
    EXPECT_GT(0.491 * Milli(Metre),
              AbsoluteError(expected_position, actual_position));
    EXPECT_GT(1.60e-7 * Metre / Second,
              AbsoluteError(expected_velocity, actual_velocity));
  }

  trajectory_->ForgetBefore(trajectory_->t_min() - step);

  Instant const forget_before_time = t0_ + 44444 * Second;
  trajectory_->ForgetBefore(forget_before_time);
  EXPECT_EQ(forget_before_time, trajectory_->t_min());
  EXPECT_EQ(t0_ + (((number_of_steps - 1) / 8) * 8 + 1) * step,
            trajectory_->t_max());
  for (Instant time = trajectory_->t_min();
       time <= trajectory_->t_max();
       time += step / number_of_substeps) {
    Position<World> const actual_position = trajectory_->EvaluatePosition(time);
    Position<World> const expected_position = position_function(time);
    Velocity<World> const actual_velocity = trajectory_->EvaluateVelocity(time);
    Velocity<World> const expected_velocity = velocity_function(time);
    EXPECT_GT(0.492 * Milli(Metre),
              AbsoluteError(expected_position, actual_position));
    EXPECT_GT(1.60e-7 * Metre / Second,
              AbsoluteError(expected_velocity, actual_velocity));
  }
}

TEST_F(ContinuousTrajectoryTest, Continuity) {
  int const number_of_steps = 100;
  Length const distance = 1 * Kilo(Metre);
  Time const initial_time = 1e9 * Second;
  Time const period = 100 * Second;
  Time const step = 1 * Milli(Second);

  auto position_function = [this, distance, period](Instant const t) {
    Angle const angle = 2 * π * Radian * (t - t0_) / period;
    return World::origin +
        Displacement<World>({
            distance * Cos(angle),
            distance * Sin(angle),
            0 * Metre});
  };
  auto velocity_function = [this, distance, period](Instant const t) {
    AngularFrequency const ω = 2 * π * Radian / period;
    Angle const angle = ω * (t - t0_);
    return Velocity<World>({
        -ω * distance * Sin(angle) / Radian,
        ω * distance * Cos(angle) / Radian,
        0 * Metre / Second});
  };

  trajectory_ = std::make_unique<ContinuousTrajectory<World>>(
                    step,
                    /*tolerance=*/1 * Milli(Metre));

  EXPECT_TRUE(trajectory_->empty());
  FillTrajectory(number_of_steps,
                 step,
                 position_function,
                 velocity_function,
                 t0_ + initial_time);
  EXPECT_FALSE(trajectory_->empty());

  // This time is exactly at the continuity point of two consecutive series.
  int const interval = 11;
  Instant const continuity_time =
      t0_ + initial_time + (8 * interval + 1) * step;

  Position<World> const p1 =
      trajectory_->EvaluatePosition(continuity_time);
  Position<World> const p2 =
      trajectory_->EvaluatePosition(continuity_time + step);
  Position<World> const p3 =
      trajectory_->EvaluatePosition(continuity_time);
  EXPECT_THAT(p1, AlmostEquals(p3, 0, 2));
}

TEST_F(ContinuousTrajectoryTest, Serialization) {
  int const number_of_steps = 20;
  int const number_of_substeps = 50;
  Time const step = 0.01 * Second;
  Length const tolerance = 0.1 * Metre;

  auto position_function =
      [this](Instant const t) {
        return World::origin +
            Displacement<World>({(t - t0_) * 3 * Metre / Second,
                                 (t - t0_) * 5 * Metre / Second,
                                 (t - t0_) * (-2) * Metre / Second});
      };
  auto velocity_function =
      [](Instant const t) {
        return Velocity<World>({3 * Metre / Second,
                                5 * Metre / Second,
                                -2 * Metre / Second});
      };

  trajectory_ = std::make_unique<ContinuousTrajectory<World>>(
                    step, tolerance);

  EXPECT_TRUE(trajectory_->empty());
  FillTrajectory(
      number_of_steps, step, position_function, velocity_function, t0_);
  serialization::ContinuousTrajectory message;
  trajectory_->WriteToMessage(&message);
  EXPECT_EQ(step / Second, message.step().magnitude());
  EXPECT_EQ(tolerance / Metre, message.tolerance().magnitude());
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
       time += step / number_of_substeps) {
    EXPECT_EQ(trajectory->EvaluateDegreesOfFreedom(time),
              trajectory_->EvaluateDegreesOfFreedom(time));
  }

  serialization::ContinuousTrajectory second_message;
  trajectory->WriteToMessage(&second_message);
  EXPECT_THAT(message, EqualsProto(second_message));
}

TEST_F(ContinuousTrajectoryTest, Checkpoint) {
  int const number_of_steps1 = 30;
  int const number_of_steps2 = 20;
  int const number_of_substeps = 50;
  Time const step = 0.01 * Second;
  Length const tolerance = 0.1 * Metre;

  auto position_function =
      [this](Instant const t) {
        return World::origin +
            Displacement<World>({(t - t0_) * 3 * Metre / Second,
                                 (t - t0_) * 5 * Metre / Second,
                                 (t - t0_) * (-2) * Metre / Second});
      };
  auto velocity_function =
      [](Instant const t) {
        return Velocity<World>({3 * Metre / Second,
                                5 * Metre / Second,
                                -2 * Metre / Second});
      };

  trajectory_ = std::make_unique<ContinuousTrajectory<World>>(
                    step, tolerance);

  EXPECT_TRUE(trajectory_->empty());

  // Fill the trajectory, get a checkpoint and fill some more.
  FillTrajectory(
      number_of_steps1, step, position_function, velocity_function, t0_);
  ContinuousTrajectory<World>::Checkpoint const checkpoint =
      trajectory_->GetCheckpoint();
  Instant const checkpoint_t_max = trajectory_->t_max();
  FillTrajectory(number_of_steps2,
                 step,
                 position_function,
                 velocity_function,
                 t0_ + number_of_steps1 * step);

  serialization::ContinuousTrajectory message;
  trajectory_->WriteToMessage(&message, checkpoint);
  EXPECT_EQ(step / Second, message.step().magnitude());
  EXPECT_EQ(tolerance / Metre, message.tolerance().magnitude());
  EXPECT_GE(message.adjusted_tolerance().magnitude(),
            message.tolerance().magnitude());
  EXPECT_TRUE(message.has_is_unstable());
  EXPECT_EQ(3, message.degree());
  EXPECT_GE(100, message.degree_age());
  EXPECT_EQ(3, message.series_size());
  EXPECT_TRUE(message.has_first_time());
  EXPECT_EQ(6, message.last_point_size());

  auto const trajectory = ContinuousTrajectory<World>::ReadFromMessage(message);
  EXPECT_EQ(trajectory->t_min(), trajectory_->t_min());
  EXPECT_EQ(trajectory->t_max(), checkpoint_t_max);
  for (Instant time = trajectory_->t_min();
       time <= checkpoint_t_max;
       time += step / number_of_substeps) {
    EXPECT_EQ(trajectory->EvaluateDegreesOfFreedom(time),
              trajectory_->EvaluateDegreesOfFreedom(time));
  }
}

}  // namespace internal_continuous_trajectory
}  // namespace physics
}  // namespace principia
