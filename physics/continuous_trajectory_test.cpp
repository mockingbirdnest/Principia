
#include "physics/continuous_trajectory.hpp"

#include <algorithm>
#include <deque>
#include <functional>
#include <limits>
#include <vector>

#include "geometry/frame.hpp"
#include "geometry/named_quantities.hpp"
#include "gtest/gtest.h"
#include "numerics/polynomial.hpp"
#include "numerics/polynomial_evaluators.hpp"
#include "numerics/чебышёв_series.hpp"
#include "physics/degrees_of_freedom.hpp"
#include "quantities/astronomy.hpp"
#include "quantities/numbers.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "serialization/geometry.pb.h"
#include "serialization/physics.pb.h"
#include "testing_utilities/almost_equals.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace physics {
namespace internal_continuous_trajectory {

using geometry::Displacement;
using geometry::Frame;
using geometry::Velocity;
using numerics::Polynomial;
using numerics::PolynomialInMonomialBasis;
using numerics::HornerEvaluator;
using quantities::Angle;
using quantities::AngularFrequency;
using quantities::Cos;
using quantities::Length;
using quantities::Sin;
using quantities::Speed;
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
using testing_utilities::IsNear;
using ::testing::Sequence;
using ::testing::SetArgReferee;
using ::testing::_;

template<typename Frame>
class TestableContinuousTrajectory : public ContinuousTrajectory<Frame> {
 public:
  using ContinuousTrajectory<Frame>::ContinuousTrajectory;

  // Mock the Newhall factory.
  not_null<std::unique_ptr<Polynomial<Displacement<Frame>, Instant>>>
  NewhallApproximationInMonomialBasis(
      int degree,
      std::vector<Displacement<Frame>> const& q,
      std::vector<Velocity<Frame>> const& v,
      Instant const& t_min,
      Instant const& t_max,
      Displacement<Frame>& error_estimate) const override;

  MOCK_CONST_METHOD7_T(
      FillNewhallApproximationInMonomialBasis,
      void(int degree,
           std::vector<Displacement<Frame>> const& q,
           std::vector<Velocity<Frame>> const& v,
           Instant const& t_min,
           Instant const& t_max,
           Displacement<Frame>& error_estimate,
           not_null<std::unique_ptr<Polynomial<Displacement<Frame>, Instant>>>&
               polynomial));

  Status LockAndComputeBestNewhallApproximation(
      Instant const& time,
      std::vector<Displacement<Frame>> const& q,
      std::vector<Velocity<Frame>> const& v);

  // Helpers to access the internal state of the Newhall optimization.
  int degree() const;
  Length adjusted_tolerance() const;
  bool is_unstable() const;
  void ResetBestNewhallApproximation();
};

template<typename Frame>
not_null<std::unique_ptr<Polynomial<Displacement<Frame>, Instant>>>
TestableContinuousTrajectory<Frame>::NewhallApproximationInMonomialBasis(
    int degree,
    std::vector<Displacement<Frame>> const& q,
    std::vector<Velocity<Frame>> const& v,
    Instant const& t_min,
    Instant const& t_max,
    Displacement<Frame>& error_estimate) const {
  using P = PolynomialInMonomialBasis<
                Displacement<Frame>, Instant, /*degree=*/1, HornerEvaluator>;
  typename P::Coefficients const coefficients = {Displacement<Frame>(),
                                                 Velocity<Frame>()};
  not_null<std::unique_ptr<Polynomial<Displacement<Frame>, Instant>>>
      polynomial = make_not_null_unique<P>(coefficients, Instant());
  FillNewhallApproximationInMonomialBasis(degree,
                                          q, v,
                                          t_min, t_max,
                                          error_estimate,
                                          polynomial);
  return polynomial;
}

template<typename Frame>
Status
TestableContinuousTrajectory<Frame>::LockAndComputeBestNewhallApproximation(
    Instant const& time,
    std::vector<Displacement<Frame>> const& q,
    std::vector<Velocity<Frame>> const& v) {
  absl::MutexLock l(&lock_);
  return ComputeBestNewhallApproximation(time, q, v);
}

template<typename Frame>
int TestableContinuousTrajectory<Frame>::degree() const {
  return this->degree_;
}

template<typename Frame>
Length TestableContinuousTrajectory<Frame>::adjusted_tolerance() const {
  return this->adjusted_tolerance_;
}

template<typename Frame>
bool TestableContinuousTrajectory<Frame>::is_unstable() const {
  return this->is_unstable_;
}

template<typename Frame>
void TestableContinuousTrajectory<Frame>::ResetBestNewhallApproximation() {
  this->degree_age_ = std::numeric_limits<int>::max();
}

class ContinuousTrajectoryTest : public testing::Test {
 protected:
  using World = Frame<serialization::Frame::TestTag,
                      serialization::Frame::TEST1, true>;

  void FillTrajectory(
      int const number_of_steps,
      Time const& step,
      std::function<Position<World>(Instant const)> const& position_function,
      std::function<Velocity<World>(Instant const)> const& velocity_function,
      Instant const& time,
      ContinuousTrajectory<World>& trajectory) {
    for (int i = 0; i < number_of_steps; ++i) {
      // We use this way of computing the time (as opposed to consecutive
      // additions) because it results in a bit of jitter in the intervals,
      // which matters for continuity.
      Instant ti = time + (i + 1) * step;
      trajectory.Append(ti,
                        DegreesOfFreedom<World>(position_function(ti),
                                                velocity_function(ti)));
    }
  }

  Instant const t0_;
};

TEST_F(ContinuousTrajectoryTest, BestNewhallApproximation) {
  Time const step = 1 * Second;
  Length const tolerance = 1 * Metre;
  Instant t = t0_;
  std::vector<Displacement<World>> const q;
  std::vector<Velocity<World>> const v;

  auto const trajectory = std::make_unique<TestableContinuousTrajectory<World>>(
                              step,
                              tolerance);
  trajectory->Append(Instant(),
                     DegreesOfFreedom<World>(Position<World>(),
                                             Velocity<World>()));

  // A case where the errors smoothly decrease.
  {
    Sequence s;
    EXPECT_CALL(*trajectory,
                FillNewhallApproximationInMonomialBasis(3, _, _, _, _, _, _))
        .WillOnce(SetArgReferee<5>(
            Displacement<World>({3 * Metre, 4 * Metre, 5 * Metre})));
    EXPECT_CALL(*trajectory,
                FillNewhallApproximationInMonomialBasis(4, _, _, _, _, _, _))
        .WillOnce(SetArgReferee<5>(
            Displacement<World>({2 * Metre, 1 * Metre, 2 * Metre})));
    EXPECT_CALL(*trajectory,
                FillNewhallApproximationInMonomialBasis(5, _, _, _, _, _, _))
        .WillOnce(SetArgReferee<5>(
            Displacement<World>({0.1 * Metre, 2 * Metre, 0 * Metre})));
    EXPECT_CALL(*trajectory,
                FillNewhallApproximationInMonomialBasis(6, _, _, _, _, _, _))
        .WillOnce(SetArgReferee<5>(
            Displacement<World>({0.5 * Metre, 0.5 * Metre, 0.1 * Metre})));
    t += step;
    trajectory->LockAndComputeBestNewhallApproximation(t, q, v);
    EXPECT_EQ(6, trajectory->degree());
    EXPECT_EQ(tolerance, trajectory->adjusted_tolerance());
    EXPECT_FALSE(trajectory->is_unstable());
    trajectory->ResetBestNewhallApproximation();
  }

  // A case where the errors increase before we have reach the desired
  // tolerance...
  {
    Sequence s;
    EXPECT_CALL(*trajectory,
                FillNewhallApproximationInMonomialBasis(3, _, _, _, _, _, _))
        .WillOnce(SetArgReferee<5>(
            Displacement<World>({3 * Metre, 4 * Metre, 5 * Metre})));
    EXPECT_CALL(*trajectory,
                FillNewhallApproximationInMonomialBasis(4, _, _, _, _, _, _))
        .WillOnce(SetArgReferee<5>(
            Displacement<World>({2 * Metre, 1 * Metre, 2 * Metre})));
    EXPECT_CALL(*trajectory,
                FillNewhallApproximationInMonomialBasis(5, _, _, _, _, _, _))
        .WillOnce(SetArgReferee<5>(
            Displacement<World>({0.1 * Metre, 2 * Metre, 0 * Metre})));
    EXPECT_CALL(*trajectory,
                FillNewhallApproximationInMonomialBasis(6, _, _, _, _, _, _))
        .WillOnce(SetArgReferee<5>(
            Displacement<World>({1 * Metre, 3 * Metre, 1 * Metre})));
    t += step;
    trajectory->LockAndComputeBestNewhallApproximation(t, q, v);
    EXPECT_EQ(5, trajectory->degree());
    EXPECT_EQ(sqrt(4.01) * Metre, trajectory->adjusted_tolerance());
    EXPECT_TRUE(trajectory->is_unstable());
  }

  // ... then the error decreases...
  {
    Sequence s;
    EXPECT_CALL(*trajectory,
                FillNewhallApproximationInMonomialBasis(5, _, _, _, _, _, _))
        .WillOnce(SetArgReferee<5>(
            Displacement<World>({0.1 * Metre, 1.5 * Metre, 0 * Metre})));
    t += step;
    trajectory->LockAndComputeBestNewhallApproximation(t, q, v);
    EXPECT_EQ(5, trajectory->degree());
    EXPECT_EQ(sqrt(4.01) * Metre, trajectory->adjusted_tolerance());
    EXPECT_TRUE(trajectory->is_unstable());
  }

  // ... then the error increases forcing us to go back to square one...
  {
    Sequence s;
    EXPECT_CALL(*trajectory,
                FillNewhallApproximationInMonomialBasis(5, _, _, _, _, _, _))
        .WillOnce(SetArgReferee<5>(
            Displacement<World>({0.1 * Metre, 2 * Metre, 0.5 * Metre})))
        .WillOnce(SetArgReferee<5>(
            Displacement<World>({1 * Metre, 2 * Metre, 1 * Metre})));
    EXPECT_CALL(*trajectory,
                FillNewhallApproximationInMonomialBasis(3, _, _, _, _, _, _))
        .WillOnce(SetArgReferee<5>(
            Displacement<World>({3 * Metre, 4 * Metre, 5 * Metre})));
    EXPECT_CALL(*trajectory,
                FillNewhallApproximationInMonomialBasis(4, _, _, _, _, _, _))
        .WillOnce(SetArgReferee<5>(
            Displacement<World>({2 * Metre, 1 * Metre, 2 * Metre})));
    EXPECT_CALL(*trajectory,
                FillNewhallApproximationInMonomialBasis(6, _, _, _, _, _, _))
        .WillOnce(SetArgReferee<5>(
            Displacement<World>({1 * Metre, 1.5 * Metre, 1 * Metre})));
    EXPECT_CALL(*trajectory,
                FillNewhallApproximationInMonomialBasis(7, _, _, _, _, _, _))
        .WillOnce(SetArgReferee<5>(
            Displacement<World>({1 * Metre, 1.2 * Metre, 1 * Metre})));
    EXPECT_CALL(*trajectory,
                FillNewhallApproximationInMonomialBasis(8, _, _, _, _, _, _))
        .WillOnce(SetArgReferee<5>(
            Displacement<World>({1 * Metre, 1.3 * Metre, 1 * Metre})));
    t += step;
    trajectory->LockAndComputeBestNewhallApproximation(t, q, v);
    EXPECT_EQ(7, trajectory->degree());
    EXPECT_EQ(sqrt(3.44) * Metre, trajectory->adjusted_tolerance());
    EXPECT_TRUE(trajectory->is_unstable());
  }

  // ... it does it again but then the computation becomes stable.
  {
    Sequence s;
    EXPECT_CALL(*trajectory,
                FillNewhallApproximationInMonomialBasis(7, _, _, _, _, _, _))
        .WillOnce(SetArgReferee<5>(
            Displacement<World>({1 * Metre, 1.3 * Metre, 1 * Metre})));
    EXPECT_CALL(*trajectory,
                FillNewhallApproximationInMonomialBasis(3, _, _, _, _, _, _))
        .WillOnce(SetArgReferee<5>(
            Displacement<World>({3 * Metre, 4 * Metre, 5 * Metre})));
    EXPECT_CALL(*trajectory,
                FillNewhallApproximationInMonomialBasis(4, _, _, _, _, _, _))
        .WillOnce(SetArgReferee<5>(
            Displacement<World>({0.1 * Metre, 0.5 * Metre, 0.2 * Metre})));
    t += step;
    trajectory->LockAndComputeBestNewhallApproximation(t, q, v);
    EXPECT_EQ(4, trajectory->degree());
    EXPECT_EQ(tolerance, trajectory->adjusted_tolerance());
    EXPECT_FALSE(trajectory->is_unstable());
    trajectory->ResetBestNewhallApproximation();
  }

  // Check that the degree is properly lowered when the age of the approximation
  // exceeds the limit.
  // First, the errors force usage of degree 6.
  {
    Sequence s;
    EXPECT_CALL(*trajectory,
                FillNewhallApproximationInMonomialBasis(3, _, _, _, _, _, _))
        .WillOnce(SetArgReferee<5>(
            Displacement<World>({3 * Metre, 3 * Metre, 3 * Metre})));
    EXPECT_CALL(*trajectory,
                FillNewhallApproximationInMonomialBasis(4, _, _, _, _, _, _))
        .WillOnce(SetArgReferee<5>(
            Displacement<World>({2 * Metre, 2 * Metre, 2 * Metre})));
    EXPECT_CALL(*trajectory,
                FillNewhallApproximationInMonomialBasis(5, _, _, _, _, _, _))
        .WillOnce(SetArgReferee<5>(
            Displacement<World>({1 * Metre, 1 * Metre, 1 * Metre})));
    EXPECT_CALL(*trajectory,
                FillNewhallApproximationInMonomialBasis(6, _, _, _, _, _, _))
        .WillOnce(SetArgReferee<5>(
            Displacement<World>({0.1 * Metre, 0.1 * Metre, 0.1 * Metre})));
    t += step;
    trajectory->LockAndComputeBestNewhallApproximation(t, q, v);
    EXPECT_EQ(6, trajectory->degree());
    EXPECT_EQ(tolerance, trajectory->adjusted_tolerance());
    EXPECT_FALSE(trajectory->is_unstable());
  }

  // Then we get low errors for a long time.
  {
    Sequence s;
    EXPECT_CALL(*trajectory,
                FillNewhallApproximationInMonomialBasis(6, _, _, _, _, _, _))
        .Times(99)
        .WillRepeatedly(SetArgReferee<5>(
            Displacement<World>({0.1 * Metre, 0.1 * Metre, 0.1 * Metre})));
    for (int i = 0; i < 99; ++i) {
      t += step;
      trajectory->LockAndComputeBestNewhallApproximation(t, q, v);
    }
    EXPECT_EQ(6, trajectory->degree());
    EXPECT_EQ(tolerance, trajectory->adjusted_tolerance());
    EXPECT_FALSE(trajectory->is_unstable());
  }

  // Finally we try all the degrees again and discover that degree 5 works.
  {
    Sequence s;
    EXPECT_CALL(*trajectory,
                FillNewhallApproximationInMonomialBasis(3, _, _, _, _, _, _))
        .WillOnce(SetArgReferee<5>(
            Displacement<World>({3 * Metre, 3 * Metre, 3 * Metre})));
    EXPECT_CALL(*trajectory,
                FillNewhallApproximationInMonomialBasis(4, _, _, _, _, _, _))
        .WillOnce(SetArgReferee<5>(
            Displacement<World>({2 * Metre, 2 * Metre, 2 * Metre})));
    EXPECT_CALL(*trajectory,
                FillNewhallApproximationInMonomialBasis(5, _, _, _, _, _, _))
        .WillOnce(SetArgReferee<5>(
            Displacement<World>({0.2 * Metre, 0.2 * Metre, 0.2 * Metre})));
    t += step;
    trajectory->LockAndComputeBestNewhallApproximation(t, q, v);
    EXPECT_EQ(5, trajectory->degree());
    EXPECT_EQ(tolerance, trajectory->adjusted_tolerance());
    EXPECT_FALSE(trajectory->is_unstable());
    trajectory->ResetBestNewhallApproximation();
  }
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

  auto const trajectory = std::make_unique<ContinuousTrajectory<World>>(
                              step,
                              /*tolerance=*/0.1 * Metre);

  EXPECT_TRUE(trajectory->empty());
  FillTrajectory(number_of_steps,
                 step,
                 position_function,
                 velocity_function,
                 t0_,
                 *trajectory);
  EXPECT_FALSE(trajectory->empty());
  EXPECT_EQ(t0_ + step, trajectory->t_min());
  EXPECT_EQ(t0_ + (((number_of_steps - 1) / 8) * 8 + 1) * step,
            trajectory->t_max());

  for (Instant time = trajectory->t_min();
       time <= trajectory->t_max();
       time += step / number_of_substeps) {
    EXPECT_THAT(trajectory->EvaluatePosition(time) - World::origin,
                AlmostEquals(position_function(time) - World::origin, 0, 11));
    EXPECT_THAT(trajectory->EvaluateVelocity(time),
                AlmostEquals(velocity_function(time), 0, 4));
    EXPECT_EQ(trajectory->EvaluateDegreesOfFreedom(time),
              DegreesOfFreedom<World>(trajectory->EvaluatePosition(time),
                                      trajectory->EvaluateVelocity(time)));
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

  auto const trajectory = std::make_unique<ContinuousTrajectory<World>>(
                              step,
                              /*tolerance=*/5 * Milli(Metre));

  EXPECT_TRUE(trajectory->empty());
  FillTrajectory(number_of_steps,
                 step,
                 position_function,
                 velocity_function,
                 t0_,
                 *trajectory);
  EXPECT_FALSE(trajectory->empty());
  EXPECT_EQ(t0_ + step, trajectory->t_min());
  EXPECT_EQ(t0_ + (((number_of_steps - 1) / 8) * 8 + 1) * step,
            trajectory->t_max());

  Length max_position_absolute_error;
  Speed max_velocity_absolute_error;
  for (Instant time = trajectory->t_min();
       time <= trajectory->t_max();
       time += step / number_of_substeps) {
    Position<World> const actual_position = trajectory->EvaluatePosition(time);
    Position<World> const expected_position = position_function(time);
    Velocity<World> const actual_velocity =
        trajectory->EvaluateVelocity(time);
    Velocity<World> const expected_velocity = velocity_function(time);
    max_position_absolute_error =
        std::max(max_position_absolute_error,
                 AbsoluteError(expected_position, actual_position));
    max_velocity_absolute_error =
        std::max(max_velocity_absolute_error,
                 AbsoluteError(expected_velocity, actual_velocity));
  }
  EXPECT_THAT(max_position_absolute_error, IsNear(31 * Milli(Metre)));
  EXPECT_THAT(max_velocity_absolute_error, IsNear(1.45e-5 * Metre / Second));

  trajectory->ForgetBefore(trajectory->t_min() - step);

  Instant const forget_before_time = t0_ + 44444 * Second;
  trajectory->ForgetBefore(forget_before_time);
  EXPECT_EQ(forget_before_time, trajectory->t_min());
  EXPECT_EQ(t0_ + (((number_of_steps - 1) / 8) * 8 + 1) * step,
            trajectory->t_max());

  max_position_absolute_error = 0 * Metre;
  max_velocity_absolute_error = 0 * Metre / Second;
  for (Instant time = trajectory->t_min();
       time <= trajectory->t_max();
       time += step / number_of_substeps) {
    Position<World> const actual_position = trajectory->EvaluatePosition(time);
    Position<World> const expected_position = position_function(time);
    Velocity<World> const actual_velocity = trajectory->EvaluateVelocity(time);
    Velocity<World> const expected_velocity = velocity_function(time);
    max_position_absolute_error =
        std::max(max_position_absolute_error,
                 AbsoluteError(expected_position, actual_position));
    max_velocity_absolute_error =
        std::max(max_velocity_absolute_error,
                 AbsoluteError(expected_velocity, actual_velocity));
  }
  EXPECT_THAT(max_position_absolute_error, IsNear(31 * Milli(Metre)));
  EXPECT_THAT(max_velocity_absolute_error, IsNear(1.45e-5 * Metre / Second));
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

  auto const trajectory = std::make_unique<ContinuousTrajectory<World>>(
                              step,
                              /*tolerance=*/1 * Milli(Metre));

  EXPECT_TRUE(trajectory->empty());
  FillTrajectory(number_of_steps,
                 step,
                 position_function,
                 velocity_function,
                 t0_ + initial_time,
                 *trajectory);
  EXPECT_FALSE(trajectory->empty());

  // This time is exactly at the continuity point of two consecutive series.
  int const interval = 11;
  Instant const continuity_time =
      t0_ + initial_time + (8 * interval + 1) * step;

  Position<World> const p1 =
      trajectory->EvaluatePosition(continuity_time);
  Position<World> const p2 =
      trajectory->EvaluatePosition(continuity_time + step);
  Position<World> const p3 =
      trajectory->EvaluatePosition(continuity_time);
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

  auto const trajectory = std::make_unique<ContinuousTrajectory<World>>(
                              step, tolerance);

  EXPECT_TRUE(trajectory->empty());
  FillTrajectory(number_of_steps,
                 step,
                 position_function,
                 velocity_function,
                 t0_,
                 *trajectory);
  serialization::ContinuousTrajectory message;
  trajectory->WriteToMessage(&message);
  EXPECT_EQ(step / Second, message.step().magnitude());
  EXPECT_EQ(tolerance / Metre, message.tolerance().magnitude());
  EXPECT_GE(message.adjusted_tolerance().magnitude(),
            message.tolerance().magnitude());
  EXPECT_TRUE(message.has_is_unstable());
  EXPECT_EQ(3, message.degree());
  EXPECT_GE(100, message.degree_age());
  EXPECT_EQ(2, message.instant_polynomial_pair_size());
  EXPECT_TRUE(message.has_first_time());
  EXPECT_EQ(4, message.last_point_size());

  auto const trajectory_read =
      ContinuousTrajectory<World>::ReadFromMessage(message);
  EXPECT_EQ(trajectory->t_min(), trajectory_read->t_min());
  EXPECT_EQ(trajectory->t_max(), trajectory_read->t_max());
  for (Instant time = trajectory->t_min();
       time <= trajectory->t_max();
       time += step / number_of_substeps) {
    EXPECT_EQ(trajectory->EvaluateDegreesOfFreedom(time),
              trajectory_read->EvaluateDegreesOfFreedom(time));
  }

  serialization::ContinuousTrajectory second_message;
  trajectory_read->WriteToMessage(&second_message);
  EXPECT_THAT(message, EqualsProto(second_message));
}

TEST_F(ContinuousTrajectoryTest, PreCohenCompatibility) {
  Time const step = 0.01 * Second;
  Length const tolerance = 0.1 * Metre;

  // Fill the basic fields of a serialized ContinuousTrajectory.
  auto const trajectory = std::make_unique<ContinuousTrajectory<World>>(
                              step, tolerance);
  serialization::ContinuousTrajectory message;
  trajectory->WriteToMessage(&message);

  // Remove the polynomials and add a single Чебышёв series of the form:
  //   T₀ - 2 * T₁ + 3 * T₂  + 4 * T₃.
  message.clear_instant_polynomial_pair();
  auto* const series = message.add_series();
  Instant t_min = Instant() - 1 * Second;
  Instant t_max = Instant() + 1 * Second;
  t_min.WriteToMessage(series->mutable_t_min());
  t_max.WriteToMessage(series->mutable_t_max());
  Displacement<World> const c0({1.0 * Metre, 1.0 * Metre, 1.0 * Metre});
  Displacement<World> const c1({-2.0 * Metre, -2.0 * Metre, -2.0 * Metre});
  Displacement<World> const c2({3.0 * Metre, 3.0 * Metre, 3.0 * Metre});
  Displacement<World> const c3({4.0 * Metre, 4.0 * Metre, 4.0 * Metre});
  c0.WriteToMessage(series->add_coefficient()->mutable_multivector());
  c1.WriteToMessage(series->add_coefficient()->mutable_multivector());
  c2.WriteToMessage(series->add_coefficient()->mutable_multivector());
  c3.WriteToMessage(series->add_coefficient()->mutable_multivector());

  // Deserialize the message and check that a polynomial was constructed and
  // that it has the form:
  //   -2 - 14 * t + 6 * t^2 + 16 * t^3.
  auto const trajectory_read =
      ContinuousTrajectory<World>::ReadFromMessage(message);
  serialization::ContinuousTrajectory message2;
  trajectory_read->WriteToMessage(&message2);
  EXPECT_EQ(1, message2.instant_polynomial_pair_size());
  EXPECT_EQ(3, message2.instant_polynomial_pair(0).polynomial().degree());
  auto const& polynomial_in_monomial_basis =
      message2.instant_polynomial_pair(0).polynomial().GetExtension(
          serialization::PolynomialInMonomialBasis::extension);
  EXPECT_EQ(-2,
            polynomial_in_monomial_basis.coefficient(0).multivector().vector().
                x().quantity().magnitude());
  EXPECT_EQ(-14,
            polynomial_in_monomial_basis.coefficient(1).multivector().vector().
                x().quantity().magnitude());
  EXPECT_EQ(6,
            polynomial_in_monomial_basis.coefficient(2).multivector().vector().
                x().quantity().magnitude());
  EXPECT_EQ(16,
            polynomial_in_monomial_basis.coefficient(3).multivector().vector().
                x().quantity().magnitude());
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

  auto const trajectory = std::make_unique<ContinuousTrajectory<World>>(
                              step, tolerance);

  EXPECT_TRUE(trajectory->empty());

  // Fill the trajectory, get a checkpoint and fill some more.
  FillTrajectory(number_of_steps1,
                 step,
                 position_function,
                 velocity_function,
                 t0_,
                 *trajectory);
  ContinuousTrajectory<World>::Checkpoint const checkpoint =
      trajectory->GetCheckpoint();
  Instant const checkpoint_t_max = trajectory->t_max();
  FillTrajectory(number_of_steps2,
                 step,
                 position_function,
                 velocity_function,
                 t0_ + number_of_steps1 * step,
                 *trajectory);

  serialization::ContinuousTrajectory message;
  trajectory->WriteToMessage(&message, checkpoint);
  EXPECT_EQ(step / Second, message.step().magnitude());
  EXPECT_EQ(tolerance / Metre, message.tolerance().magnitude());
  EXPECT_GE(message.adjusted_tolerance().magnitude(),
            message.tolerance().magnitude());
  EXPECT_TRUE(message.has_is_unstable());
  EXPECT_EQ(3, message.degree());
  EXPECT_GE(100, message.degree_age());
  EXPECT_EQ(3, message.instant_polynomial_pair_size());
  EXPECT_TRUE(message.has_first_time());
  EXPECT_EQ(6, message.last_point_size());

  auto const trajectory_read =
      ContinuousTrajectory<World>::ReadFromMessage(message);
  EXPECT_EQ(trajectory_read->t_min(), trajectory->t_min());
  EXPECT_EQ(trajectory_read->t_max(), checkpoint_t_max);
  for (Instant time = trajectory->t_min();
       time <= checkpoint_t_max;
       time += step / number_of_substeps) {
    EXPECT_EQ(trajectory_read->EvaluateDegreesOfFreedom(time),
              trajectory->EvaluateDegreesOfFreedom(time));
  }
}

}  // namespace internal_continuous_trajectory
}  // namespace physics
}  // namespace principia
