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
#include "testing_utilities/approximate_quantity.hpp"
#include "testing_utilities/is_near.hpp"
#include "testing_utilities/matchers.hpp"
#include "testing_utilities/numerics.hpp"

namespace principia {
namespace physics {
namespace internal_continuous_trajectory {

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
using testing_utilities::operator""_;
using ::testing::Sequence;
using ::testing::SetArgReferee;
using ::testing::_;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_named_quantities;

template<typename Frame>
class TestableContinuousTrajectory : public ContinuousTrajectory<Frame> {
 public:
  using ContinuousTrajectory<Frame>::ContinuousTrajectory;

  // Fake the Newhall factory.
  not_null<std::unique_ptr<Polynomial<Position<Frame>, Instant>>>
  NewhallApproximationInMonomialBasis(
      int degree,
      std::vector<Position<Frame>> const& q,
      std::vector<Velocity<Frame>> const& v,
      Instant const& t_min,
      Instant const& t_max,
      Displacement<Frame>& error_estimate) const override;

  MOCK_METHOD(
      void,
      FillNewhallApproximationInMonomialBasis,
      (int degree,
       std::vector<Position<Frame>> const& q,
       std::vector<Velocity<Frame>> const& v,
       Instant const& t_min,
       Instant const& t_max,
       Displacement<Frame>& error_estimate,
       (not_null<std::unique_ptr<Polynomial<Position<Frame>, Instant>>> &
        polynomial)),
      (const));

  absl::Status LockAndComputeBestNewhallApproximation(
      Instant const& time,
      std::vector<Position<Frame>> const& q,
      std::vector<Velocity<Frame>> const& v);

  // Helpers to access the internal state of the Newhall optimization.
  int degree() const;
  Length adjusted_tolerance() const;
  bool is_unstable() const;
  void ResetBestNewhallApproximation();
};

template<typename Frame>
not_null<std::unique_ptr<Polynomial<Position<Frame>, Instant>>>
TestableContinuousTrajectory<Frame>::NewhallApproximationInMonomialBasis(
    int degree,
    std::vector<Position<Frame>> const& q,
    std::vector<Velocity<Frame>> const& v,
    Instant const& t_min,
    Instant const& t_max,
    Displacement<Frame>& error_estimate) const {
  using P = PolynomialInMonomialBasis<
                Position<Frame>, Instant, /*degree=*/1, HornerEvaluator>;
  typename P::Coefficients const coefficients = {Position<Frame>(),
                                                 Velocity<Frame>()};
  not_null<std::unique_ptr<Polynomial<Position<Frame>, Instant>>>
      polynomial = make_not_null_unique<P>(coefficients, Instant());
  FillNewhallApproximationInMonomialBasis(degree,
                                          q, v,
                                          t_min, t_max,
                                          error_estimate,
                                          polynomial);
  return polynomial;
}

template<typename Frame>
absl::Status
TestableContinuousTrajectory<Frame>::LockAndComputeBestNewhallApproximation(
    Instant const& time,
    std::vector<Position<Frame>> const& q,
    std::vector<Velocity<Frame>> const& v) {
  absl::MutexLock l(&this->lock_);
  return this->ComputeBestNewhallApproximation(time, q, v);
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
                      Inertial,
                      Handedness::Right,
                      serialization::Frame::TEST>;

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
      EXPECT_OK(
          trajectory.Append(ti,
                            DegreesOfFreedom<World>(position_function(ti),
                                                    velocity_function(ti))));
    }
  }

  // Fills the pre-Grassmann fields and clear the post-Grassmann fields.
  serialization::ContinuousTrajectory MakePreGrassmann(
      serialization::ContinuousTrajectory const& message,
      std::optional<Instant> const& checkpoint_time) {
    serialization::ContinuousTrajectory pre_grassmann = message;
    EXPECT_EQ(1, pre_grassmann.checkpoint_size());
    auto checkpoint = pre_grassmann.checkpoint(0);
    *pre_grassmann.mutable_adjusted_tolerance() =
        checkpoint.adjusted_tolerance();
    pre_grassmann.set_is_unstable(checkpoint.is_unstable());
    pre_grassmann.set_degree(checkpoint.degree());
    pre_grassmann.set_degree_age(checkpoint.degree_age());
    pre_grassmann.mutable_last_point()->Swap(checkpoint.mutable_last_point());
    if (checkpoint_time.has_value()) {
      // This is post-Fatou.
      checkpoint_time->WriteToMessage(pre_grassmann.mutable_checkpoint_time());
    }
    pre_grassmann.clear_checkpoint();
    return pre_grassmann;
  }

  // Fills the pre-Gröbner fields and clear the post-Gröbner fields.
  serialization::ContinuousTrajectory MakePreGröbner(
      serialization::ContinuousTrajectory const& message) {
    serialization::ContinuousTrajectory pre_gröbner = message;
    for (int i = 0; i < pre_gröbner.instant_polynomial_pair_size(); ++i) {
      auto const coefficient0 = pre_gröbner.instant_polynomial_pair(i).
          polynomial().
          GetExtension(serialization::PolynomialInMonomialBasis::extension).
              coefficient(0).point().multivector();
      *pre_gröbner.mutable_instant_polynomial_pair(i)
           ->mutable_polynomial()->MutableExtension(
               serialization::PolynomialInMonomialBasis::extension)
           ->mutable_coefficient(0)->mutable_multivector() = coefficient0;
    }
    return pre_gröbner;
  }

  Instant const t0_;
};

TEST_F(ContinuousTrajectoryTest, BestNewhallApproximation) {
  Time const step = 1 * Second;
  Length const tolerance = 1 * Metre;
  Instant t = t0_;
  std::vector<Position<World>> const q;
  std::vector<Velocity<World>> const v;

  auto const trajectory = std::make_unique<TestableContinuousTrajectory<World>>(
                              step,
                              tolerance);
  EXPECT_OK(trajectory->Append(
      Instant(), DegreesOfFreedom<World>(World::origin, World::unmoving)));

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
    EXPECT_OK(trajectory->LockAndComputeBestNewhallApproximation(t, q, v));
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
    EXPECT_OK(trajectory->LockAndComputeBestNewhallApproximation(t, q, v));
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
    EXPECT_OK(trajectory->LockAndComputeBestNewhallApproximation(t, q, v));
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
    EXPECT_OK(trajectory->LockAndComputeBestNewhallApproximation(t, q, v));
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
    EXPECT_OK(trajectory->LockAndComputeBestNewhallApproximation(t, q, v));
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
    EXPECT_OK(trajectory->LockAndComputeBestNewhallApproximation(t, q, v));
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
      EXPECT_OK(trajectory->LockAndComputeBestNewhallApproximation(t, q, v));
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
    EXPECT_OK(trajectory->LockAndComputeBestNewhallApproximation(t, q, v));
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

  // Check that the positions and velocities match the ones given by the
  // functions above.
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

#if PRINCIPIA_CONTINUOUS_TRAJECTORY_SUPPORTS_PIECEWISE_POISSON_SERIES
  // Now check that it can be converted to a piecewise Poisson series.
  Instant const t_min = trajectory->t_min() + 2 * step / number_of_substeps;
  Instant const t_max = trajectory->t_max() - 3 * step / number_of_substeps;
  EXPECT_EQ(3, trajectory->PiecewisePoissonSeriesDegree(t_min, t_max));
  auto const piecewise_poisson_series =
      trajectory->ToPiecewisePoissonSeries<3, 0>(t_min, t_max);
  EXPECT_EQ(t_min, piecewise_poisson_series.t_min());
  EXPECT_EQ(t_max, piecewise_poisson_series.t_max());
  for (Instant time = t_min; time <= t_max; time += step / number_of_substeps) {
    EXPECT_THAT(
        piecewise_poisson_series(time),
        AlmostEquals(trajectory->EvaluatePosition(time) - World::origin, 0, 0));
  }
#endif
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
  EXPECT_THAT(max_position_absolute_error, IsNear(31_(1) * Milli(Metre)));
  EXPECT_THAT(max_velocity_absolute_error,
              IsNear(1.40e-5_(1) * Metre / Second));
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
  [[maybe_unused]] Position<World> const p2 =
      trajectory->EvaluatePosition(continuity_time + step);
  Position<World> const p3 =
      trajectory->EvaluatePosition(continuity_time);
  EXPECT_THAT(p1, AlmostEquals(p3, 0, 2));
}

TEST_F(ContinuousTrajectoryTest, Prepend) {
  int const number_of_steps1 = 20;
  int const number_of_steps2 = 15;
  int const number_of_substeps = 50;
  Time const step = 0.01 * Second;
  Length const tolerance = 0.1 * Metre;

  // Construct two trajectories with different functions.

  Instant const t1 = t0_;
  auto position_function1 =
      [t1](Instant const t) {
        return World::origin +
            Displacement<World>({(t - t1) * 3 * Metre / Second,
                                 (t - t1) * 5 * Metre / Second,
                                 (t - t1) * (-2) * Metre / Second});
      };
  auto velocity_function1 =
      [](Instant const t) {
        return Velocity<World>({3 * Metre / Second,
                                5 * Metre / Second,
                                -2 * Metre / Second});
      };
  auto trajectory1 =
      std::make_unique<ContinuousTrajectory<World>>(step, tolerance);
  FillTrajectory(number_of_steps1,
                 step,
                 position_function1,
                 velocity_function1,
                 t1,
                 *trajectory1);
  EXPECT_EQ(t1 + step, trajectory1->t_min());
  EXPECT_EQ(t1 + (((number_of_steps1 - 1) / 8) * 8 + 1) * step,
            trajectory1->t_max());

  Instant const t2 = trajectory1->t_max();
  auto position_function2 =
      [&position_function1, t2](Instant const t) {
        return position_function1(t2) +
               Displacement<World>({(t - t2) * 6 * Metre / Second,
                                    (t - t2) * 1.5 * Metre / Second,
                                    (t - t2) * 7 * Metre / Second});
      };
  auto velocity_function2 =
      [](Instant const t) {
        return Velocity<World>({6 * Metre / Second,
                                1.5 * Metre / Second,
                                7 * Metre / Second});
      };
  auto trajectory2 =
      std::make_unique<ContinuousTrajectory<World>>(step, tolerance);
  FillTrajectory(number_of_steps2 + 1,
                 step,
                 position_function2,
                 velocity_function2,
                 t2 - step,  // First point at t2.
                 *trajectory2);
  EXPECT_EQ(t2, trajectory2->t_min());
  EXPECT_EQ(t2 + (number_of_steps2 / 8) * 8 * step,
            trajectory2->t_max());

  // Prepend one trajectory to the other.
  trajectory2->Prepend(std::move(*trajectory1));

  // Verify the resulting trajectory.
  EXPECT_EQ(t1 + step, trajectory2->t_min());
  EXPECT_EQ(t2 + (number_of_steps2 / 8) * 8 * step,
            trajectory2->t_max());
  for (Instant time = trajectory2->t_min();
       time <= t2;
       time += step / number_of_substeps) {
    EXPECT_THAT(trajectory2->EvaluatePosition(time),
                AlmostEquals(position_function1(time), 0, 10)) << time;
    EXPECT_THAT(trajectory2->EvaluateVelocity(time),
                AlmostEquals(velocity_function1(time), 0, 4)) << time;
  }
  for (Instant time = t2 + step / number_of_substeps;
       time <= trajectory2->t_max();
       time += step / number_of_substeps) {
    EXPECT_THAT(trajectory2->EvaluatePosition(time),
                AlmostEquals(position_function2(time), 0, 2842)) << time;
    EXPECT_THAT(trajectory2->EvaluateVelocity(time),
                AlmostEquals(velocity_function2(time), 0, 34)) << time;
  }
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

  // Take a checkpoint and verify that the checkpointed data is properly
  // serialized.
  trajectory->WriteToCheckpoint(trajectory->t_max());
  serialization::ContinuousTrajectory message;
  trajectory->WriteToMessage(&message);
  EXPECT_EQ(step / Second, message.step().magnitude());
  EXPECT_EQ(tolerance / Metre, message.tolerance().magnitude());
  EXPECT_EQ(2, message.instant_polynomial_pair_size());
  EXPECT_TRUE(message.has_first_time());

  EXPECT_EQ(1, message.checkpoint_size());
  EXPECT_FALSE(message.has_adjusted_tolerance());
  EXPECT_FALSE(message.has_is_unstable());
  EXPECT_FALSE(message.has_degree());
  EXPECT_FALSE(message.has_degree_age());
  EXPECT_EQ(0, message.last_point_size());

  auto const& checkpoint = message.checkpoint(0);
  EXPECT_GE(checkpoint.adjusted_tolerance().magnitude(),
            message.tolerance().magnitude());
  EXPECT_TRUE(checkpoint.has_is_unstable());
  EXPECT_EQ(3, checkpoint.degree());
  EXPECT_GE(100, checkpoint.degree_age());
  EXPECT_EQ(4, checkpoint.last_point_size());

  auto const trajectory_read = ContinuousTrajectory<World>::ReadFromMessage(
      /*desired_t_min=*/InfiniteFuture,
      message);
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
  int const number_of_steps = 110;
  Time const step = 0.01 * Second;
  Length const tolerance = 0.1 * Metre;

  // Fill a ContinuousTrajectory and take a checkpoint.
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
  FillTrajectory(number_of_steps,
                 step,
                 position_function,
                 velocity_function,
                 t0_,
                 *trajectory);
  trajectory->WriteToCheckpoint(trajectory->t_max());

  serialization::ContinuousTrajectory message;
  trajectory->WriteToMessage(&message);

  // Revert the message to pre-Fatou, pre-Grassmann.
  message = MakePreGrassmann(message, /*checkpoint_time=*/std::nullopt);

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
  auto const trajectory_read = ContinuousTrajectory<World>::ReadFromMessage(
      /*desired_t_min=*/InfiniteFuture,
      message);
  serialization::ContinuousTrajectory message2;
  trajectory_read->WriteToMessage(&message2);
  EXPECT_EQ(1, message2.instant_polynomial_pair_size());
  EXPECT_EQ(3, message2.instant_polynomial_pair(0).polynomial().degree());
  auto const& polynomial_in_monomial_basis =
      message2.instant_polynomial_pair(0).polynomial().GetExtension(
          serialization::PolynomialInMonomialBasis::extension);
  EXPECT_EQ(-2,
            polynomial_in_monomial_basis.coefficient(0).point().multivector().
                vector().x().quantity().magnitude());
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

TEST_F(ContinuousTrajectoryTest, PreGrassmannCompatibility) {
  int const number_of_steps = 30;
  Time const step = 0.01 * Second;
  Length const tolerance = 0.1 * Metre;

  // Fill a ContinuousTrajectory and take a checkpoint.
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

  auto const trajectory1 = std::make_unique<ContinuousTrajectory<World>>(
                              step, tolerance);
  FillTrajectory(number_of_steps,
                 step,
                 position_function,
                 velocity_function,
                 t0_,
                 *trajectory1);
  Instant const checkpoint_time = trajectory1->t_max();
  trajectory1->WriteToCheckpoint(checkpoint_time);

  serialization::ContinuousTrajectory message1;
  trajectory1->WriteToMessage(&message1);

  serialization::ContinuousTrajectory const pre_grassmann =
      MakePreGrassmann(MakePreGröbner(message1), checkpoint_time);

  // Read from the pre-Grassmann message, write to a second message, and check
  // that we get the same result.
  auto const trajectory2 = ContinuousTrajectory<World>::ReadFromMessage(
      /*desired_t_min=*/InfiniteFuture,
      pre_grassmann);
  serialization::ContinuousTrajectory message2;
  trajectory2->WriteToMessage(&message2);

  EXPECT_THAT(message2, EqualsProto(message1));
}

TEST_F(ContinuousTrajectoryTest, PreGröbnerCompatibility) {
  int const number_of_steps = 30;
  Time const step = 0.01 * Second;
  Length const tolerance = 0.1 * Metre;

  // Fill a ContinuousTrajectory and take a checkpoint.
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

  auto const trajectory1 = std::make_unique<ContinuousTrajectory<World>>(
                              step, tolerance);
  FillTrajectory(number_of_steps,
                 step,
                 position_function,
                 velocity_function,
                 t0_,
                 *trajectory1);
  Instant const checkpoint_time = trajectory1->t_max();
  trajectory1->WriteToCheckpoint(checkpoint_time);

  serialization::ContinuousTrajectory message1;
  trajectory1->WriteToMessage(&message1);

  serialization::ContinuousTrajectory const pre_gröbner =
      MakePreGröbner(message1);

  // Read from the pre-Gröbner message, write to a second message, and check
  // that we get the same result.
  auto const trajectory2 = ContinuousTrajectory<World>::ReadFromMessage(
      /*desired_t_min=*/InfiniteFuture,
      pre_gröbner);
  serialization::ContinuousTrajectory message2;
  trajectory2->WriteToMessage(&message2);

  EXPECT_THAT(message2, EqualsProto(message1));
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

  // Fill the trajectory, create a checkpoint and fill some more.
  FillTrajectory(number_of_steps1,
                 step,
                 position_function,
                 velocity_function,
                 t0_,
                 *trajectory);
  EXPECT_EQ(t0_ + (((number_of_steps1 - 1) / 8) * 8 + 1) * step,
            trajectory->t_max());
  Instant const checkpoint_time = trajectory->t_max();
  trajectory->WriteToCheckpoint(checkpoint_time);
  FillTrajectory(number_of_steps2,
                 step,
                 position_function,
                 velocity_function,
                 t0_ + number_of_steps1 * step,
                 *trajectory);
  EXPECT_EQ(
      t0_ + (((number_of_steps1 + number_of_steps2 - 1) / 8) * 8 + 1) * step,
      trajectory->t_max());

  serialization::ContinuousTrajectory message;
  trajectory->WriteToMessage(&message);
  EXPECT_EQ(step / Second, message.step().magnitude());
  EXPECT_EQ(tolerance / Metre, message.tolerance().magnitude());
  EXPECT_EQ(1, message.checkpoint_size());
  EXPECT_EQ(3, message.instant_polynomial_pair_size());
  EXPECT_TRUE(message.has_first_time());

  auto const& checkpoint = message.checkpoint(0);
  EXPECT_GE(checkpoint.adjusted_tolerance().magnitude(),
            message.tolerance().magnitude());
  EXPECT_TRUE(checkpoint.has_is_unstable());
  EXPECT_EQ(3, checkpoint.degree());
  EXPECT_GE(100, checkpoint.degree_age());
  EXPECT_EQ(6, checkpoint.last_point_size());

  // Read the trajectory and check that everything is identical up to the
  // checkpoint.
  auto const trajectory_read = ContinuousTrajectory<World>::ReadFromMessage(
      /*desired_t_min=*/InfiniteFuture,
      message);
  EXPECT_EQ(trajectory_read->t_min(), trajectory->t_min());
  EXPECT_EQ(trajectory_read->t_max(), checkpoint_time);
  for (Instant time = trajectory->t_min();
       time <= checkpoint_time;
       time += step / number_of_substeps) {
    EXPECT_EQ(trajectory_read->EvaluateDegreesOfFreedom(time),
              trajectory->EvaluateDegreesOfFreedom(time));
  }

  // Extend the trajectory that was just read.
  FillTrajectory(number_of_steps2,
                 step,
                 position_function,
                 velocity_function,
                 t0_ + number_of_steps1 * step,
                 *trajectory_read);
  EXPECT_EQ(
      t0_ + (((number_of_steps1 + number_of_steps2 - 1) / 8) * 8 + 1) * step,
      trajectory->t_max());

  // Reset to the checkpoint and check that the polynomials were truncated.
  EXPECT_OK(trajectory_read->ReadFromCheckpointAt(
      checkpoint_time, trajectory_read->MakeCheckpointerReader()));
  EXPECT_EQ(trajectory_read->t_max(), checkpoint_time);
}

}  // namespace internal_continuous_trajectory
}  // namespace physics
}  // namespace principia
