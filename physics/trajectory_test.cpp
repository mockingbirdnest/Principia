#include "trajectory.hpp"

#include "body.hpp"
#include "geometry/grassmann.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"

using principia::geometry::Vector;
using principia::quantities::Length;
using principia::quantities::Mass;
using principia::quantities::Speed;
using principia::quantities::Time;
using principia::quantities::SIUnit;
using principia::si::Metre;
using principia::si::Second;

namespace principia {
namespace physics {

class World;

class TrajectoryTest : public testing::Test {
 protected:
  void SetUp() override {
    body_.reset(new Body(SIUnit<Mass>()));
    trajectory_.reset(new Trajectory<World>(body_.get()));
  }

  std::unique_ptr<Body> body_;
  std::unique_ptr<Trajectory<World>> trajectory_;
};

TEST_F(TrajectoryTest, Append) {
  trajectory_->Append(Vector<Length, World>({1 * Metre, 2 * Metre, 3 * Metre}),
                      Vector<Speed, World>({4 * Metre / Second,
                                            5 * Metre / Second,
                                            6 * Metre / Second}),
                      7 * Second);
}

}  // namespace physics
}  // namespace principia
