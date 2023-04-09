#include "ksp_plugin/celestial.hpp"

#include "geometry/instant.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/massive_body.hpp"
#include "quantities/si.hpp"

namespace principia {
namespace ksp_plugin {

using namespace principia::base::_not_null;
using namespace principia::geometry::_instant;
using namespace principia::ksp_plugin::_celestial;
using namespace principia::physics::_continuous_trajectory;
using namespace principia::physics::_massive_body;
using namespace principia::physics::_rotating_body;
using namespace principia::quantities::_si;

class CelestialTest : public testing::Test {
 protected:
  CelestialTest()
      : body_(MassiveBody::Parameters(1 * Kilogram),
              RotatingBody<Barycentric>::Parameters(
                  /*mean_radius=*/1 * Metre,
                  /*reference_angle=*/0 * Degree,
                  /*reference_instant=*/astronomy::J2000,
                  /*angular_frequency=*/1 * Radian / Second,
                  /*right_ascension_of_pole=*/0 * Degree,
                  /*declination_of_pole=*/90 * Degree)),
        celestial_(make_not_null_unique<Celestial>(&body_)),
        trajectory_(1 * Second, 1 * Metre) {}

  RotatingBody<Barycentric> body_;
  not_null<std::unique_ptr<Celestial>> celestial_;
  ContinuousTrajectory<Barycentric> trajectory_;
};

using CelestialDeathTest = CelestialTest;

TEST_F(CelestialDeathTest, Uninitialized) {
  EXPECT_DEATH({celestial_->trajectory();}, "is_initialized");
  EXPECT_DEATH({celestial_->current_degrees_of_freedom(Instant());},
               "is_initialized");
  EXPECT_DEATH({celestial_->current_position(Instant());}, "is_initialized");
  EXPECT_DEATH({celestial_->current_velocity(Instant());}, "is_initialized");
}

TEST_F(CelestialDeathTest, OverlyInitialized) {
  EXPECT_DEATH({
    celestial_->set_trajectory(&trajectory_);
    celestial_->set_trajectory(&trajectory_);
  },
  "!is_initialized");
}

TEST_F(CelestialTest, Initialization) {
  EXPECT_FALSE(celestial_->is_initialized());
  ContinuousTrajectory<Barycentric> trajectory(1 * Second, 1 * Metre);
  celestial_->set_trajectory(&trajectory);
  EXPECT_TRUE(celestial_->is_initialized());
}

}  // namespace ksp_plugin
}  // namespace principia
