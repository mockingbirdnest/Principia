
#include "ksp_plugin/equator_relevance_threshold.hpp"

#include <string>

#include "astronomy/frames.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "physics/solar_system.hpp"
#include "quantities/astronomy.hpp"
#include "testing_utilities/is_near.hpp"

namespace principia {
namespace ksp_plugin {
namespace internal_equator_relevance_threshold {

using astronomy::ICRS;
using base::not_null;
using physics::SolarSystem;
using quantities::astronomy::JovianEquatorialRadius;
using quantities::astronomy::SolarRadius;
using quantities::astronomy::TerrestrialEquatorialRadius;
using testing_utilities::IsNear;

class EquatorRelevanceThresholdTest : public testing::Test {
 protected:
  EquatorRelevanceThresholdTest()
      : solar_system_j2000_(
            SOLUTION_DIR / "astronomy" / "sol_gravity_model.proto.txt",
            SOLUTION_DIR / "astronomy" /
                "sol_initial_state_jd_2451545_000000000.proto.txt",
            /*ignore_frame=*/true){};

  not_null<std::unique_ptr<RotatingBody<Barycentric>>> MakeBody(
      std::string const& name) {
    return solar_system_j2000_.MakeRotatingBody(
        solar_system_j2000_.gravity_model_message(name));
  }

  Length mean_radius(std::string const& name) {
    return MakeBody(name)->mean_radius();
  }

  SolarSystem<Barycentric> solar_system_j2000_;
};

TEST_F(EquatorRelevanceThresholdTest, Planets) {
  // See the discussion on #1841.
  EXPECT_THAT(EquatorRelevanceThreshold(*MakeBody("Sun")),
              IsNear(58 * SolarRadius));
  EXPECT_THAT(EquatorRelevanceThreshold(*MakeBody("Mercury")),
              IsNear(158 * mean_radius("Mercury")));
  EXPECT_THAT(EquatorRelevanceThreshold(*MakeBody("Venus")),
              IsNear(403 * mean_radius("Venus")));
  EXPECT_THAT(EquatorRelevanceThreshold(*MakeBody("Earth")),
              IsNear(233 * TerrestrialEquatorialRadius));
  EXPECT_THAT(EquatorRelevanceThreshold(*MakeBody("Mars")),
              IsNear(314 * mean_radius("Mars")));
  EXPECT_THAT(EquatorRelevanceThreshold(*MakeBody("Jupiter")),
              IsNear(860 * JovianEquatorialRadius));
  EXPECT_THAT(EquatorRelevanceThreshold(*MakeBody("Saturn")),
              IsNear(938 * mean_radius("Saturn")));
  EXPECT_THAT(EquatorRelevanceThreshold(*MakeBody("Neptune")),
              IsNear(423 * mean_radius("Neptune")));
  EXPECT_THAT(EquatorRelevanceThreshold(*MakeBody("Uranus")),
              IsNear(424 * mean_radius("Uranus")));
}

}  // namespace internal_equator_relevance_threshold
}  // namespace ksp_plugin
}  // namespace principia
