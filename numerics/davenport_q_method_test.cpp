#include "numerics/davenport_q_method.hpp"

#include <random>
#include <vector>

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace numerics {

using geometry::Frame;
using geometry::Quaternion;
using geometry::Vector;
using quantities::Length;
using quantities::si::Metre;
using testing_utilities::AlmostEquals;

constexpr int number_of_test_vectors = 100;

class DavenportQMethodTest : public ::testing::Test {
 protected:
  using World1 = Frame<enum class World1Tag>;
  using World2 = Frame<enum class World2Tag>;

  DavenportQMethodTest() {
    std::mt19937_64 random(42);
    std::uniform_int_distribution<std::int64_t> coordinate_distribution(-5, 5);
    for (int i = 0; i < number_of_test_vectors; ++i) {
      double const x = coordinate_distribution(random);
      double const y = coordinate_distribution(random);
      double const z = coordinate_distribution(random);
      Vector<double, World1> const v({x, y, z});
      vectors1_.push_back(v / v.Norm());
      weights_.push_back(1 * Metre);
    }
  }

  std::vector<Vector<double, World1>> vectors1_;
  std::vector<Length> weights_;
};

TEST_F(DavenportQMethodTest, Identity) {
  EXPECT_THAT(DavenportQMethod(vectors1_, vectors1_, weights_),
              AlmostEquals(Quaternion(1), 0));
}

}  // namespace numerics
}  // namespace principia
