#include "numerics/davenport_q_method.hpp"

#include <random>
#include <vector>

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "geometry/r3_element.hpp"
#include "geometry/rotation.hpp"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/almost_equals.hpp"

namespace principia {
namespace numerics {

using geometry::Bivector;
using geometry::Frame;
using geometry::Quaternion;
using geometry::R3Element;
using geometry::Rotation;
using geometry::Vector;
using quantities::Length;
using quantities::si::Metre;
using quantities::si::Radian;
using testing_utilities::AlmostEquals;

constexpr int number_of_test_vectors = 100;

class DavenportQMethodTest : public ::testing::Test {
 protected:
  using World1 = Frame<struct World1Tag>;
  using World2 = Frame<struct World2Tag>;

  DavenportQMethodTest()
      : random_(42) {
    std::uniform_int_distribution<std::int64_t> coordinate_distribution(-5, 5);
    for (int i = 0; i < number_of_test_vectors; ++i) {
      double const x = coordinate_distribution(random_);
      double const y = coordinate_distribution(random_);
      double const z = coordinate_distribution(random_);
      Vector<double, World1> const v({x, y, z});
      vectors1_.push_back(v / v.Norm());
      weights_.push_back(1 * Metre);
    }
  }

  std::mt19937_64 random_;
  std::vector<Vector<double, World1>> vectors1_;
  std::vector<Length> weights_;
};

TEST_F(DavenportQMethodTest, Identity) {
  auto const rotation = Rotation<World1, World2>::Identity();
  EXPECT_THAT(DavenportQMethod(vectors1_, vectors1_, weights_),
              AlmostEquals(rotation, 0));
}

TEST_F(DavenportQMethodTest, FarFromIdentity) {
  Quaternion const q(2, R3Element<double>({1, -3, -2}));
  auto const normalized_q = q / q.Norm();
  Rotation<World1, World2> const rotation(normalized_q);

  std::vector<Vector<double, World2>> vectors2;
  for (auto const& vector1 : vectors1_) {
    vectors2.push_back(rotation(vector1));
  }

  EXPECT_THAT(DavenportQMethod(vectors1_, vectors2, weights_),
              AlmostEquals(rotation , 2, 16));
}

TEST_F(DavenportQMethodTest, Perturbed) {
  Quaternion const q(2, R3Element<double>({1, -3, -2}));

  std::vector<Vector<double, World2>> vectors2;
  std::uniform_real_distribution<double> quaternion_distribution(-1e-6, 1e-6);
  for (auto const& vector1 : vectors1_) {
    auto const perturbed_q =
        q + Quaternion(quaternion_distribution(random_),
                       R3Element<double>({quaternion_distribution(random_),
                                          quaternion_distribution(random_),
                                          quaternion_distribution(random_)}));
    Rotation<World1, World2> const perturbed_rotation(perturbed_q /
                                                      perturbed_q.Norm());
    vectors2.push_back(perturbed_rotation(vector1));
  }

  EXPECT_THAT(
      DavenportQMethod(vectors1_, vectors2, weights_),
      AlmostEquals(
          Rotation<World1, World2>(q / q.Norm()), 361'747'092, 745'100'359));
}

TEST_F(DavenportQMethodTest, PerturbedWeighted) {
  Quaternion const q(2, R3Element<double>({1, -3, -2}));
  auto const normalized_q = q / q.Norm();
  Rotation<World1, World2> const rotation(normalized_q);

  std::vector<Vector<double, World2>> vectors2;
  std::vector<double> weights;

  // The first entry is unperturbed and has a greater weight.
  vectors2.push_back(rotation(vectors1_[0]));
  weights.push_back(100);

  std::uniform_real_distribution<double> quaternion_distribution(-1e-6, 1e-6);
  for (int i = 1; i < vectors1_.size(); ++i) {
    auto const& vector1 = vectors1_[i];
    auto const perturbed_q =
        q + Quaternion(quaternion_distribution(random_),
                       R3Element<double>({quaternion_distribution(random_),
                                          quaternion_distribution(random_),
                                          quaternion_distribution(random_)}));
    Rotation<World1, World2> const perturbed_rotation(perturbed_q /
                                                      perturbed_q.Norm());
    vectors2.push_back(perturbed_rotation(vector1));
    weights.push_back(1);
  }

  EXPECT_THAT(DavenportQMethod(vectors1_, vectors2, weights),
              AlmostEquals(rotation, 112'728'156, 473'381'658));
}

}  // namespace numerics
}  // namespace principia
