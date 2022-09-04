#include "numerics/nearest_neighbour.hpp"

#include <random>
#include <vector>

#include "base/not_null.hpp"
#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "quantities/quantities.hpp"

namespace principia {
namespace numerics {

using base::not_null;
using geometry::Frame;
using geometry::Vector;
using quantities::Infinity;
using ::testing::Pointee;

class PrincipalComponentPartitioningTreeTest : public ::testing::Test {
 protected:
  using World = Frame<enum class WorldTag>;

  using V = Vector<double, World>;
};

TEST_F(PrincipalComponentPartitioningTreeTest, XYPlane) {
  V const v1({-1, -1, 0});
  V const v2({-1, 1, 0});
  V const v3({1, -1, -0});
  V const v4({2, 2, 0});

  // The centroid is at (1/4, 1/4, 0).  The first principal axis is
  // {1/√2, 1/√2, 0} and separates (v1, v2) from (v3, v4).  (Note that v2 and v3
  // are on the separator plane.)  The other principal axes are more complicated
  // because of the position of the centroid, but the first query entails an
  // "other side" search on (v1, v2).
  PrincipalComponentPartitioningTree<V> tree({&v1, &v2, &v3, &v4},
                                             /*max_values_per_cell=*/1);

  EXPECT_THAT(tree.FindNearestNeighbour(V({-0.5, -0.5, 0})), Pointee(v1));
  EXPECT_THAT(tree.FindNearestNeighbour(V({-0.5, 0.5, 0})), Pointee(v2));
  EXPECT_THAT(tree.FindNearestNeighbour(V({1, 1, 0})), Pointee(v4));
  EXPECT_THAT(tree.FindNearestNeighbour(V({3, 3, 0})), Pointee(v4));

  // Check that the points of the tree can be retrieved.
  EXPECT_THAT(tree.FindNearestNeighbour(v1), Pointee(v1));
  EXPECT_THAT(tree.FindNearestNeighbour(v2), Pointee(v2));
  EXPECT_THAT(tree.FindNearestNeighbour(v3), Pointee(v3));
  EXPECT_THAT(tree.FindNearestNeighbour(v4), Pointee(v4));
}

TEST_F(PrincipalComponentPartitioningTreeTest, Random) {
  static constexpr int points_in_tree = 100;
  static constexpr int points_to_test = 100;
  std::mt19937_64 random(42);
  std::uniform_real_distribution<double> coordinate_distribution(-10, 10);

  // Build two trees with the same points but different leaf sizes.
  std::vector<V> tree_points;
  std::vector<not_null<V const*>> tree_pointers;
  tree_points.reserve(points_in_tree);  // To avoid pointer invalidation below.
  for (int i = 0; i < points_in_tree; ++i) {
    tree_points.push_back(V({coordinate_distribution(random),
                             coordinate_distribution(random),
                             coordinate_distribution(random)}));
    tree_pointers.push_back(&tree_points.back());
  }
  PrincipalComponentPartitioningTree<V> tree1(tree_pointers,
                                              /*max_values_per_cell=*/1);
  PrincipalComponentPartitioningTree<V> tree3(tree_pointers,
                                              /*max_values_per_cell=*/3);

  for (int i = 0; i < points_to_test; ++i) {
    auto const query_point = V({coordinate_distribution(random),
                                coordinate_distribution(random),
                                coordinate_distribution(random)});

    // Compute the nearest point using the naive algorithm.
    V const* nearest = nullptr;
    double nearest_distance = Infinity<double>;
    for (auto const& tree_point : tree_points) {
      double const distance = (tree_point - query_point).Norm();
      if (distance < nearest_distance) {
        nearest_distance = distance;
        nearest = &tree_point;
      }
    }

    auto const nearest1 = tree1.FindNearestNeighbour(query_point);
    auto const nearest3 = tree3.FindNearestNeighbour(query_point);

    EXPECT_EQ(nearest1, nearest);
    EXPECT_EQ(nearest3, nearest);
  }
}

}  // namespace numerics
}  // namespace principia
