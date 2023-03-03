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

using ::testing::Eq;
using ::testing::Pointee;
using namespace principia::base::_not_null;
using namespace principia::geometry::_frame;
using namespace principia::geometry::_grassmann;
using namespace principia::quantities::_quantities;

class PrincipalComponentPartitioningTreeTest : public ::testing::Test {
 protected:
  using World = Frame<struct WorldTag>;
  using V = Vector<double, World>;

  // Computes the nearest point using the brute force algorithm.
  V const* BruteForceNearestNeighbour(
      V const& query_value,
      std::vector<V> const& values,
      PrincipalComponentPartitioningTree<V>::Filter const& filter = nullptr) {
    V const* nearest = nullptr;
    double nearest_distance = Infinity<double>;
    for (auto const& value : values) {
      double const distance = (value - query_value).Norm();
      if (distance < nearest_distance &&
          (filter == nullptr || filter(&value))) {
        nearest_distance = distance;
        nearest = &value;
      }
    }
    return nearest;
  }

  // Fills the vectors with |number_of_values| randomly generated values.
  void MakeValues(
      int const number_of_values,
      std::vector<V>& values,
      std::vector<not_null<V const*>>& pointers,
      std::mt19937_64& random,
      std::uniform_real_distribution<double>& coordinate_distribution) {
    values.clear();
    pointers.clear();
    values.reserve(number_of_values);  // To avoid pointer invalidation below.
    for (int i = 0; i < number_of_values; ++i) {
      values.push_back(V({coordinate_distribution(random),
                          coordinate_distribution(random),
                          coordinate_distribution(random)}));
      pointers.push_back(&values.back());
    }
  }
};

// Points in the x-y plane.
TEST_F(PrincipalComponentPartitioningTreeTest, XYPlaneConstructor) {
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

// Same as the previous test, but the tree is initially empty and points are
// added using |Add|.
TEST_F(PrincipalComponentPartitioningTreeTest, XYPlaneAdd) {
  V const v1({-1, -1, 0});
  V const v2({-1, 1, 0});
  V const v3({1, -1, -0});
  V const v4({2, 2, 0});

  PrincipalComponentPartitioningTree<V> tree({},
                                             /*max_values_per_cell=*/1);
  tree.Add(&v1);
  tree.Add(&v2);
  tree.Add(&v3);
  tree.Add(&v4);

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

// Random points, validated against the brute force algorithm.
TEST_F(PrincipalComponentPartitioningTreeTest, RandomConstructor) {
  static constexpr int points_in_tree = 100;
  static constexpr int points_to_test = 100;
  std::mt19937_64 random(42);
  std::uniform_real_distribution<double> coordinate_distribution(-10, 10);

  // Build two trees with the same points but different leaf sizes.
  std::vector<V> tree_points;
  std::vector<not_null<V const*>> tree_pointers;
  MakeValues(points_in_tree,
             tree_points,
             tree_pointers,
             random,
             coordinate_distribution);
  PrincipalComponentPartitioningTree<V> tree1(tree_pointers,
                                              /*max_values_per_cell=*/1);
  PrincipalComponentPartitioningTree<V> tree3(tree_pointers,
                                              /*max_values_per_cell=*/3);

  for (int i = 0; i < points_to_test; ++i) {
    auto const query_point = V({coordinate_distribution(random),
                                coordinate_distribution(random),
                                coordinate_distribution(random)});

    auto* const nearest = BruteForceNearestNeighbour(query_point, tree_points);
    auto* const nearest1 = tree1.FindNearestNeighbour(query_point);
    auto* const nearest3 = tree3.FindNearestNeighbour(query_point);

    EXPECT_THAT(nearest1, Eq(nearest)) << *nearest1 << " " << *nearest;
    EXPECT_THAT(nearest3, Eq(nearest)) << *nearest3 << " " << *nearest;
  }
}

// Same as the previous test, but the trees are initially empty and points are
// added using |Add|.
TEST_F(PrincipalComponentPartitioningTreeTest, RandomAdd) {
  static constexpr int points_in_tree = 100;
  static constexpr int points_to_test = 100;
  std::mt19937_64 random(42);
  std::uniform_real_distribution<double> coordinate_distribution(-10, 10);

  // Build two trees with the same points but different leaf sizes.
  std::vector<V> tree_points;
  std::vector<not_null<V const*>> tree_pointers;
  MakeValues(points_in_tree,
             tree_points,
             tree_pointers,
             random,
             coordinate_distribution);
  PrincipalComponentPartitioningTree<V> tree1({},
                                              /*max_values_per_cell=*/1);
  PrincipalComponentPartitioningTree<V> tree3({},
                                              /*max_values_per_cell=*/3);
  for (auto const tree_pointer : tree_pointers) {
    tree1.Add(tree_pointer);
    tree3.Add(tree_pointer);
  }

  for (int i = 0; i < points_to_test; ++i) {
    auto const query_point = V({coordinate_distribution(random),
                                coordinate_distribution(random),
                                coordinate_distribution(random)});

    auto* const nearest = BruteForceNearestNeighbour(query_point, tree_points);
    auto* const nearest1 = tree1.FindNearestNeighbour(query_point);
    auto* const nearest3 = tree3.FindNearestNeighbour(query_point);

    EXPECT_THAT(nearest1, Eq(nearest)) << *nearest1 << " " << *nearest;
    EXPECT_THAT(nearest3, Eq(nearest)) << *nearest3 << " " << *nearest;
  }
}

// Random points with a filter, validated against the brute force algorithm.
TEST_F(PrincipalComponentPartitioningTreeTest, RandomFilter) {
  static constexpr int points_in_tree = 100;
  static constexpr int points_to_test = 10;
  std::mt19937_64 random(42);
  std::uniform_real_distribution<double> coordinate_distribution(-10, 10);

  // Build two trees with the same points but different leaf sizes.
  std::vector<V> tree_points;
  std::vector<not_null<V const*>> tree_pointers;
  MakeValues(points_in_tree,
             tree_points,
             tree_pointers,
             random,
             coordinate_distribution);
  PrincipalComponentPartitioningTree<V> tree1(tree_pointers,
                                              /*max_values_per_cell=*/1);
  PrincipalComponentPartitioningTree<V> tree3(tree_pointers,
                                              /*max_values_per_cell=*/3);

  const PrincipalComponentPartitioningTree<V>::Filter filter =
      [](V const* const point) {
    return point->Norm²() < 100;
  };

  bool filtering_was_effective = false;
  for (int i = 0; i < points_to_test; ++i) {
    auto const query_point = V({coordinate_distribution(random),
                                coordinate_distribution(random),
                                coordinate_distribution(random)});

    auto* const nearest =
        BruteForceNearestNeighbour(query_point, tree_points, filter);
    auto* const nearest1 = tree1.FindNearestNeighbour(query_point, filter);
    auto* const nearest3 = tree3.FindNearestNeighbour(query_point, filter);

    EXPECT_THAT(nearest1, Eq(nearest)) << *nearest1 << " " << *nearest;
    EXPECT_THAT(nearest3, Eq(nearest)) << *nearest3 << " " << *nearest;

    filtering_was_effective |=
        nearest != BruteForceNearestNeighbour(query_point, tree_points);
  }
  EXPECT_TRUE(filtering_was_effective) << "Filtering did nothing";
}

}  // namespace numerics
}  // namespace principia
