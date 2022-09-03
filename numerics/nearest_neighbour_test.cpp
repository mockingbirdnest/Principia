#include "numerics/nearest_neighbour.hpp"

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace principia {
namespace numerics {

using geometry::Frame;
using geometry::Vector;
using ::testing::Optional;

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
  PrincipalComponentPartitioningTree<V> tree({v1, v2, v3, v4},
                                             /*max_values_per_cell=*/1);

  EXPECT_THAT(tree.FindNearestNeighbour(V({-0.5, -0.5, 0})), Optional(v1));
  EXPECT_THAT(tree.FindNearestNeighbour(V({-0.5, 0.5, 0})), Optional(v2));
  EXPECT_THAT(tree.FindNearestNeighbour(V({1, 1, 0})), Optional(v4));
  EXPECT_THAT(tree.FindNearestNeighbour(V({3, 3, 0})), Optional(v4));
}

}  // namespace numerics
}  // namespace principia
