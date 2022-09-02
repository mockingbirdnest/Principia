#include "numerics/nearest_neighbour.hpp"

#include "geometry/frame.hpp"
#include "geometry/grassmann.hpp"
#include "gtest/gtest.h"

namespace principia {
namespace numerics {

using geometry::Frame;
using geometry::Vector;

class NearestNeighbourTest : public ::testing::Test {
protected:
  using World = Frame<enum class WorldTag>;

  NearestNeighbourTest() : tree_({}, /*max_values_per_cell=*/3) {}

  PrincipalComponentPartitioningTree<Vector<double, World>> tree_;
};

}  // namespace numerics
}  // namespace principia
