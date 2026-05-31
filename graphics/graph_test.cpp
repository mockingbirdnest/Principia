#include "graphics/graph.hpp"

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "quantities/numbers.hpp"
#include "numerics/elementary_functions.hpp"
#include "testing_utilities/golden_graphs.hpp"

namespace principia {
namespace graphics {
namespace _graph {

using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::geometry::_interval;
using namespace principia::numerics::_elementary_functions;
using ::testing::Eq;

class GraphTest : public ::testing::Test {};

TEST_F(GraphTest, Sinusoid) {
  Graph<Angle, double> graph(162,
                             100,
                             {-2 * π * Radian, 2 * π * Radian},
                             {-1.0, 1.0},
                             {.colour = {}, .alpha = 255});
  graph.PlotHorizontalLine(0, {255, 255, 255});
  graph.PlotVerticalLine(0 * Radian, {255, 255, 255});
  graph.PlotVerticalLine(π / 2 * Radian, {255, 255, 255}, {{0, 1}});
  graph.PlotVerticalLine(π * Radian, {255, 255, 255}, {{-1, 0}});
  graph.Plot(Sin,
             {-2 * π * Radian, 2 * π * Radian},
             {.red = 255, .green = 0, .blue = 0});
  graph.Plot(
      Cos, {-π * Radian, π * Radian}, {.red = 0, .green = 0, .blue = 255});
  EXPECT_GOLDEN_GRAPH(graph, "sinusoids");
}

}  // namespace _graph
}  // namespace graphics
}  // namespace principia