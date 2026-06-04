#include "graphics/graph.hpp"

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "quantities/numbers.hpp"
#include "numerics/elementary_functions.hpp"
#include "testing_utilities/golden_graphs.hpp"  // 🧙 for EXPECT_GOLDEN_GRAPH.
#include "graphics/colours.hpp"

namespace principia {
namespace graphics {
namespace _graph {

using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;
using namespace principia::geometry::_interval;
using namespace principia::numerics::_elementary_functions;
using namespace principia::graphics::_colours;

class GraphTest : public ::testing::Test {};

TEST_F(GraphTest, Sinusoid) {
  Graph<Angle, double> graph(162,
                             100,
                             {-2 * π * Radian, 2 * π * Radian},
                             {-1.0, 1.0},
                             /*background=*/Opaque(xkcd::black));
  graph.PlotHorizontalLine(0, xkcd::white);
  graph.PlotVerticalLine(0 * Radian, xkcd::white);
  graph.PlotVerticalLine(π / 2 * Radian, xkcd::white, {{0, 1}});
  graph.PlotVerticalLine(π * Radian, xkcd::white, {{-1, 0}});
  graph.Plot(Sin, {-2 * π * Radian, 2 * π * Radian}, xkcd::red);
  graph.Plot(Cos, {-π * Radian, π * Radian}, xkcd::blue);
  graph.ListPointPlot(
      std::views::iota(-24, 12) | std::views::transform([](double const x) {
        return std::pair{(x / 6) * Radian, Sin((x / 6 + 2) * Radian)};
      }),
      xkcd::green);
  EXPECT_GOLDEN_GRAPH(graph, "sinusoids");
}

}  // namespace _graph
}  // namespace graphics
}  // namespace principia