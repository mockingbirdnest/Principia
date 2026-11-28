#include "graphics/graph.hpp"

#include <algorithm>

#include "geometry/interval.hpp"
#include "gmock/gmock.h"
#include "graphics/colours.hpp"
#include "gtest/gtest.h"
#include "numerics/elementary_functions.hpp"
#include "quantities/quantities.hpp"
#include "quantities/si.hpp"
#include "testing_utilities/golden_graphs.hpp"  // 🧙 for EXPECT_GOLDEN_GRAPH.

namespace principia {
namespace graphics {
namespace _graph {

using namespace principia::base::_not_null;
using namespace principia::geometry::_interval;
using namespace principia::graphics::_colours;
using namespace principia::numerics::_elementary_functions;
using namespace principia::quantities::_quantities;
using namespace principia::quantities::_si;

class GraphTest : public ::testing::Test {};


TEST_F(GraphTest, Frame) {
  // This test draws lines aligned with the centres of the edge pixels, using a
  // square pixel-sized brush: the resulting image does not change with
  // anti-aliasing.
  for (std::int64_t supersampling = 1; supersampling <= 8; ++supersampling) {
    Graph<Angle, double> graph(
        std::round(100 * φ),
        100,
        {-2 * π * Radian, 2 * π * Radian},
        {-1.0, 1.0},
        /*background=*/Opaque(xkcd::red),
        {.anti_aliasing = {.supersampling = supersampling}});
    graph.PlotHorizontalLine(-1.0, xkcd::blue);
    graph.PlotHorizontalLine(1.0, xkcd::blue);
    graph.PlotVerticalLine(-2 * π * Radian, xkcd::blue);
    graph.PlotVerticalLine(2 * π * Radian, xkcd::blue);
    graph.AntiAlias();
    EXPECT_GOLDEN_GRAPH(graph, "frame");
  }
}

TEST_F(GraphTest, Sinusoid) {
  Graph<Angle, double> aliased_graph(std::round(100 * φ),
                             100,
                             {-2 * π * Radian, 2 * π * Radian},
                             {-1.0, 1.0},
                             /*background=*/Opaque(xkcd::black));
  Graph<Angle, double> anti_aliased_graph(std::round(100 * φ),
                             100,
                             {-2 * π * Radian, 2 * π * Radian},
                             {-1.0, 1.0},
                             /*background=*/Opaque(xkcd::black),
                             {.anti_aliasing={.supersampling = 8}});
  for (not_null<Graph<Angle, double>*> const graph :
       {&aliased_graph, &anti_aliased_graph}) {
    graph->PlotHorizontalLine(0, xkcd::white);
    graph->PlotVerticalLine(0 * Radian, xkcd::white);
    graph->PlotVerticalLine(π / 2 * Radian, xkcd::white, {{0, 1}});
    graph->PlotVerticalLine(π * Radian, xkcd::white, {{-1, 0}});
    graph->Plot(Sin, {-2 * π * Radian, 2 * π * Radian}, xkcd::red);
    graph->Plot(Cos, {-π * Radian, π * Radian}, xkcd::blue);
    graph->ListPointPlot(
        std::views::iota(-24, 12) | std::views::transform([](double const x) {
          return std::pair{(x / 6) * Radian, Sin((x / 6 + 2) * Radian)};
        }),
        xkcd::green);
    graph->ListPointPlot(
        std::views::iota(-24, 12) | std::views::transform([](double const x) {
          return std::pair{(x / 6) * Radian, -Sin((x / 6 + 2) * Radian)};
        }),
        xkcd::green,
        {.joined = true});
    graph->Plot([](Angle const α) { return -Sin(α + 1 * Radian); },
                {-3 * Radian, 0 * Radian},
                xkcd::green);
    graph->AntiAlias();
  }
  EXPECT_GOLDEN_GRAPH(aliased_graph, "pixelated_sinusoids");
  EXPECT_GOLDEN_GRAPH(anti_aliased_graph, "sinusoids");
}

}  // namespace _graph
}  // namespace graphics
}  // namespace principia
