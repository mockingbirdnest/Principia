#pragma once

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <vector>
#include <string_view>

#include "base/macros.hpp"  // 🧙 For FILESYSTEM_STRING.
#include "gmock/gmock.h"
#include "graphics/graph.hpp"
#include "gtest/gtest.h"
#include "lodepng/lodepng.h"

#define EXPECT_GOLDEN_GRAPH(graph, suffix)                            \
  (::principia::testing_utilities::_golden_graphs::ExpectGoldenGraph( \
      (graph), FILESYSTEM_STRING(suffix), FILESYSTEM_STRING(__FILE__)))
namespace principia {
namespace testing_utilities {
namespace _golden_graphs {
namespace internal {

using namespace principia::graphics::_graph;
using ::testing::Eq;

template<typename Abscissa, typename Ordinate, typename Character>
void ExpectGoldenGraph(Graph<Abscissa, Ordinate> const& graph,
                       Character const* const suffix,
                       Character const* const test_file) {
  auto const image_path = std::filesystem::path(test_file)
                              .replace_extension()
                              .concat("_")
                              .concat(suffix)
                              .replace_extension(".png");
  std::vector<std::uint8_t> golden;
  std::ifstream in(image_path, std::ios::binary | std::ios::in);
  while (in.good()) {
    in.read(reinterpret_cast<char*>(&golden.emplace_back()), 1);
    if (in.eof()) {
      golden.pop_back();
      break;
    }
  }
  std::vector<std::uint8_t> actual;
  lodepng::encode(actual,
                  reinterpret_cast<std::uint8_t const*>(graph.pixels().data()),
                  graph.width(),
                  graph.height());
  EXPECT_THAT(actual, Eq(golden)) << image_path << " has changed";
  std::ofstream(image_path, std::ios::binary | std::ios::out)
      .write(reinterpret_cast<char const*>(actual.data()), actual.size());
}

}  // namespace internal

using internal::ExpectGoldenGraph;

}  // namespace _golden_graphs
}  // namespace testing_utilities
}  // namespace principia
