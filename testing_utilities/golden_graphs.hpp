#pragma once

#include <string_view>

#include "base/macros.hpp"  // 🧙 For FILESYSTEM_STRING.
#include "graphics/graph.hpp"

#define EXPECT_GOLDEN_GRAPH(graph, suffix)                            \
  (::principia::testing_utilities::_golden_graphs::ExpectGoldenGraph( \
      (graph),                                                        \
      FILESYSTEM_STRING_VIEW(suffix),                                 \
      FILESYSTEM_STRING_VIEW(__FILE__)))

namespace principia {
namespace testing_utilities {
namespace _golden_graphs {
namespace internal {

using namespace principia::graphics::_graph;

// This is templatized on the character type because it needs to
// platform-appropriate for paths (UTF-16 on Windows, UTF-8 on *nix).  We take
// pointers rather than `basic_string_view<Character>` so that we can deduce
// `Character` from the literals.
template<typename Abscissa, typename Ordinate, typename Character>
void ExpectGoldenGraph(Graph<Abscissa, Ordinate> const& graph,
                       std::basic_string_view<Character> const suffix,
                       std::basic_string_view<Character> const test_file);

}  // namespace internal

using internal::ExpectGoldenGraph;

}  // namespace _golden_graphs
}  // namespace testing_utilities
}  // namespace principia

#include "testing_utilities/golden_graphs_body.hpp"
