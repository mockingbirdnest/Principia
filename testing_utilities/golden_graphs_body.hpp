#pragma once

#include "testing_utilities/golden_graphs.hpp"

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <vector>
#include <span>
#include <string_view>

#include "gtest/gtest.h"
#include "lodepng/lodepng.h"

namespace principia {
namespace testing_utilities {
namespace _golden_graphs {
namespace internal {

template<typename Abscissa, typename Ordinate, typename Character>
void ExpectGoldenGraph(Graph<Abscissa, Ordinate> const& graph,
                       std::basic_string_view<Character> const suffix,
                       std::basic_string_view<Character> const test_file) {
  auto const image_path = std::filesystem::path(test_file)
                              .replace_extension()
                              .concat("_")
                              .concat(suffix)
                              .replace_extension(".png");
  std::vector<std::uint8_t> golden;
  std::ifstream in(image_path, std::ios::binary | std::ios::in);
  while (in.good()) {
    char byte;
    in.read(&byte, 1);
    if (in.eof()) {
      break;
    }
    golden.push_back(byte);
  }
  std::uint8_t* actual_data;
  std::size_t actual_size;
  lodepng_encode32(&actual_data,
                   &actual_size,
                   reinterpret_cast<std::uint8_t const*>(graph.pixels().data()),
                   graph.width(),
                   graph.height());
  EXPECT_TRUE(std::equal(
      actual_data, actual_data + actual_size, golden.begin(), golden.end()))
      << image_path << " has changed; golden size: " << golden.size()
      << " B, actual size: " << actual_size << " B";
  std::ofstream(image_path, std::ios::binary | std::ios::out)
      .write(reinterpret_cast<char const*>(actual_data), actual_size);
  std::free(actual_data);
}

}  // namespace internal
}  // namespace _golden_graphs
}  // namespace testing_utilities
}  // namespace principia
