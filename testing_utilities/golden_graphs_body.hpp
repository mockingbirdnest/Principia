#pragma once

#include "testing_utilities/golden_graphs.hpp"

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <vector>
#include <span>
#include <string_view>

#include "absl/strings/ascii.h"
#include "gtest/gtest.h"
#include "lodepng/lodepng.h"
#include "numerics/fma.hpp"

namespace principia {
namespace testing_utilities {
namespace _golden_graphs {
namespace internal {

using namespace principia::numerics::_fma;

 std::vector<std::uint8_t> ReadFile(std::filesystem::path const& path) {
  std::vector<std::uint8_t> result;
  std::ifstream in(path, std::ios::binary | std::ios::in);
  while (in.good()) {
    char byte;
    in.read(&byte, 1);
    if (in.eof()) {
      break;
    }
    result.push_back(byte);
  }
}

template<typename Abscissa, typename Ordinate, typename Character>
void ExpectGoldenGraph(Graph<Abscissa, Ordinate> const& graph,
                       std::basic_string_view<Character> const suffix,
                       std::basic_string_view<Character> const test_file) {
  auto const image_path = std::filesystem::path(test_file)
                              .replace_extension()
                              .concat("_")
                              .concat(suffix)
                              .replace_extension(".png");
  auto const platform_image_path =
      OS_WIN && CanUseHardwareFMA
          ? image_path
          : std::filesystem::path(image_path)
                .replace_extension()
                .concat("_")
                .concat(absl::AsciiStrToLower(base::OperatingSystem))
                .concat(CanUseHardwareFMA ? "_fma" : "")
                .replace_extension(".png");
  auto const primary_golden = ReadFile(image_path);
  auto const platform_specific_golden = ReadFile(platform_image_path);
  std::uint8_t* actual_data;
  std::size_t actual_size;
  lodepng_encode32(&actual_data,
                   &actual_size,
                   reinterpret_cast<std::uint8_t const*>(graph.pixels().data()),
                   graph.width(),
                   graph.height());
  bool const matches_primary = std::equal(actual_data,
                                         actual_data + actual_size,
                                         primary_golden.begin(),
                                         primary_golden.end());
  if (matches_primary) {
    EXPECT_EQ(platform_specific_golden.size(), 0)
        << image_path << " matches, platform-specific override "
        << platform_image_path << " should be removed";
    std::filesystem::remove(platform_image_path);
  } else {
    bool const maches_platform_specific =
        std::equal(actual_data,
                   actual_data + actual_size,
                   platform_specific_golden.begin(),
                   platform_specific_golden.end());
    EXPECT_TRUE(maches_platform_specific)
        << platform_image_path
        << " has changed; golden size: " << platform_specific_golden.size()
        << " B, actual size: " << actual_size << " B";
    std::ofstream(platform_image_path, std::ios::binary | std::ios::out)
        .write(reinterpret_cast<char const*>(actual_data), actual_size);
    std::free(actual_data);
  }
}

}  // namespace internal
}  // namespace _golden_graphs
}  // namespace testing_utilities
}  // namespace principia
