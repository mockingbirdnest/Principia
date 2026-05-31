#pragma once

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <vector>

#include "graphics/graph.hpp"
#include "lodepng/lodepng.h"

#define EXPECT_GOLDEN_GRAPH(graph, suffix)                                   \
  [&] {                                                                      \
    auto const image_path = std::filesystem::path(__FILE__)                  \
                                .replace_extension()                         \
                                .concat("_" suffix)                          \
                                .replace_extension(".png");                  \
    std::vector<std::uint8_t> golden;                                        \
    std::ifstream in(image_path, std::ios::binary | std::ios::in);           \
    while (in.good()) {                                                      \
      in.read(reinterpret_cast<char*>(&golden.emplace_back()), 1);           \
      if (in.eof()) {                                                        \
        golden.pop_back();                                                   \
        break;                                                               \
      }                                                                      \
    }                                                                        \
    std::vector<std::uint8_t> actual;                                        \
    lodepng::encode(                                                         \
        actual,                                                              \
        reinterpret_cast<std::uint8_t const*>((graph).pixels().data()),      \
        (graph).width(),                                                     \
        (graph).height());                                                   \
    EXPECT_THAT(actual, ::testing::Eq(golden))                               \
        << image_path << " has changed";                                     \
    std::ofstream(image_path, std::ios::binary | std::ios::out)              \
        .write(reinterpret_cast<char const*>(actual.data()), actual.size()); \
  }()
