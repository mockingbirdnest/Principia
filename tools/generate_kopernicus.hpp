#pragma once

#include <string>

namespace principia {
namespace tools {
namespace _generate_kopernicus {
namespace internal {

void GenerateKopernicusForSlippist1(
    std::string const& gravity_model_stem,
    std::string const& initial_state_stem);

}  // namespace internal

using internal::GenerateKopernicusForSlippist1;

}  // namespace _generate_kopernicus
}  // namespace tools
}  // namespace principia

namespace principia::tools {
using namespace principia::tools::_generate_kopernicus;
}  // namespace principia::tools
