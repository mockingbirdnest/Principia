#pragma once

namespace principia {
namespace tools {
namespace _generate_profiles {
namespace internal {

void GenerateProfiles();

}  // namespace internal

using internal::GenerateProfiles;

}  // namespace _generate_profiles
}  // namespace tools
}  // namespace principia

namespace principia::tools {
using namespace principia::tools::_generate_profiles;
}  // namespace principia::tools
