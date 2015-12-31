#include "tools/generate_profiles.hpp"

#include <experimental/filesystem>  // NOLINT
#include <fstream>

#include "glog/logging.h"
#include "google/protobuf/descriptor.h"
#include "serialization/journal.pb.h"
#include "tools/journal_proto_processor.hpp"

namespace principia {
namespace tools {

void GenerateProfiles() {
  JournalProtoProcessor processor;
  processor.ProcessMethods();

  // Now write the output.
  std::experimental::filesystem::path const directory =
      SOLUTION_DIR / "journal";
  std::ofstream profiles_hpp(directory / "profiles.gen.hpp");
  CHECK(profiles_hpp.good());
  for (auto const& cpp_method_type : processor.GetCppMethodTypes()) {
    profiles_hpp << cpp_method_type;
  }

  std::ofstream profiles_cpp(directory / "profiles.gen.cpp");
  CHECK(profiles_cpp.good());
  for (auto const& cpp_method_implementation :
          processor.GetCppMethodImplementations()) {
    profiles_cpp << cpp_method_implementation;
  }
}

}  // namespace tools
}  // namespace principia
