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
  std::ofstream profiles_hpp(directory / "profiles.hpp");
  CHECK(profiles_hpp.good());

  profiles_hpp << "#pragma once\n";
  profiles_hpp << "\n";
  profiles_hpp << "#include \"base/not_null.hpp\"\n";
  profiles_hpp << "#include \"journal/player.hpp\"\n";
  profiles_hpp << "#include \"ksp_plugin/interface.hpp\"\n";
  profiles_hpp << "#include \"serialization/journal.pb.h\"\n";
  profiles_hpp << "\n";
  profiles_hpp << "namespace principia {\n";
  profiles_hpp << "\n";
  profiles_hpp << "using base::not_null;\n";
  profiles_hpp << "using ksp_plugin::KSPPart;\n";
  profiles_hpp << "using ksp_plugin::LineAndIterator;\n";
  profiles_hpp << "using ksp_plugin::NavigationFrame;\n";
  profiles_hpp << "using ksp_plugin::Plugin;\n";
  profiles_hpp << "using ksp_plugin::QP;\n";
  profiles_hpp << "using ksp_plugin::WXYZ;\n";
  profiles_hpp << "using ksp_plugin::XYZ;\n";
  profiles_hpp << "using ksp_plugin::XYZSegment;\n";
  profiles_hpp << "\n";
  profiles_hpp << "namespace journal {\n";
  profiles_hpp << "\n";

  for (auto const& cpp_method_type : processor.GetCppMethodTypes()) {
    profiles_hpp << cpp_method_type;
  }

  profiles_hpp << "}  // namespace journal\n";
  profiles_hpp << "}  // namespace principia\n";
}

}  // namespace tools
}  // namespace principia
