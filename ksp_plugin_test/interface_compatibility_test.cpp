
#include "base/file.hpp"
#include "base/pull_serializer.hpp"
#include "base/push_deserializer.hpp"
#include "glog/logging.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "ksp_plugin/interface.hpp"
#include "ksp_plugin/plugin.hpp"
#include "serialization/ksp_plugin.pb.h"
#include "testing_utilities/serialization.hpp"

namespace principia {
namespace interface {

using base::OFStream;
using base::PullSerializer;
using base::PushDeserializer;
using ksp_plugin::Plugin;
using testing_utilities::ReadLinesFromBase64File;
using testing_utilities::ReadLinesFromHexadecimalFile;
using ::testing::NotNull;

const char preferred_compressor[] = "gipfeli";
const char preferred_encoder[] = "base64";

class InterfaceCompatibilityTest : public testing::Test {
 protected:

  void CheckSaveCompatibility(std::filesystem::path const& filename,
                              std::string_view const compressor,
                              std::string_view const encoder) {
    Plugin const* plugin = nullptr;

    // Read a plugin from a file containing only the "serialized_plugin = "
    // lines, with "serialized_plugin = " dropped.
    {
      PushDeserializer* deserializer = nullptr;
      auto const lines =
          encoder == "hexadecimal" ? ReadLinesFromHexadecimalFile(filename)
          : encoder == "base64"    ? ReadLinesFromBase64File(filename)
                                   : std::vector<std::string>{};
      CHECK(!lines.empty());
      LOG(ERROR) << "Deserialization starting";
      int i = 0;
      for (std::string const& line : lines) {
        principia__DeserializePlugin(line.c_str(),
                                     &deserializer,
                                     &plugin,
                                     compressor.data(),
                                     encoder.data());
      }
      principia__DeserializePlugin("",
                                   &deserializer,
                                   &plugin,
                                   compressor.data(),
                                   encoder.data());
      LOG(ERROR) << "Deserialization complete";
    }
    EXPECT_THAT(plugin, NotNull());

    // Write that plugin back to another file with the preferred format.
    {
      OFStream file(TEMP_DIR / "serialized_plugin.proto.b64");
      PullSerializer* serializer = nullptr;
      char const* b64 = nullptr;
      for (;;) {
        b64 = principia__SerializePlugin(plugin,
                                         &serializer,
                                         preferred_compressor,
                                         preferred_encoder);
        if (b64 == nullptr) {
          break;
        }
        file << b64 << "\n";
        principia__DeleteString(&b64);
      }
      LOG(ERROR) << "Serialization complete";
    }
    principia__DeletePlugin(&plugin);

    // Read the plugin from the new file to make sure that it's fine.
    {
      PushDeserializer* deserializer = nullptr;
      auto const lines =
          ReadLinesFromBase64File(TEMP_DIR / "serialized_plugin.proto.b64");
      for (std::string const& line : lines) {
        principia__DeserializePlugin(line.c_str(),
                                     &deserializer,
                                     &plugin,
                                     preferred_compressor,
                                     preferred_encoder);
      }
      principia__DeserializePlugin("",
                                   &deserializer,
                                   &plugin,
                                   preferred_compressor,
                                   preferred_encoder);
      LOG(ERROR) << "Deserialization complete";
    }

    EXPECT_THAT(plugin, NotNull());
    principia__DeletePlugin(&plugin);
  }
};

TEST_F(InterfaceCompatibilityTest, PreCohen) {
  CheckSaveCompatibility(
      SOLUTION_DIR / "ksp_plugin_test" / "pre_cohen.proto.hex",
      /*compressor=*/"",
      /*decoder=*/"hexadecimal");
}

// Use for debugging saves given by users.
TEST_F(InterfaceCompatibilityTest, DISABLED_SECULAR_Debug) {
  CheckSaveCompatibility(
      R"(P:\Public Mockingbird\Principia\Saves\2685\five-minute-scene-change-neptune.txt)",
      /*compressor=*/"gipfeli",
      /*decoder=*/"base64");
}

}  // namespace interface
}  // namespace principia
