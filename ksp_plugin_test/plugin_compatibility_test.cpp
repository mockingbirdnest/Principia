
#include <memory>
#include <string>
#include <vector>

#include "base/file.hpp"
#include "base/not_null.hpp"
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

using base::not_null;
using base::OFStream;
using base::PullSerializer;
using base::PushDeserializer;
using ksp_plugin::Plugin;
using testing_utilities::ReadLinesFromBase64File;
using testing_utilities::ReadLinesFromHexadecimalFile;
using ::testing::NotNull;

const char preferred_compressor[] = "gipfeli";
const char preferred_encoder[] = "base64";

class PluginCompatibilityTest : public testing::Test {
 protected:
  PluginCompatibilityTest() {
    google::SetStderrLogging(google::WARNING);
  }

  // Reads a plugin from a file containing only the "serialized_plugin = "
  // lines, with "serialized_plugin = " dropped.
  not_null<std::unique_ptr<Plugin const>> ReadPluginFromFile(
      std::filesystem::path const& filename,
      std::string_view const compressor,
      std::string_view const encoder) {
    Plugin const* plugin = nullptr;

    PushDeserializer* deserializer = nullptr;
    auto const lines =
        encoder == "hexadecimal" ? ReadLinesFromHexadecimalFile(filename)
        : encoder == "base64"    ? ReadLinesFromBase64File(filename)
                                  : std::vector<std::string>{};
    CHECK(!lines.empty());

    LOG(ERROR) << "Deserialization starting";
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

    return std::unique_ptr<Plugin const>(plugin);
  }

  // Writes a plugin to a file.
  void WritePluginToFile(std::filesystem::path const& filename,
                         std::string_view const compressor,
                         std::string_view const encoder,
                         not_null<std::unique_ptr<Plugin const>> plugin) {
    OFStream file(filename);
    PullSerializer* serializer = nullptr;
    char const* b64 = nullptr;

    LOG(ERROR) << "Serialization starting";
    for (;;) {
      b64 = principia__SerializePlugin(plugin.get(),
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

    Plugin const* released_plugin = plugin.release();
    principia__DeletePlugin(&released_plugin);
  }

  void CheckSaveCompatibility(std::filesystem::path const& filename,
                              std::string_view const compressor,
                              std::string_view const encoder) {
    // Read a plugin from the given file.
    auto plugin1 = ReadPluginFromFile(filename, compressor, encoder);

    // Write that plugin back to another file with the preferred format.
    WritePluginToFile(TEMP_DIR / "serialized_plugin.proto.b64",
                      preferred_compressor,
                      preferred_encoder,
                      std::move(plugin1));

    // Read the plugin from the new file to make sure that it's fine.
    auto plugin2 = ReadPluginFromFile(TEMP_DIR / "serialized_plugin.proto.b64",
                                      preferred_compressor,
                                      preferred_encoder);
  }
};

TEST_F(PluginCompatibilityTest, PreCartan) {
  // This space for rent.
}

TEST_F(PluginCompatibilityTest, PreCohen) {
  CheckSaveCompatibility(
      SOLUTION_DIR / "ksp_plugin_test" / "pre_cohen.proto.hex",
      /*compressor=*/"",
      /*decoder=*/"hexadecimal");
}

// Use for debugging saves given by users.
TEST_F(PluginCompatibilityTest, DISABLED_SECULAR_Debug) {
  CheckSaveCompatibility(
      R"(P:\Public Mockingbird\Principia\Saves\2685\five-minute-scene-change-neptune.txt)",
      /*compressor=*/"gipfeli",
      /*decoder=*/"base64");
}

}  // namespace interface
}  // namespace principia
